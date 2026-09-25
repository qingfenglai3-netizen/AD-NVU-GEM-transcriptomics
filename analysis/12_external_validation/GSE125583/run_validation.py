from pathlib import Path
import csv,gzip,io,json,tarfile,hashlib,platform
import numpy as np,pandas as pd,sklearn
from scipy import stats
from scipy.special import expit
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score,roc_curve
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

W=Path(__file__).resolve().parent
P=W.parents[2]
import argparse
parser=argparse.ArgumentParser()
parser.add_argument('--input-dir',type=Path,default=W)
parser.add_argument('--output-dir',type=Path,default=W/'recomputed')
args=parser.parse_args()
I=args.input_dir
O=args.output_dir;O.mkdir(parents=True,exist_ok=True)
assert sklearn.__version__=='1.5.2'
spec=json.loads((W/'analysis_spec_frozen.json').read_text())
mapping=json.loads((W/'gene_mapping.json').read_text())
train=pd.read_csv(P/'data/discovery_GSE132903_labelblind_mean_features.csv')
genes=train.columns[2:].tolist(); assert len(genes)==18 and 'GEM' not in genes
rows=[r for r in csv.reader(io.StringIO(gzip.decompress((W/'GSE125583_series_matrix.txt.gz').read_bytes()).decode()),delimiter='\t') if r]
ids=next(r[1:] for r in rows if r[0]=='!Sample_geo_accession')
meta=pd.DataFrame({'sample_id':ids})
for r in rows:
    if r[0]=='!Sample_characteristics_ch1':
        key=r[1].split(': ',1)[0]
        meta[key]=[v.split(': ',1)[1] for v in r[1:]]
    elif r[0]=='!Sample_description':meta['repeated_GSE95587']=[v.split(': ',1)[1] for v in r[1:]]
assert len(meta)==289 and meta['case id'].is_unique and meta.sample_id.is_unique
meta['group']=meta.diagnosis.map({"Alzheimer's disease":'AD','control':'CN'});assert meta.group.notna().all()
meta['age']=pd.to_numeric(meta.age,errors='coerce')
meta['sex']=meta['Sex'].replace('NA',np.nan)
values={};counts={};input_hashes={}
with tarfile.open(I/'GSE125583_RAW.tar') as tar:
    for member in tar.getmembers():
        if not member.isfile():continue
        gsm=member.name.split('_')[0];b=tar.extractfile(member).read()
        input_hashes[gsm]=hashlib.sha256(b).hexdigest()
        df=pd.read_csv(io.BytesIO(gzip.decompress(b)),sep='\t',comment='#',dtype={'ID_REF':str})
        assert df.ID_REF.is_unique
        df=df.set_index('ID_REF');sel=df.loc[list(mapping.values()),['VALUE','count']]
        assert np.isfinite(sel.to_numpy(float)).all() and (sel.to_numpy(float)>=0).all()
        values[gsm]={g:float(df.loc[e,'VALUE']) for g,e in mapping.items()}
        counts[gsm]={g:float(df.loc[e,'count']) for g,e in mapping.items()}
assert set(values)==set(ids)
expr=pd.DataFrame.from_dict(values,orient='index').loc[ids];cnt=pd.DataFrame.from_dict(counts,orient='index').loc[ids]
expr.index.name=cnt.index.name='sample_id';expr.to_csv(O/'module19_nRPKM.csv');cnt.to_csv(O/'module19_counts.csv')
meta.to_csv(O/'donor_metadata.csv',index=False)
pd.DataFrame({'symbol':list(mapping),'entrez_id':list(mapping.values()),'nonzero_nRPKM':(expr>0).sum().values,'nonzero_counts':(cnt>0).sum().values}).to_csv(O/'gene_coverage.csv',index=False)
x=train[genes].to_numpy(float);y=train.group.eq('AD').to_numpy(int)
scaler=StandardScaler().fit(x)
model=LogisticRegression(penalty='l2',C=1,solver='liblinear',max_iter=2000,random_state=20260608).fit(scaler.transform(x),y)
# Reproduce the archived GSE118553 predictions with the same fixed model.
bench=pd.read_csv(P/'data/GSE118553_TC_donor_features.csv')
bp=model.predict_proba(scaler.transform(bench[genes].to_numpy(float)))[:,1]
bench_auc=roc_auc_score(bench.group.eq('AD'),bp)
assert abs(bench_auc-0.8194285714285714)<0.002,bench_auc
arch=pd.read_csv(P/'data/GSE118553_TC_predictions.csv');arch=arch[(arch.model=='Logistic_L2')&(arch.selector=='all')]
aligned=arch.set_index('sample_id').loc[bench.sample_id,'prediction'].to_numpy()
assert np.allclose(bp,aligned,rtol=0,atol=1e-12)
print('Existing GSE118553 probabilities reproduced; AUC',bench_auc,flush=True)
logexpr=np.log2(expr+1);z=logexpr[genes].to_numpy()
strict=model.decision_function(scaler.transform(z))
adapted=model.decision_function(StandardScaler().fit_transform(z))
pred=meta[['sample_id','case id','group','repeated_GSE95587']].copy()
pred['strict_decision_score']=strict;pred['strict_probability']=expit(strict)
pred['adapted_decision_score']=adapted;pred['adapted_probability']=expit(adapted)
pred.to_csv(O/'external_predictions.csv',index=False)
pd.DataFrame({'gene':genes,'coefficient':model.coef_[0],'discovery_mean':scaler.mean_,'discovery_sd':scaler.scale_,'RNAseq_log_mean':z.mean(axis=0),'RNAseq_log_sd':z.std(axis=0)}).to_csv(O/'fixed_model_and_scale_shift.csv',index=False)
yy=meta.group.eq('AD').to_numpy(int);populations={'all_289':np.ones(len(meta),bool),'exclude_GSE95587_overlap':meta.repeated_GSE95587.eq('NA').to_numpy()}
aucrows=[]
for pop,mask in populations.items():
    for name,score in [('strict_frozen_preprocessing',strict),('label_blind_cohort_z_adapted',adapted)]:
        a=score[mask];b=yy[mask];rng=np.random.default_rng(20260925);p=np.where(b==1)[0];n=np.where(b==0)[0]
        boots=[]
        for _ in range(2000):
            ix=np.r_[rng.choice(p,len(p)),rng.choice(n,len(n))];boots.append(roc_auc_score(b[ix],a[ix]))
        ci=np.quantile(boots,[.025,.975]);auc=roc_auc_score(b,a)
        aucrows.append(dict(population=pop,method=name,n_AD=int(sum(b)),n_CN=int(sum(b==0)),auc=auc,ci_low=ci[0],ci_high=ci[1],unique_scores=len(np.unique(a)),probability_auc=roc_auc_score(b,expit(a))))
aucdf=pd.DataFrame(aucrows);aucdf.to_csv(O/'AUC_results.csv',index=False)
gemrows=[];gem=logexpr.GEM.to_numpy();frame=meta.copy();frame['GEM_log2_nRPKM_plus1']=gem;frame['AD']=yy;frame['male']=frame.sex.map({'M':1,'F':0})
frame.to_csv(O/'GEM_donor_values.csv',index=False)
for pop,mask in populations.items():
    a=gem[mask&(yy==1)];b=gem[mask&(yy==0)];diff=a.mean()-b.mean();se=np.sqrt(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))
    df=(a.var(ddof=1)/len(a)+b.var(ddof=1)/len(b))**2/((a.var(ddof=1)/len(a))**2/(len(a)-1)+(b.var(ddof=1)/len(b))**2/(len(b)-1))
    ci=diff+np.array([-1,1])*stats.t.ppf(.975,df)*se
    gemrows.append(dict(population=pop,method='unadjusted_Welch',n_AD=len(a),n_CN=len(b),effect=diff,ci_low=ci[0],ci_high=ci[1],p=stats.ttest_ind(a,b,equal_var=False).pvalue,mannwhitney_p=stats.mannwhitneyu(a,b,alternative='two-sided').pvalue))
    sub=frame.loc[mask].dropna(subset=['age','male']);fit=sm.OLS(sub.GEM_log2_nRPKM_plus1,sm.add_constant(sub[['AD','age','male']])).fit(cov_type='HC3')
    ci=fit.conf_int().loc['AD'];gemrows.append(dict(population=pop,method='age_sex_adjusted_OLS_HC3',n_AD=int(sub.AD.sum()),n_CN=int(sum(sub.AD==0)),effect=fit.params['AD'],ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues['AD']))
gemdf=pd.DataFrame(gemrows);gemdf.to_csv(O/'GEM_results.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(10,4.2))
for name,score,color in [('Frozen preprocessing',strict,'#254b80'),('Cohort z-score adaptation',adapted,'#a54936')]:
    fpr,tpr,_=roc_curve(yy,score);axs[0].plot(fpr,tpr,label=f'{name}: AUC {roc_auc_score(yy,score):.3f}',color=color)
axs[0].plot([0,1],[0,1],'--',color='grey',lw=.8);axs[0].set(xlabel='False positive rate',ylabel='True positive rate',title='Fixed 18-gene model: GSE125583');axs[0].legend(fontsize=8,loc='lower right')
rng=np.random.default_rng(20260925)
for k,(group,color) in enumerate([('CN','#254b80'),('AD','#a54936')]):
    v=gem[meta.group.eq(group)];axs[1].scatter(k+rng.uniform(-.16,.16,len(v)),v,s=10,alpha=.45,color=color)
    axs[1].plot([k-.2,k+.2],[v.mean(),v.mean()],color='black',lw=2)
axs[1].set_xticks([0,1],['Control (70)','AD (219)']);axs[1].set(ylabel='GEM log2(nRPKM + 1)',title='Bulk fusiform gyrus GEM expression')
for ax in axs:ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(O/'validation_summary.png',dpi=200);fig.savefig(O/'validation_summary.pdf');plt.close(fig)
qc={'python':platform.python_version(),'sklearn':sklearn.__version__,'input_tar_sha256':hashlib.sha256((I/'GSE125583_RAW.tar').read_bytes()).hexdigest(),'frozen_spec_sha256':hashlib.sha256((W/'analysis_spec_frozen.json').read_bytes()).hexdigest(),'all19genes_all289samples_complete':True,'n_unique_donors':len(meta),'existing_benchmark_max_probability_difference':float(np.max(abs(bp-aligned))),'GSE118553_auc':bench_auc,'analysis_scope':'Additional exploratory cross-platform and brain-region transfer; donor overlap with discovery not established from sample IDs alone'}
(O/'verification.json').write_text(json.dumps(qc,indent=2));(O/'sample_file_sha256.json').write_text(json.dumps(input_hashes,indent=2))
print(aucdf.to_string(index=False),flush=True);print(gemdf.to_string(index=False),flush=True)
