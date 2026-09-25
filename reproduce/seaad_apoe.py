"""APOE models of donor-mean GEM expression in SEA-AD vascular subsets."""
from pathlib import Path
import sys,json,re
R=Path(__file__).resolve().parents[1];
(R/'results').mkdir(exist_ok=True)
import numpy as np,pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
spec={'outcome':'Donor mean normalized GEM expression from the SEA-AD analysis','contrast':'Public dementia versus normal labels; not AD versus control','base_covariates':['age','sex','PMI','RIN'],'genotype':'Official full APOE genotype; epsilon2 and epsilon4 dosage','complete_case':'Identical donor rows for baseline, dosage, carrier and interaction models in each region/cell subset; no imputation','inference':'OLS with HC3 heteroskedasticity-robust SE and t-based 95% CI; BH across four subsets separately for each model','strata':'Fit each epsilon4 stratum only if >=5 donors in both diagnostic groups and full-rank design','overlap':'Report pairwise donor intersections; regions are not independent cohorts','recorded_before_model_evaluation':'2026-09-10'}
(R/'results/SEAAD_APOE_frozen_spec.json').write_text(json.dumps(spec,indent=2))
meta=pd.read_excel(R/'data/SEAAD_official_donor_metadata.xlsx').rename(columns={'Donor ID':'donor_id','Age at Death':'age','Sex':'metadata_sex','APOE Genotype':'genotype','Cognitive Status':'cognitive_status'})
meta=meta[['donor_id','age','metadata_sex','genotype','cognitive_status','PMI','RIN']]
assert meta.donor_id.is_unique
def dosage(g,a):
 if not isinstance(g,str) or not re.fullmatch(r'[234]/[234]',g.strip()):return np.nan
 return g.split('/').count(str(a))
meta['e2']=meta.genotype.apply(lambda g:dosage(g,2));meta['e4']=meta.genotype.apply(lambda g:dosage(g,4))
sea=pd.read_csv(R/'data/SEAAD_GEM_donor_values.csv');d=sea.merge(meta,on='donor_id',how='left',validate='m:1',indicator=True)
d['diagnosis']=d.disease.eq('dementia').astype(int);d['carrier']=(d.e4>0).astype(float);d.loc[d.e4.isna(),'carrier']=np.nan
d['sex_match']=d.sex.str.lower().eq(d.metadata_sex.str.lower())
d['diagnosis_match']=d.cognitive_status.map({'Dementia':'dementia','No dementia':'normal'}).eq(d.disease)
d.to_csv(R/'results/SEAAD_APOE_metadata_matching.csv',index=False)
matched=d[d._merge.eq('both')]
print('sex mismatch',matched.loc[~matched.sex_match,['donor_id','sex','metadata_sex']].drop_duplicates().to_dict('records'))
print('diagnosis labels',matched[['disease','cognitive_status']].value_counts().to_string())
assert matched.sex_match.all() and matched.diagnosis_match.all()
cc=d.dropna(subset=['expr','age','sex','PMI','RIN','e2','e4']).copy();cc['age_c']=cc.age-cc.age.mean();cc['PMI_c']=cc.PMI-cc.PMI.mean();cc['RIN_c']=cc.RIN-cc.RIN.mean()
cc.to_csv(R/'results/SEAAD_APOE_complete_case_values.csv',index=False)
ct=cc.groupby(['dataset','genotype','disease']).size().rename('n').reset_index();ct.to_csv(R/'results/SEAAD_APOE_genotype_crosstab.csv',index=False)
cc.groupby(['dataset','carrier','disease']).size().rename('n').reset_index().to_csv(R/'results/SEAAD_APOE4_crosstab.csv',index=False)
sets={k:set(x.donor_id) for k,x in cc.groupby('dataset')};overlap=pd.DataFrame([{'subset1':a,'subset2':b,'shared_donors':len(sets[a]&sets[b]),'n1':len(sets[a]),'n2':len(sets[b])} for a in sets for b in sets]);overlap.to_csv(R/'results/SEAAD_APOE_donor_overlap.csv',index=False)
base='age_c + C(sex) + PMI_c + RIN_c';forms={'diagnosis_only':'expr ~ diagnosis','base':'expr ~ '+base+' + diagnosis','APOE_dosage':'expr ~ '+base+' + e2 + e4 + diagnosis','APOE4_carrier':'expr ~ '+base+' + carrier + diagnosis','interaction':'expr ~ '+base+' + diagnosis * carrier','noncarrier':'expr ~ '+base+' + diagnosis','carrier':'expr ~ '+base+' + diagnosis'}
rows=[];infl=[]
for dataset,z in cc.groupby('dataset'):
 for model,formula in forms.items():
  zz=z if model not in ['noncarrier','carrier'] else z[z.carrier.eq(int(model=='carrier'))]
  if zz.groupby('disease').size().min()<5 or zz.disease.nunique()<2:
   rows.append(dict(dataset=dataset,model=model,status='Stratum size below prespecified minimum'));continue
  fit=smf.ols(formula,zz).fit(cov_type='HC3',use_t=True);rank=np.linalg.matrix_rank(fit.model.exog)
  if rank<fit.model.exog.shape[1]:
   rows.append(dict(dataset=dataset,model=model,status='Rank-deficient design'));continue
  term='diagnosis:carrier' if model=='interaction' else 'diagnosis';ci=fit.conf_int().loc[term]
  rows.append(dict(dataset=dataset,model=model,status='estimated',term=term,n_dementia=int(zz.diagnosis.sum()),n_normal=int((zz.diagnosis==0).sum()),coefficient=fit.params[term],SE_HC3=fit.bse[term],CI_low=ci.iloc[0],CI_high=ci.iloc[1],PValue=fit.pvalues[term],design_rank=rank,n_columns=fit.model.exog.shape[1]))
  pd.DataFrame(fit.model.exog,columns=fit.model.exog_names,index=zz.donor_id).to_csv(R/'results'/f'SEAAD_design_{dataset}_{model}.csv',index_label='donor_id')
  if model=='APOE_dosage':
   ols=smf.ols(formula,zz).fit();v=ols.get_influence();infl.extend(dict(dataset=dataset,donor_id=donor,Cooks_distance=cook,hat=hat) for donor,cook,hat in zip(zz.donor_id,v.cooks_distance[0],v.hat_matrix_diag))
out=pd.DataFrame(rows);out['q_BH_4_subsets']=out.groupby('model').PValue.transform(lambda s:multipletests(s,method='fdr_bh')[1])
out.to_csv(R/'results/SEAAD_APOE_nested_models.csv',index=False);pd.DataFrame(infl).to_csv(R/'results/SEAAD_APOE_influence.csv',index=False)
print(out[out.model.isin(['base','APOE_dosage','interaction'])].to_string(index=False));print(cc.groupby(['dataset','disease']).size().to_string());print('Unique complete-case donors',cc.donor_id.nunique())
