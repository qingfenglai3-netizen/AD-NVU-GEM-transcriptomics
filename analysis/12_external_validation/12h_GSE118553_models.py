"""Donor-level evaluation of the fixed pipelines in GSE118553 temporal cortex."""
from pathlib import Path
import json,gzip,csv,hashlib,importlib.util,warnings
import numpy as np,pandas as pd,sklearn
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.metrics import roc_auc_score
assert sklearn.__version__=='1.5.2'
R=Path(__file__).resolve().parents[2];O=R/'intermediate/GSE118553';spec=json.loads((O/'frozen_spec.json').read_text(encoding='utf-8'));assert spec['region']=='Temporal_Cortex only'
tr=pd.read_csv(R/'results/discovery_GSE132903_labelblind_mean_features.csv');genes=tr.columns[2:].tolist()
meta=pd.read_csv(R/'results/GSE118553_TC_metadata_all.csv');meta=meta[meta['disease state'].isin(['AD','control'])].copy();meta['group']=meta['disease state'].map({'AD':'AD','control':'CN'})
assert meta.groupby('donor_id').group.nunique().max()==1
mapping=pd.read_csv(R/'results/ML_probe_lineage.csv');assert set(mapping.symbol)==set(genes) and mapping.probe_id.is_unique
probe_set=set(mapping.probe_id);data=[]
with gzip.open(O/'GSE118553_series_matrix.txt.gz','rt',encoding='utf-8') as fh:
 for line in fh:
  if line.startswith('!series_matrix_table_begin'):break
 reader=csv.reader(fh,delimiter='\t');header=next(reader);idx=[header.index(s) for s in meta.sample_id]
 for row in reader:
  if row[0].startswith('!series_matrix_table_end'):break
  if row[0] in probe_set:data.append([row[0]]+[float(row[i]) for i in idx])
ex=pd.DataFrame(data,columns=['probe']+meta.sample_id.tolist()).set_index('probe');assert set(ex.index)==probe_set
X=pd.DataFrame({g:ex.loc[mapping.loc[mapping.symbol.eq(g),'probe_id']].mean(axis=0) for g in genes})
sample=meta.join(X,on='sample_id');sample.to_csv(R/'results/GSE118553_TC_sample_features.csv',index=False)
d=sample.groupby(['donor_id','group'],sort=True)[genes].mean().reset_index();d['sample_id']='GSE118553_'+d.donor_id.str.replace(' ','_',regex=False);d['n_arrays']=d.donor_id.map(meta.donor_id.value_counts())
d.to_csv(R/'results/GSE118553_TC_donor_features.csv',index=False)
assert d.donor_id.is_unique and np.isfinite(d[genes]).all().all() and d.group.value_counts().min()>=10
src=Path(r'<WORKDIR>/AD-NVU-GEM-transcriptomics/analysis/11_translational/11a_bulk_multimodel_cross_validation.py');sp=importlib.util.spec_from_file_location('original',src);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
x=tr[genes].to_numpy(float);y=tr.group.eq('AD').to_numpy(int);xt=d[genes].to_numpy(float);yy=d.group.eq('AD').to_numpy(int);rng=np.random.default_rng(20260910);rows=[];prs=[]
for name,est in m.model_zoo().items():
 for sel in ['all','top3','top5','top8','top12','top16']:
  steps=[('scale',StandardScaler())]
  if sel!='all':steps.append(('select',SelectKBest(f_classif,k=int(sel[3:]))))
  steps.append(('model',clone(est)));pipe=Pipeline(steps).fit(x,y);pred=pipe.predict_proba(xt)[:,1]
  assert np.array_equal(pipe.classes_,[0,1]);ci=[np.nan,np.nan]
  if name=='Logistic_L2' and sel=='all':
   boot=[]
   for _ in range(2000):
    ix=np.r_[rng.choice(np.where(yy==0)[0],sum(yy==0)),rng.choice(np.where(yy==1)[0],sum(yy==1))];boot.append(roc_auc_score(yy[ix],pred[ix]))
   ci=np.quantile(boot,[.025,.975])
  rows.append(dict(model=name,selector=sel,region='TC',cohort='GSE118553',n_AD=int(sum(yy)),n_CN=int(sum(yy==0)),auc=roc_auc_score(yy,pred),ci_low=ci[0],ci_high=ci[1]))
  prs.extend(dict(model=name,selector=sel,region='TC',cohort='GSE118553',sample_id=s,group=g,prediction=float(p)) for s,g,p in zip(d.sample_id,d.group,pred))
 print(name,flush=True)
pd.DataFrame(rows).to_csv(R/'results/GSE118553_TC_model_performance.csv',index=False);pd.DataFrame(prs).to_csv(R/'results/GSE118553_TC_predictions.csv',index=False)
(O/'execution.json').write_text(json.dumps({'sklearn':sklearn.__version__,'training_sha256':hashlib.sha256((R/'results/discovery_GSE132903_labelblind_mean_features.csv').read_bytes()).hexdigest(),'n_arrays':len(meta),'n_donors':len(d),'n_AD':int(sum(yy)),'n_CN':int(sum(yy==0)),'n_probes':len(mapping),'n_genes':len(genes),'aggregation':'all corresponding probe means, then all repeated same-donor temporal arrays averaged; no diagnostic choice of array','fit_scope':'195 training subjects only; no test preprocessing fitting or hyperparameter selection'},indent=2))
print(pd.DataFrame(rows).query("model=='Logistic_L2' and selector=='all'").to_string(index=False))
