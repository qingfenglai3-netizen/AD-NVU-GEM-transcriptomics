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
R=Path(__file__).resolve().parents[1];OUT=R/'results';OUT.mkdir(exist_ok=True)
tr=pd.read_csv(R/'data/discovery_GSE132903_labelblind_mean_features.csv');genes=tr.columns[2:].tolist()
d=pd.read_csv(R/'data/GSE118553_TC_donor_features.csv')
src=R/'analysis/11_translational/11a_bulk_multimodel_cross_validation.py';sp=importlib.util.spec_from_file_location('original',src);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
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
pd.DataFrame(rows).to_csv(OUT/'GSE118553_TC_model_performance.csv',index=False);pd.DataFrame(prs).to_csv(OUT/'GSE118553_TC_predictions.csv',index=False)
