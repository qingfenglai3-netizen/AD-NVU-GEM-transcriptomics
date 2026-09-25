from pathlib import Path
import json,sys,importlib.util
import sklearn
assert sklearn.__version__=='1.5.2', 'Use the pinned reproducible sklearn 1.5.2 environment.'
import numpy as np,pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.metrics import roc_auc_score
from sklearn.base import clone
sys.stdout.reconfigure(encoding='utf-8')
R=Path(__file__).resolve().parents[1];OUT=R/'results';OUT.mkdir(exist_ok=True)
src=R/'analysis/11_translational/11a_bulk_multimodel_cross_validation.py'
spec=importlib.util.spec_from_file_location('original_models',src);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
tr=pd.read_csv(R/'data/discovery_GSE132903_labelblind_mean_features.csv');te=pd.read_csv(R/'data/external_GSE36980_features.csv')
genes=tr.columns[2:].tolist();X=tr[genes].to_numpy();y=tr.group.eq('AD').astype(int).to_numpy()
assert set(te.sample_id).isdisjoint(tr.sample_id)
assert len(tr)==195 and y.sum()==97
assert set(tr.group)=={'AD','CN'} and set(te.group)=={'AD','CN'}
assert tr.sample_id.is_unique and te.sample_id.is_unique
assert np.isfinite(X).all() and np.isfinite(te[genes].to_numpy()).all()
rows=[];predrows=[];rng=np.random.default_rng(20260910)
for name,est in mod.model_zoo().items():
 for sel in ['all','top3','top5','top8','top12','top16']:
  steps=[('scale',StandardScaler())]
  if sel!='all':steps.append(('select',SelectKBest(f_classif,k=int(sel[3:]))))
  steps.append(('model',clone(est)));pipe=Pipeline(steps).fit(X,y)
  for reg in ['TC']:
   d=te[te.region.eq(reg)].copy();yy=d.group.eq('AD').astype(int).to_numpy();pred=pipe.predict_proba(d[genes].to_numpy())[:,1]
   auc=roc_auc_score(yy,pred);boot=[]
   if name=='Logistic_L2' and sel=='all':
    for i in range(2000):
     idx=np.r_[rng.choice(np.where(yy==0)[0],sum(yy==0)),rng.choice(np.where(yy==1)[0],sum(yy==1))];boot.append(roc_auc_score(yy[idx],pred[idx]))
   ci=np.quantile(boot,[.025,.975]) if boot else [np.nan,np.nan]
   rows.append(dict(model=name,selector=sel,region=reg,n_AD=int(sum(yy)),n_CN=int(sum(yy==0)),auc=auc,ci_low=ci[0],ci_high=ci[1]))
   predrows.extend(dict(model=name,selector=sel,region=reg,sample_id=s,group=g,prediction=float(p)) for s,g,p in zip(d.sample_id,d.group,pred))
  print(name,sel,'done',flush=True)
pd.DataFrame(rows).to_csv(OUT/'external_model_performance.csv',index=False)
pd.DataFrame(predrows).to_csv(OUT/'external_model_predictions.csv',index=False)
print(pd.DataFrame(rows).query("model=='Logistic_L2' and selector=='all'").to_string(index=False))
