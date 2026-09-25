"""Out-of-fold performance of the 72 pipelines with probes averaged per gene without labels."""
from pathlib import Path
import importlib.util,warnings,json
import numpy as np,pandas as pd,sklearn
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
assert sklearn.__version__=='1.5.2'
R=Path(__file__).resolve().parents[1];OUT=R/'results';OUT.mkdir(exist_ok=True)
src=R/'analysis/11_translational/11a_bulk_multimodel_cross_validation.py'
sp=importlib.util.spec_from_file_location('original',src);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
d=pd.read_csv(R/'data/discovery_GSE132903_labelblind_mean_features.csv');X=d.iloc[:,2:].to_numpy(float);y=d.group.eq('AD').to_numpy(int)
assert d.sample_id.is_unique and set(d.group)=={'AD','CN'} and np.isfinite(X).all()
cv=StratifiedKFold(10,shuffle=True,random_state=m.SEED);rows=[];prs=[];ws=[]
for name,est in m.model_zoo().items():
 for sel in ['all','top3','top5','top8','top12','top16']:
  steps=[('scale',StandardScaler())]
  if sel!='all':steps.append(('select',SelectKBest(f_classif,k=int(sel[3:]))))
  steps.append(('model',clone(est)))
  with warnings.catch_warnings(record=True) as ww:
   warnings.simplefilter('always');pred=cross_val_predict(Pipeline(steps),X,y,cv=cv,method='predict_proba',n_jobs=1)[:,1]
  ws.extend(dict(pipeline=sel+'+'+name,message=str(w.message)) for w in ww if issubclass(w.category,ConvergenceWarning))
  rows.append(dict(model=name,selector=sel,pipeline=sel+'+'+name,auc=roc_auc_score(y,pred),cv='10fold_labelblind_probe_mean'))
  prs.extend(dict(sample_id=s,group=g,pipeline=sel+'+'+name,score=p) for s,g,p in zip(d.sample_id,d.group,pred))
 print(name,flush=True)
pd.DataFrame(rows).to_csv(OUT/'ML_labelblind_OOF_summary.csv',index=False)
pd.DataFrame(prs).to_csv(OUT/'ML_labelblind_OOF_predictions.csv',index=False)
print(pd.DataFrame(rows).auc.agg(['min','max']));print('convergence warnings',len(ws))
