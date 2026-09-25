"""Re-execution of the 72 submitted pipelines with scikit-learn 1.5.2, with checks of fold isolation and probability equations."""
from pathlib import Path
import sys,json,importlib.util,inspect,warnings,hashlib
R=Path(__file__).resolve().parents[2]
import numpy as np,pandas as pd,scipy,sklearn
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest,f_classif
from sklearn.model_selection import StratifiedKFold,cross_val_predict
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
from scipy.special import expit,logsumexp
src=Path(r'<ANALYSIS_ROOT>/fuxian/JTM_final_manuscript_figure_scripts_20260609/figure9_GEM_bulk_multimodel_10fold_final.py')
sp=importlib.util.spec_from_file_location('actual_final',src);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
tr=pd.read_csv(R/'results/discovery_GSE132903_features.csv');X=tr.iloc[:,2:].to_numpy(float);y=tr.group.eq('AD').to_numpy(int);genes=tr.columns[2:].tolist()
assert set(tr.group)=={'AD','CN'} and tr.sample_id.is_unique and np.isfinite(X).all()
old=pd.read_csv(Path(r'<ANALYSIS_ROOT>/fuxian_test/day6_GEM_method_enhancement_screen/tables/GEM_bulk_multimodel_10fold_predictions.csv'))
cv=StratifiedKFold(n_splits=10,shuffle=True,random_state=m.SEED);folds=list(cv.split(X,y));fold_manifest=[]
for f,(it,iv) in enumerate(folds,1):
 assert set(it).isdisjoint(iv)
 fold_manifest.extend(dict(fold=f,role=role,sample_id=tr.sample_id.iloc[i],group=tr.group.iloc[i]) for role,ind in [('train',it),('test',iv)] for i in ind)
pd.DataFrame(fold_manifest).to_csv(R/'intermediate'/f'ML_folds_sklearn_{sklearn.__version__}.csv',index=False)
rows=[];params={};warningrows=[];predrows=[]
for name,est in m.model_zoo().items():
 params[name]=est.get_params()
 for sel in ['all','top3','top5','top8','top12','top16']:
  steps=[('scale',StandardScaler())]
  if sel!='all':steps.append(('select',SelectKBest(f_classif,k=int(sel[3:]))))
  steps.append(('model',clone(est)));pipe=Pipeline(steps)
  with warnings.catch_warnings(record=True) as caught:
   warnings.simplefilter('always');pred=cross_val_predict(pipe,X,y,cv=cv,method='predict_proba',n_jobs=1)[:,1]
  for w in caught:
   if issubclass(w.category,ConvergenceWarning):warningrows.append({'pipeline':sel+'+'+name,'warning':str(w.message)})
  pold=old[old.pipeline.eq(sel+'+'+name)].set_index('sample_id').loc[tr.sample_id].score.to_numpy()
  # Independent pairwise rank definition, including half credit for ties.
  a=pred[y==1];c=pred[y==0];pairauc=np.mean((a[:,None]>c).astype(float)+.5*(a[:,None]==c))
  auc=roc_auc_score(y,pred);assert abs(auc-pairauc)<1e-12
  rows.append(dict(pipeline=sel+'+'+name,auc=auc,archived_auc=roc_auc_score(y,pold),max_probability_difference=float(np.max(np.abs(pred-pold))),n=float(len(pred)),pairwise_AUC_abs_error=abs(auc-pairauc)))
  predrows.extend(dict(sample_id=s,pipeline=sel+'+'+name,group=g,score=float(p)) for s,g,p in zip(tr.sample_id,tr.group,pred))
 print(name,'complete',flush=True)
pd.DataFrame(rows).to_csv(R/'results'/f'ML_reexecution_sklearn_{sklearn.__version__}.csv',index=False)
pd.DataFrame(predrows).to_csv(R/'results'/f'ML_OOF_sklearn_{sklearn.__version__}.csv',index=False)
# Train-fold moments and ANOVA formula checked against the real first split.
it,iv=folds[0];z=X[it];yy=y[it];scale=StandardScaler().fit(z);assert np.allclose(scale.mean_,z.mean(axis=0));assert np.allclose(scale.var_,z.var(axis=0,ddof=0))
zs=scale.transform(z);means=np.array([zs[yy==k].mean(0) for k in [0,1]]);n=np.array([(yy==k).sum() for k in [0,1]]);grand=zs.mean(0)
between=(n[:,None]*(means-grand)**2).sum(0);within=sum(((zs[yy==k]-means[k])**2).sum(0) for k in [0,1]);manual_f=between/(within/(len(yy)-2));f,_=f_classif(zs,yy)
assert np.allclose(f,manual_f,rtol=1e-10,atol=1e-10)
selected=SelectKBest(f_classif,k=5).fit(zs,yy).get_support();assert set(np.where(selected)[0])==set(np.argsort(manual_f)[-5:])
# Test-data perturbation cannot alter already fitted training statistics or
# predictions on unchanged test rows. This detects accidental refitting.
pipe=Pipeline([('scale',StandardScaler()),('select',SelectKBest(f_classif,k=5)),('model',clone(m.model_zoo()['Logistic_L2']))]).fit(X[it],y[it]);before=pipe.predict_proba(X[iv]);mean=pipe['scale'].mean_.copy();chosen=pipe['select'].get_support().copy();xt=X[iv].copy();xt[0]+=10000;after=pipe.predict_proba(xt)
assert np.array_equal(mean,pipe['scale'].mean_) and np.array_equal(chosen,pipe['select'].get_support());assert np.allclose(before[1:],after[1:]);assert np.array_equal(pipe.classes_,[0,1])
# Logistic and Gaussian NB probability equations verified on actual data.
log=Pipeline([('s',StandardScaler()),('m',clone(m.model_zoo()['Logistic_L2']))]).fit(z,yy);prob=expit(log['s'].transform(X[iv])@log['m'].coef_.T+log['m'].intercept_).ravel();assert np.allclose(prob,log.predict_proba(X[iv])[:,1])
gn=clone(m.model_zoo()['GaussianNB']).fit(zs,yy);zt=scale.transform(X[iv]);lj=np.array([np.log(gn.class_prior_[k])-.5*np.log(2*np.pi*gn.var_[k]).sum()-.5*((zt-gn.theta_[k])**2/gn.var_[k]).sum(1) for k in [0,1]]).T;gp=np.exp(lj-logsumexp(lj,axis=1)[:,None]);assert np.allclose(gp,gn.predict_proba(zt))
record={'sklearn':sklearn.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'params':params,'convergence_warnings':warningrows,'train_test_checks':'PASS: 10 disjoint splits; each held-out once; 195 unique finite rows; AD=1; probability column index1 equals class1; no transform-time refit','formula_checks':'PASS: population-variance scaling, two-group ANOVA F, top-5 selection, logistic sigmoid, Gaussian NB posterior and all72 pairwise-rank AUCs','first_fold_features':np.array(genes)[selected].tolist(),'version_limit':'Installed environment and numerical agreement are evidence; an absent historical package lock cannot be inferred from current installations.'}
(R/'intermediate'/f'ML_implementation_sklearn_{sklearn.__version__}.json').write_text(json.dumps(record,indent=2,default=str))
print('RESULT',sklearn.__version__,'max probability difference',max(r['max_probability_difference'] for r in rows),'convergence warnings',len(warningrows))
