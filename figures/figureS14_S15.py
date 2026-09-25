"""S14 Fig (external classification) and S15 Fig (APOE sensitivity and spatial scores) from the tables in data/."""
from pathlib import Path
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve
R=Path(__file__).resolve().parents[1];OUT=R/'results'/'figures';OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 for ext in ['png','tif','pdf','svg']:
  kw={'dpi':600,'bbox_inches':'tight','facecolor':'white'}
  if ext=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
  fig.savefig(OUT/f'{name}.{ext}',**kw)
 plt.close(fig)
perf=pd.read_csv(R/'data/external_model_performance.csv');pred=pd.read_csv(R/'data/external_model_predictions.csv')
perf['cohort']='GSE36980';pred['cohort']='GSE36980'
perf=pd.concat([perf,pd.read_csv(R/'data/GSE118553_TC_model_performance.csv')],ignore_index=True)
pred=pd.concat([pred,pd.read_csv(R/'data/GSE118553_TC_predictions.csv')],ignore_index=True)
models=['Logistic_L2','Logistic_L1','ElasticNet_logistic','Linear_SVM','RBF_SVM','RandomForest','ExtraTrees','GradientBoosting','DecisionTree','KNN_3','KNN_5','GaussianNB']
selectors=['all','top3','top5','top8','top12','top16']
fig=plt.figure(figsize=(6.7,6.0));gs=fig.add_gridspec(2,2,height_ratios=[1,.8],hspace=.52,wspace=.55)
for j,cohort in enumerate(['GSE118553','GSE36980']):
 ax=fig.add_subplot(gs[0,j]);mat=perf[perf.cohort.eq(cohort)].pivot(index='model',columns='selector',values='auc').reindex(index=models,columns=selectors)
 sns.heatmap(mat,ax=ax,cmap='RdYlBu_r',vmin=0,vmax=1,annot=True,fmt='.2f',annot_kws={'fontsize':5.8},cbar=(j==1),cbar_kws={'label':'External AUC','shrink':.6},linewidths=.2)
 ax.set_xlabel('');ax.set_ylabel('');ax.set_title(('A  ' if j==0 else 'B  ')+cohort+' temporal cortex',loc='left',fontweight='bold');ax.tick_params(axis='x',rotation=35);ax.set_yticklabels([x.replace('_',' ') for x in models],rotation=0)
 ax=fig.add_subplot(gs[1,j]);z=pred[pred.cohort.eq(cohort)&pred.model.eq('Logistic_L2')&pred.selector.eq('all')];f,t,_=roc_curve(z.group.eq('AD'),z.prediction);st=perf[perf.cohort.eq(cohort)&perf.model.eq('Logistic_L2')&perf.selector.eq('all')].iloc[0]
 ax.plot(f,t,color='#E64B35' if j==0 else '#4DBBD5',lw=1.5);ax.plot([0,1],[0,1],ls='--',color='.65',lw=.7)
 ax.set(xlim=(0,1),ylim=(0,1.02),xlabel='False positive rate',ylabel='True positive rate');ax.set_title(('C' if j==0 else 'D')+'  Fixed all-gene L2 model',loc='left',fontweight='bold');ax.text(.98,.05,f'AUC {st.auc:.3f}\n95% CI {st.ci_low:.3f}–{st.ci_high:.3f}\n{st.n_AD} AD, {st.n_CN} controls',transform=ax.transAxes,ha='right',fontsize=7)
save(fig,'S14_Fig')
# APOE estimates: unshrunk coefficients with approximate QL Wald intervals.
x=pd.read_csv(R/'data/APOE_nested_stratified_GEM.csv')
fig=plt.figure(figsize=(6.7,7.8));gs=fig.add_gridspec(2,1,height_ratios=[1,2.0],hspace=.40,left=.41,right=.96)
ax=fig.add_subplot(gs[0]);labels=['Diagnosis','Sex','Sex + APOE dosage','Sex + ε4 carrier','AD × ε4 interaction','ε4 noncarriers','ε4 carriers']
y=np.arange(len(x));ax.errorbar(x.unshrunk_logFC,y,xerr=np.vstack([x.unshrunk_logFC-x.approx_QL_CI_low,x.approx_QL_CI_high-x.unshrunk_logFC]),fmt='o',ms=3,color='#4DBBD5',ecolor='#3C5488',elinewidth=.9,capsize=2)
ax.axvline(0,color='.5',ls='--',lw=.6);ax.set_yticks(y);ax.set_yticklabels(labels);ax.invert_yaxis();ax.set_xlabel('Unshrunk log2 coefficient (approximate 95% CI)');ax.set_title('A  Complete-case GEM sensitivity',loc='left',fontweight='bold')
z=pd.read_csv(R/'data/spatial_axis_section_scores.csv')
abbr={'Astro_Homeostatic':'Ast.H','Astro_Intermediate':'Ast.I','Astro_Reactive':'Ast.R','CapEC_ActHigh':'Cap.H','CapEC_ActMid':'Cap.M','CapEC_ActLow':'Cap.L','Cerebrovascular_Pericyte':'Per','Excitatory':'Exc','Inhibitory':'Inh','OPCs':'OPC','Oligodendrocytes':'Olig','Micro_DAM':'Mic.DAM'}
key=['source','target','LR_Pair'];z['axis']=z.source.replace(abbr)+'→'+z.target.replace(abbr)+' | '+z.LR_Pair
scorecol='mean_spatial_score'
samplecol=next(c for c in ['sample','section','sample_id'] if c in z.columns)
mat=z.pivot(index='axis',columns=samplecol,values=scorecol)
mat=mat.reindex(columns=[c for c in ['1-1','18-64','2-5','2-3','2-8','T4857'] if c in mat.columns])
scaled=mat.sub(mat.min(axis=1),axis=0).div((mat.max(axis=1)-mat.min(axis=1)).replace(0,1),axis=0)
ax=fig.add_subplot(gs[1]);sns.heatmap(scaled,ax=ax,cmap='YlGnBu',vmin=0,vmax=1,cbar_kws={'label':'Within-axis scale','shrink':.4},linewidths=.3);ax.set_ylabel('');ax.set_xlabel('CN sections              AD sections');ax.tick_params(axis='y',labelsize=6);ax.tick_params(axis='x',rotation=30);ax.set_title('B  Spatial scores of Fig. 6 contexts',loc='left',fontweight='bold')
save(fig,'S15_Fig')
