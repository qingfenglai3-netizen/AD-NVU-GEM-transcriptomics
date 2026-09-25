"""Spatial co-expression scores for the prioritized ligand-receptor axes.
Receptor complexes require every subunit and use the minimum subunit expression.
"""
from pathlib import Path
import sys,json,re,itertools
import h5py,numpy as np,pandas as pd
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
sys.stdout.reconfigure(encoding='utf-8')
R=Path(r'<WORKDIR>');P=Path(r'<ANALYSIS_ROOT>\fuxian_test');src=Path(r'<WORKDIR>\AD-NVU-GEM-transcriptomics\derived_data')
c=pd.read_csv(src/'Figure_6C_v3_NVU_BBB_focused_AD_enriched_LR.csv');d=pd.read_csv(src/'Figure_6D_v3_activated_capillary_LR.csv');d['source']='CapEC_ActHigh'
c['panel']='6C';d['panel']='6D';axes=pd.concat([c,d],ignore_index=True)[['source','target','ligand','receptor','panel']].drop_duplicates(['source','target','ligand','receptor'])
axes['LR_Pair']=axes.ligand+':'+axes.receptor.str.replace('_','/',regex=False)
axes.to_csv(R/'results/spatial_prespecified_axes.csv',index=False)
propall=pd.read_csv(P/'day4_1_spatial/tables_rctd/TableS4_spot_level_RCTD_proportions.csv')
def decode(x):return np.array([v.decode() if isinstance(v,bytes) else str(v) for v in x])
def scale(v):
 lo,hi=np.percentile(v,[5,95]);return np.clip((v-lo)/(hi-lo),0,1) if hi>lo else v
def ctcols(label):
 if label.startswith('CapEC_'):return ['Endo_Capillary']
 if label.startswith('Cerebrovascular_'):return ['Endo_'+label.split('_',1)[1]]
 if label in ['Astrocytes','Astrocyte']:return ['Astro_Homeostatic','Astro_Intermediate','Astro_Reactive']
 if label in ['Microglia']:return ['Micro_DAM','Micro_Homeostatic']
 return [label]
records=[]
for sample,group in {'1-1':'CN','18-64':'CN','2-5':'CN','2-3':'AD','2-8':'AD','T4857':'AD'}.items():
 root=Path(r'<ANALYSIS_ROOT>\spatial_data\GSE220442\counts_and_images')/sample
 with h5py.File(root/'filtered_feature_bc_matrix.h5') as f:
  m=f['matrix'];X=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=m['shape'][:]);genes=decode(m['features/name'][:]);bar=decode(m['barcodes'][:])
 keep=(X>0).sum(axis=0).A1>=200;X=X[:,keep];bar=bar[keep];kg=(X>0).sum(axis=1).A1>=3;X=X[kg];genes=genes[kg];libs=X.sum(axis=0).A1
 pos=pd.read_csv(root/'spatial/tissue_positions_list.csv',header=None,names=['barcode','in_tissue','array_row','array_col','pixel_row','pixel_col']).set_index('barcode').loc[bar]
 idx=cKDTree(pos[['pixel_col','pixel_row']].to_numpy()).query(pos[['pixel_col','pixel_row']].to_numpy(),k=9)[1][:,1:]
 prop=propall[propall.sample_id.eq(sample)].set_index('spot_id').reindex(bar).fillna(0)
 cache={}
 def vector(g):
  if g not in cache:
   ii=np.where(genes==g)[0];cache[g]=np.log1p(X[ii[0]].toarray().ravel()/libs*1e4) if len(ii) else np.zeros(len(bar))
  return cache[g]
 for a in axes.itertuples():
  ls=re.split('[_/]',a.ligand);rs=re.split('[_/]',a.receptor);lg=np.minimum.reduce([vector(g) for g in ls]);rg=np.minimum.reduce([vector(g) for g in rs])
  sc=ctcols(a.source);tc=ctcols(a.target);assert all(k in prop for k in sc+tc),(a.source,a.target,sc,tc)
  local=scale(lg)*scale(rg)*(.5*prop[sc].sum(axis=1).to_numpy()+.5*prop[tc].sum(axis=1).to_numpy())
  local=scale(scale(.6*local+.4*local[idx].mean(axis=1)))
  records.append(dict(sample=sample,group=group,source=a.source,target=a.target,LR_Pair=a.LR_Pair,panel=a.panel,n_spots=len(bar),ligand_detected_fraction=float(np.mean(lg>0)),receptor_detected_fraction=float(np.mean(rg>0)),all_subunits_available=all(g in genes for g in ls+rs),missing_subunits=';'.join(g for g in ls+rs if g not in genes),mean_spatial_score=float(np.mean(local))))
 print(sample,'axes',len(axes),'spots',len(bar),flush=True)
out=pd.DataFrame(records);out.to_csv(R/'results/spatial_axis_section_scores.csv',index=False)
ss=[]
for keys,z in out.groupby(['source','target','LR_Pair','panel']):
 delta=z.loc[z.group.eq('AD'),'mean_spatial_score'].mean()-z.loc[z.group.eq('CN'),'mean_spatial_score'].mean();v=z.mean_spatial_score.to_numpy()
 perms=[v[list(i)].mean()-np.delete(v,list(i)).mean() for i in itertools.combinations(range(6),3)]
 p=np.mean(np.abs(perms)>=abs(delta)-1e-12)
 ss.append(dict(zip(['source','target','LR_Pair','panel'],keys),delta_AD_CN=delta,permutation_p=p,n_sections=6,min_receptor_detection=z.receptor_detected_fraction.min(),all_subunits_available=bool(z.all_subunits_available.all())))
ss=pd.DataFrame(ss);ss['q_BH']=multipletests(ss.permutation_p,method='fdr_bh')[1];ss.to_csv(R/'results/spatial_axis_summary.csv',index=False)
print(ss[ss.LR_Pair.str.startswith(('ADM:','TGFB2:','TGFB1:'))].to_string(index=False))
