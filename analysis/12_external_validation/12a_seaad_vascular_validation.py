import h5py, os, pandas as pd, numpy as np, scipy.sparse as sp
from scipy.stats import mannwhitneyu, spearmanr
files={
 'MTG_Endothelial':r'/path/to/public_data/SEAAD_cellxgene_small_vascular/SEAAD_MTG_Endothelial.h5ad',
 'DLPFC_Endothelial':r'/path/to/public_data/SEAAD_cellxgene_small_vascular/SEAAD_DLPFC_Endothelial.h5ad',
 'MTG_VLMC':r'/path/to/public_data/SEAAD_cellxgene_small_vascular/SEAAD_MTG_VLMC.h5ad',
 'DLPFC_VLMC':r'/path/to/public_data/SEAAD_cellxgene_small_vascular/SEAAD_DLPFC_VLMC.h5ad'}
outdir=r'/path/to/project/results/12_seaad_validation'
os.makedirs(outdir,exist_ok=True)
def read_enc(group,key):
    g=group[key]
    if isinstance(g,h5py.Group) and 'categories' in g and 'codes' in g:
        cats=np.array([x.decode() if isinstance(x,bytes) else str(x) for x in g['categories'][:]],dtype=object); codes=g['codes'][:]
        return np.array([cats[i] if i>=0 else None for i in codes],dtype=object)
    arr=g[:]; return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in arr],dtype=object)
def symbols(f):
    for c in ['feature_name','gene_name','_index']:
        if c in f['var']:
            arr=read_enc(f['var'],c)
            if any(g in set(arr) for g in ['GEM','PDGFRB','CLDN5']): return arr
    return read_enc(f['var'],'_index')
def expr_cols(f,idx):
    X=f['X']
    if isinstance(X,h5py.Group):
        mat=sp.csr_matrix((X['data'][:],X['indices'][:],X['indptr'][:]), shape=tuple(X.attrs['shape']))
        return mat[:,idx].toarray()
    return X[:,idx]
allrows=[]; allpb=[]
for label,p in files.items():
    with h5py.File(p,'r') as f:
        obs_cols=['ADNC','Braak stage','CERAD score','Thal phase','Cognitive status','disease','donor_id','sex','Age at death','PMI','Subclass','Supertype','cell_type','tissue']
        obs={c:read_enc(f['obs'],c) for c in obs_cols if c in f['obs']}
        obsdf=pd.DataFrame(obs)
        sym=symbols(f); genes=['GEM','NOTCH3','PDGFRB','PTH1R','CLDN5','FLT1','VWF','PECAM1','ACTA2','RGS5','KCNJ8','ABCC9','PDGFB','TEK','TIE1','EDN1','ADM']
        gidx=[np.where(sym==g)[0][0] for g in genes if np.any(sym==g)]; gp=[sym[i] for i in gidx]
        E=expr_cols(f,gidx)
    for gi,g in enumerate(gp):
        tmp=obsdf[['donor_id','ADNC','Braak stage','CERAD score','Thal phase','disease','sex']].copy(); tmp['expr']=E[:,gi]
        pb=tmp.groupby(['donor_id','ADNC','Braak stage','CERAD score','Thal phase','disease','sex']).agg(expr=('expr','mean'),n_cells=('expr','size')).reset_index()
        pb=pb[pb.n_cells>=3].copy(); pb['dataset']=label; pb['gene']=g
        allpb.append(pb)
        comparisons={
            'ADNC_High_vs_NotADRef': (pb['ADNC'].eq('High'), pb['ADNC'].isin(['Not AD','Reference'])),
            'ADNC_HighInt_vs_NotADRef': (pb['ADNC'].isin(['High','Intermediate']), pb['ADNC'].isin(['Not AD','Reference'])),
            'dementia_vs_normal': (pb['disease'].eq('dementia'), pb['disease'].eq('normal')),
        }
        for cname,(case,ctrl) in comparisons.items():
            ncase=case.sum(); nctrl=ctrl.sum()
            if ncase>=3 and nctrl>=3:
                y=pb.loc[case,'expr']; x=pb.loc[ctrl,'expr']; pv=mannwhitneyu(y,x,alternative='two-sided').pvalue
                allrows.append(dict(dataset=label,gene=g,comparison=cname,n_case=int(ncase),n_control=int(nctrl),mean_case=float(y.mean()),mean_control=float(x.mean()),delta=float(y.mean()-x.mean()),p_wilcox=float(pv)))
        # ADNC ordinal Spearman
        scoremap={'Reference':0,'Not AD':0,'Low':1,'Intermediate':2,'High':3}
        score=pb['ADNC'].map(scoremap)
        ok=score.notna()
        if ok.sum()>=6:
            rho,pv=spearmanr(score[ok], pb.loc[ok,'expr'])
            allrows.append(dict(dataset=label,gene=g,comparison='ADNC_ordinal_spearman',n_case=int(ok.sum()),n_control=0,mean_case=np.nan,mean_control=np.nan,delta=float(rho),p_wilcox=float(pv)))
stats=pd.DataFrame(allrows)
if len(stats):
    stats['q_by_gene_comparison']=np.nan
    for key,grp in stats.groupby(['gene','comparison']):
        p=grp.p_wilcox.fillna(1).values; order=np.argsort(p); m=len(p); spv=p[order]
        bh=np.minimum.accumulate((spv*m/np.arange(1,m+1))[::-1])[::-1]; arr=np.empty(m); arr[order]=np.minimum(bh,1)
        for idx,val in zip(grp.index,arr): stats.loc[idx,'q_by_gene_comparison']=val
stats.to_csv(os.path.join(outdir,'SEAAD_ADNC_vascular_GEM_module_stats.csv'),index=False)
pbdf=pd.concat(allpb,ignore_index=True); pbdf.to_csv(os.path.join(outdir,'SEAAD_ADNC_vascular_sample_expression.csv'),index=False)
print(stats[stats.gene.eq('GEM')].sort_values('p_wilcox').to_string(index=False))
