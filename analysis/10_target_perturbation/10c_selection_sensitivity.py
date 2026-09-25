from pathlib import Path
import sys,json,re,hashlib
import pandas as pd,numpy as np
from sklearn.metrics import roc_auc_score
R=Path(r'<WORKDIR>');P=Path(r'<ANALYSIS_ROOT>\fuxian_test');S=Path(r'<SUBMISSION_ARCHIVE>\14_PLOS_ONE_AD_GEM_upload_20260712');repo=S/'06_Code_Data_Repository/AD-NVU-GEM-transcriptomics'
sys.stdout.reconfigure(encoding='utf-8')
dd=repo/'derived_data';pr=pd.read_csv(dd/'GEM_Figure9_Table_bulk_multimodel_10fold_predictions.csv');print('prediction columns',pr.columns.tolist())
print(pr.head(2).to_string(index=False))
for stem in ['glia_to_vascular_ligand_activities','glia_to_vascular_target_priority']:
 a=pd.read_csv(dd/f'{stem}.csv');b=pd.read_csv(R/'results/nichenet'/f'{stem}.csv');pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=1e-8);print(stem,'REPRODUCED')
la=pd.read_csv(R/'results/nichenet/glia_to_vascular_ligand_activities.csv');la['AUPR_rank']=np.arange(1,len(la)+1)
prior=set(pd.read_csv(R/'results/nichenet/prior_ligand_universe.csv').ligand)
full=pd.read_csv(dd/'Figure_6_v3_full_filtered_LIANA_LR_source_data.csv')
rows=[]
for g in ['TGFB1','TGFB2','COL1A1','HBEGF','CXCL10','SLIT2','ADM','VEGFC','ANGPT2','EGF','FGF1','VEGFA','FGF2','ANGPT1','PDGFB','SELPLG','IL1B','TNFSF10','IL18']:
 z=full[full.ligand.eq(g)];f=z[(z.delta_AD_vs_CN>0)&z.source_family.isin(['Astrocyte','Microglia'])&z.target_family.eq('Vascular')];q=la[la.test_ligand.eq(g)]
 rows.append(dict(ligand=g,n_all_edges=len(z),n_glial_to_vascular_AD_increased_edges=len(f),in_prior=g in prior,eligible=bool(len(f) and g in prior),AUPR_rank=int(q.AUPR_rank.iloc[0]) if len(q) else np.nan,aupr=float(q.aupr.iloc[0]) if len(q) else np.nan,source_families=q.source_families.iloc[0] if len(q) else '',in_top20=g in set(la.head(20).test_ligand)))
au=pd.DataFrame(rows);au.to_csv(R/'results/ligand_selection.csv',index=False);print(au.to_string(index=False))
summary=[]
for c in ['Cerebrovascular','Astrocytes','Microglia']:
 for model in ['diagnosis_only','cohort_adjusted','D1_APOE_sex_adjusted']:
  x=pd.read_csv(R/'results'/f'{c}_MiloR_{model}.csv')
  for sub,z in x.groupby('Subtype'):
   summary.append(dict(compartment=c,Subtype=sub,model=model,n_neighborhoods=len(z),n_AD_FDR10=int(((z.FDR<.1)&(z.logFC>0)).sum()),n_CN_FDR10=int(((z.FDR<.1)&(z.logFC<0)).sum()),min_FDR=z.FDR.min()))
pd.DataFrame(summary).to_csv(R/'results/MiloR_subtype_sensitivity.csv',index=False)
pd.read_csv(P/'JTM_manuscript/sandbox/figure8_target_reselection/NicheNet_target_side_topn_sensitivity_interest.csv').to_csv(R/'results/NicheNet_topN_sensitivity.csv',index=False)
pd.read_csv(P/'JTM_manuscript/sandbox/figure8_target_reselection/NicheNet_target_side_metric_sensitivity_interest.csv').to_csv(R/'results/NicheNet_metric_sensitivity.csv',index=False)
