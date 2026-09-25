from pathlib import Path
import sys,json
import pandas as pd,numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
R=Path(r'<WORKDIR>');P=Path(r'<ANALYSIS_ROOT>\fuxian_test')
sys.stdout.reconfigure(encoding='utf-8')
d=pd.read_csv(P/'day3_2_sctour/final_submission_package/source_data/Figure_9D_source_data.csv')
dm=pd.read_csv(R/'results/donor_APOE_manifest.csv');d=d.merge(dm,left_on='sample_id',right_on='SampleID',suffixes=('','_manifest'),validate='m:1')
rows=[]
for k,x in d.groupby('trajectory'):
 a=x.loc[x.Group.eq('AD'),'median_pseudotime'].to_numpy();c=x.loc[x.Group.eq('CN'),'median_pseudotime'].to_numpy()
 p=mannwhitneyu(a,c,alternative='two-sided').pvalue
 rng=np.random.default_rng(20260910);delta=np.array([np.median(rng.choice(a,len(a)))-np.median(rng.choice(c,len(c))) for _ in range(5000)])
 rows.append(dict(trajectory=k,n_AD=len(a),n_CN=len(c),median_AD=np.median(a),median_CN=np.median(c),delta=np.median(a)-np.median(c),p_value=p,ci_low=np.quantile(delta,.025),ci_high=np.quantile(delta,.975),n_cells_AD=int(x.loc[x.Group.eq('AD'),'n_cells'].sum()),n_cells_CN=int(x.loc[x.Group.eq('CN'),'n_cells'].sum())))
stats=pd.DataFrame(rows);stats['q_BH']=multipletests(stats.p_value,method='fdr_bh')[1];stats.to_csv(R/'results/trajectory_statistics.csv',index=False);d.to_csv(R/'results/trajectory_donor_values.csv',index=False)
print(stats.to_string(index=False))
sr=[]
for k,x in d.groupby('trajectory'):
 for model in ['cohort_adjusted','APOE_sex_adjusted']:
  z=x if model=='cohort_adjusted' else x[x.Dataset.eq('D1')]
  formula='median_pseudotime ~ C(Group) + C(Dataset)' if model=='cohort_adjusted' else 'median_pseudotime ~ C(Group) + APOE2_dosage + APOE4_dosage + C(sex)'
  f=smf.ols(formula,z).fit(cov_type='HC3');term='C(Group)[T.CN]';sr.append(dict(trajectory=k,model=model,n=len(z),AD_minus_CN=-f.params[term],p_value=f.pvalues[term]))
pd.DataFrame(sr).to_csv(R/'results/trajectory_covariate_sensitivity.csv',index=False)
sea=pd.read_csv(P/'SEAAD_ADNC_vascular_GEM_validation/SEAAD_ADNC_vascular_sample_expression.csv');sea=sea[sea.gene.eq('GEM')];rs=[]
for label,x in sea.groupby('dataset'):
 for comp,cases,controls in [('dementia_vs_normal',['dementia'],['normal']),('ADNC_High_vs_NotADRef',['High'],['Not AD','Reference'])]:
  col='disease' if comp=='dementia_vs_normal' else 'ADNC';a=x[x[col].isin(cases)].expr.to_numpy();c=x[x[col].isin(controls)].expr.to_numpy()
  rs.append(dict(dataset=label,comparison=comp,n_case=len(a),n_control=len(c),mean_case=np.mean(a),mean_control=np.mean(c),delta=np.mean(a)-np.mean(c),p_value=mannwhitneyu(a,c,alternative='two-sided').pvalue))
rs=pd.DataFrame(rs);rs['q_BH']=rs.groupby('comparison').p_value.transform(lambda p:multipletests(p,method='fdr_bh')[1]);rs.to_csv(R/'results/SEAAD_GEM_statistics.csv',index=False);sea.to_csv(R/'results/SEAAD_GEM_donor_values.csv',index=False);print(rs.to_string(index=False))
