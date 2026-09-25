from pathlib import Path
import pandas as pd,numpy as np,json,re,sys
sys.stdout.reconfigure(encoding='utf-8')
R=Path(__file__).resolve().parents[2];P=Path('<ANALYSIS_ROOT>/fuxian_test');J=P/'JTM_manuscript';src=J/'source_data';day=P/'day3_1_liana_go_reactome';t=day/'final_submission_package/tables';out=R/'results';A=R/'intermediate'
keys=['source','target','ligand','receptor'];steps=[]
def count(name,d):steps.append({'step':name,'n_rows':len(d),'n_ligands':d.ligand.nunique(),'n_pairs':len(d[['ligand','receptor']].drop_duplicates())})
ad=pd.read_csv(t/'Day3_1_State_LIANA_raw_AD.csv');cn=pd.read_csv(t/'Day3_1_State_LIANA_raw_CN.csv');count('LIANA AD pooled output',ad);count('LIANA CN pooled output',cn)
for d in [ad,cn]:assert d.aggregate_rank.notna().all() and np.allclose(d.communication_strength,-np.log10(np.maximum(d.aggregate_rank,1e-12)))
z=ad[keys+['communication_strength']].merge(cn[keys+['communication_strength']],on=keys,how='outer',suffixes=('_AD','_CN')).fillna(0);count('Outer join AD and CN',z)
z=z[(z.source!='Cerebrovascular_Arterial')&(z.target!='Cerebrovascular_Arterial')].copy();count('Exclude arterial contexts',z)
for group in ['AD','CN']:
 co=pd.read_csv(t/f'Day3_1_State_LIANA_cellcounts_after_downsampling_{group}.csv').set_index('label').n_cells
 for side in ['source','target']:z[side+'_n_'+group]=z[side].map(co).fillna(0)
z=z[z[[side+'_n_'+g for side in ['source','target'] for g in ['AD','CN']]].min(axis=1)>=100].copy();count('At least100 source and receiver nuclei in both groups',z)
f=pd.read_csv(src/'Figure_6_v3_full_filtered_LIANA_LR_source_data.csv');m=f.merge(z,on=keys,suffixes=('_archive','_recomputed'));assert len(m)==len(f)==len(z)
assert np.allclose(m.communication_strength_AD_archive,m.communication_strength_AD_recomputed) and np.allclose(m.communication_strength_CN_archive,m.communication_strength_CN_recomputed)
g=f[(f.delta_AD_vs_CN>0)&f.source_family.isin(['Astrocyte','Microglia'])&(f.target_family=='Vascular')];count('Positive AD minus CN glia to vascular edges',g)
prior=set(pd.read_csv(out/'nichenet/prior_ligand_universe.csv').ligand);g=g[g.ligand.isin(prior)];count('Ligands covered by NicheNet prior',g)
ca=pd.read_csv(out/'nichenet/glia_to_vascular_candidate_ligands.csv');la=pd.read_csv(out/'nichenet/glia_to_vascular_ligand_activities.csv');ta=pd.read_csv(out/'nichenet/glia_to_vascular_target_priority.csv');assert set(g.ligand)==set(ca.ligand_clean);la['AUPR_rank']=range(1,len(la)+1)
disp=la.assign(source_class=la.source_families.map({'Astrocyte':'Astrocyte','Microglia':'Microglia'})).dropna(subset=['source_class']).groupby('source_class',sort=False).head(5)
ca=ca.merge(la[['ligand_clean','aupr','AUPR_rank']],on='ligand_clean');ca['input_top20']=ca.AUPR_rank<=20;ca['displayed_Fig8A']=ca.ligand_clean.isin(disp.ligand_clean)
ca.to_csv(out/'NicheNet_all186_ligand_provenance.csv',index=False);g.to_csv(out/'NicheNet_eligible_LIANA_edges.csv',index=False)
for panel in ['B','C','D']:
 name={'B':'Figure_6B_v3_global_cross_compartment_top_LR.csv','C':'Figure_6C_v3_NVU_BBB_focused_AD_enriched_LR.csv','D':'Figure_6D_v3_activated_capillary_LR.csv'}[panel]
 d=pd.read_csv(src/name)
 if panel=='D':d['source']='CapEC_ActHigh'
 d['panel']='6'+panel
 d['ligand_eligible_NicheNet']=d.ligand.isin(ca.ligand_clean);d['AUPR_rank']=d.ligand.map(la.set_index('ligand_clean').AUPR_rank);d['ligand_displayed_Fig8A']=d.ligand.isin(disp.ligand_clean)
 d['exact_edge_eligible']=pd.MultiIndex.from_frame(d[keys]).isin(pd.MultiIndex.from_frame(g[keys]))
 if panel=='B':cross=d
 else:cross=pd.concat([cross,d],ignore_index=True)
cross.to_csv(out/'Fig6_display_to_NicheNet_crosswalk.csv',index=False)
sp=P/'day4_2_HoloNet_Nature/full_NVU_landscape_extension/tables';old=pd.read_csv(sp/'Table_FULL00_all_NVU_candidates_HoloNet_DB_classification.csv');frames=[];inputs=[]
for name in ['Day3_1_State_LIANA_Endo_to_NVU_UC.csv','Day3_1_State_LIANA_Endo_to_NVU_NC.csv','Day3_1_LIANA_Endo_to_NVU_UC.csv','Day3_1_LIANA_Endo_to_NVU_NC.csv']:
 p=day/'tables'/name;inputs.append({'path':str(p),'exists_now':p.exists()})
 if p.exists():
  x=pd.read_csv(p);inputs[-1]['n_rows']=len(x)
  if 'communication_strength' not in x:x['communication_strength']=x['sca.LRscore'] if 'sca.LRscore' in x else x.get('connectome.weight_sc',1.)
  x['input_file']=name;frames.append(x[keys+['communication_strength','input_file']])
oldre=pd.concat(frames).drop_duplicates(keys).sort_values('communication_strength',ascending=False).head(80)
oldre['archived_top80']=pd.MultiIndex.from_frame(oldre[keys]).isin(pd.MultiIndex.from_frame(old[keys]));oldre.to_csv(out/'Spatial_top80_inputs.csv',index=False)
db=pd.read_csv(A/'HoloNet_current_human_interactions.csv',index_col=0);cx=pd.read_csv(A/'HoloNet_current_human_complexes.csv',index_col=0)
def canon(s):
 if s in cx.index:return '+'.join(sorted(str(v).upper() for v in cx.loc[s].dropna()))
 return '+'.join(sorted(x.upper() for x in re.split(r'[_+/ ()]+',str(s)) if x))
db['canonical_pair']=[canon(l)+':'+canon(r) for l,r in zip(db.ligand,db.receptor)]
db['exact_pair']=db.interaction_name_2.str.replace(' - ',':',regex=False).str.upper()
oldp=old.groupby(['ligand','receptor','LR_Pair'],as_index=False).agg(n_contexts=('source','size'),original_exact_match=('holonet_db_exact','max'))
oldp['canonical_pair']=[canon(l)+':'+canon(r) for l,r in zip(oldp.ligand,oldp.receptor)]
oldp['current_exact_match']=oldp.LR_Pair.isin(db.exact_pair);oldp['current_canonical_match']=oldp.canonical_pair.isin(db.canonical_pair)
oldp['matching_db_interactions']=oldp.canonical_pair.map(db.groupby('canonical_pair').interaction_name.agg(';'.join))
oldp['in_current_Fig6_BCD_pair']=pd.MultiIndex.from_frame(oldp[['ligand','receptor']]).isin(pd.MultiIndex.from_frame(cross[['ligand','receptor']]))
fig7=(Path('<ANALYSIS_ROOT>/fuxian/JTM_final_manuscript_figure_scripts_20260609/figure7_rebuild_nature_style_v3_final.py')).read_text(encoding='utf-8')
import ast
show13=ast.literal_eval(re.search(r'    order=(\["VWF:LRP1".*?\])',fig7)[1]);show11=ast.literal_eval(re.search(r'    keep=(\["VWF:LRP1".*?\])',fig7)[1])
oldp['in_S11_display13']=oldp.LR_Pair.isin(show13);oldp['in_Fig7D_display11']=oldp.LR_Pair.isin(show11)
oldp.to_csv(out/'Spatial28_database_crosswalk.csv',index=False)
pd.read_csv(sp/'Table_FULL03_sample_level_full_NVU_spatial_LR_scores.csv').to_csv(out/'Spatial28_context_scores.csv',index=False)
current=cross[cross.panel.isin(['6C','6D'])].drop_duplicates(keys).copy();current['canonical_pair']=[canon(l)+':'+canon(r) for l,r in zip(current.ligand,current.receptor)]
current['in_original_spatial_pool']=pd.MultiIndex.from_frame(current[['ligand','receptor']]).isin(pd.MultiIndex.from_frame(oldp[['ligand','receptor']]))
current['current_HoloNet_canonical_match']=current.canonical_pair.isin(db.canonical_pair)
current.to_csv(out/'Fig6_CD_to_spatial_crosswalk.csv',index=False)
pd.DataFrame(steps).to_csv(out/'Communication_stepwise_filter_counts.csv',index=False)
report={'filters':steps,'spatial_input_files':inputs,'spatial_top80_reproduced_keys':bool(len(oldre)==len(old)==80 and oldre.archived_top80.all()),'original_spatial_unique_pairs':len(oldp),'original_HoloNet_exact_pairs':oldp.loc[oldp.original_exact_match,'LR_Pair'].tolist(),'current_additional_canonical_matches':oldp.loc[oldp.current_canonical_match&~oldp.current_exact_match,'LR_Pair'].tolist(),'Fig6CD_unique_contexts':len(current),'Fig6CD_in_old_spatial_pairs':current.loc[current.in_original_spatial_pool,'LR_Pair'].tolist() if 'LR_Pair'in current else current.loc[current.in_original_spatial_pool,['ligand','receptor']].to_dict('records'),'display10':disp[['ligand_clean','source_families','AUPR_rank']].to_dict('records'),'top8targets':ta.head(8)[['target_clean','rank_in_vascular_receiver','logFC']].to_dict('records')}
report.update(S11_display13=show13,Fig7D_display11=show11,main_only_pairs=sorted(set(show11)-set(show13)))
(A/'communication_filter_counts.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8');print(json.dumps(report,ensure_ascii=False,indent=2))
