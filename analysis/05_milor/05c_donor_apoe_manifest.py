from pathlib import Path
import pandas as pd, csv,json,sys
sys.stdout.reconfigure(encoding='utf-8')
R=Path(r'<WORKDIR>')
m=pd.read_csv(r'<RAW_DATA>_test\D1\sample_GSM_mapping.csv')
s=pd.read_csv(r'<RAW_DATA>_test\SraRunTable.csv')
with open(r'<RAW_DATA>_test\GSE237718_series_matrix.txt',encoding='utf-8') as f:
 rows=[]
 for line in f:
  if line.startswith('!series_matrix_table_begin'):break
  rows.append(next(csv.reader([line],delimiter='\t')))
data={}
for row in rows:
 if row and row[0]=='!Sample_geo_accession':data['Library Name']=row[1:]
 if row and row[0]=='!Sample_characteristics_ch1':
  k=row[1].split(': ',1)[0];data[k]=[v.split(': ',1)[-1] for v in row[1:]]
g=pd.DataFrame(data).rename(columns={'pathological diagnosis':'pathological_diagnosis','Sex':'sex'})
assert not g['Library Name'].duplicated().any()
m=m.merge(g,left_on='GSM_Accession',right_on='Library Name',how='left',validate='1:1')
d=pd.read_csv(r'<ANALYSIS_ROOT>\fuxian_test\day1_1_merge\tables\metadata_final.csv',usecols=['orig.ident','SampleID','Dataset','Group','SubjectID']).drop_duplicates().merge(m[['SampleID','GSM_Accession','genotype','pathological_diagnosis','sex']],on='SampleID',how='left',validate='1:1')
d['APOE4_dosage']=d.genotype.str.count('E4');d['APOE2_dosage']=d.genotype.str.count('E2')
d.to_csv(R/'results/donor_APOE_manifest.csv',index=False)
print(pd.crosstab(d.genotype.fillna('Unavailable'),d.Group).to_string())
assert (d.loc[d.Dataset.eq('D1'),'Group']==d.loc[d.Dataset.eq('D1'),'pathological_diagnosis'].replace({'Control':'CN','NC':'CN','non-AD':'CN'})).all()
