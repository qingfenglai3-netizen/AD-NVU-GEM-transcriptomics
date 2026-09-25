from pathlib import Path
import csv,re,json
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'intermediate/GSE118553';rows=list(csv.reader((O/'headers.txt').read_text(encoding='utf-8').splitlines(),delimiter='\t'))
get=lambda key:next(r[1:] for r in rows if r and r[0]==key)
d=pd.DataFrame({'sample_id':get('!Sample_geo_accession'),'title':get('!Sample_title')})
for row in rows:
 if row and row[0]=='!Sample_characteristics_ch1':
  keys=[x.split(': ',1)[0] for x in row[1:]];assert len(set(keys))==1
  d[keys[0]]=[x.split(': ',1)[1] for x in row[1:]]
d['donor_id']=d.individual.str.replace(r'\s+',' ',regex=True).str.strip()
t=d[d.tissue.eq('Temporal_Cortex')].copy();t.to_csv(R/'results/GSE118553_TC_metadata_all.csv',index=False)
print('Temporal samples',len(t),'donors',t.donor_id.nunique());print(t.groupby('disease state').agg(samples=('sample_id','size'),donors=('donor_id','nunique')))
print('Repeated',t[t.duplicated('donor_id',keep=False)].to_string(index=False));print('labels',t['disease state'].unique())
print('age',t.groupby('disease state').age.agg(['min','max']))
assert t.sample_id.is_unique and t.groupby('donor_id')['disease state'].nunique().max()==1
