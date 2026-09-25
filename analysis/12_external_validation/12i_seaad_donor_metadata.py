from pathlib import Path
import sys,json,urllib.request,hashlib
R=Path(__file__).resolve().parents[2]
import h5py,pandas as pd,numpy as np
url='https://cdn.prod.website-files.com/689cfbd308fa7373b604d290/68debdfdd1b8e9f8fd64dab0_sea-ad_cohort_donor_metadata_072524.xlsx'
p=R/'intermediate/SEAAD_official_donor_metadata.xlsx'
if not p.exists():urllib.request.urlretrieve(url,p)
xl=pd.ExcelFile(p);print(xl.sheet_names)
for sh in xl.sheet_names:
 d=pd.read_excel(p,sheet_name=sh);print(sh,d.shape,d.columns.tolist());print(d.head(3).to_string(index=False)[:4500]);d.to_csv(R/'intermediate'/f'SEAAD_metadata_{sh}.csv',index=False)
def read_col(obj):
 if isinstance(obj,h5py.Group):
  cats=obj['categories'][:];cats=[x.decode() if isinstance(x,bytes) else x for x in cats];codes=obj['codes'][:];return [cats[c] if c>=0 else np.nan for c in codes]
 a=obj[:];return [x.decode() if isinstance(x,bytes) else x for x in a]
rows=[]
for reg in ['DLPFC','MTG']:
 for ct in ['Endothelial','VLMC']:
  with h5py.File(Path(r'<RAW_DATA>/SEAAD_cellxgene_small_vascular')/f'SEAAD_{reg}_{ct}.h5ad') as f:
   keys=['donor_id','APOE4 status','Age at death','PMI','sex','disease','ADNC','Neurotypical reference'];d=pd.DataFrame({k:read_col(f['obs'][k]) for k in keys if k in f['obs']}).drop_duplicates()
   assert d.donor_id.is_unique;d['dataset']=f'{reg}_{ct}';rows.append(d)
out=pd.concat(rows);out.to_csv(R/'results/SEAAD_cached_donor_metadata.csv',index=False)
print(out.groupby(['dataset','APOE4 status','disease']).size().to_string())
(R/'intermediate/SEAAD_metadata_source.json').write_text(json.dumps({'landing_page':'https://brain-map.org/consortia/sea-ad/our-data','download':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'accessed':'2026-09-10','cached_h5ad':'Read obs donor-level fields only; expression unchanged'},indent=2))
