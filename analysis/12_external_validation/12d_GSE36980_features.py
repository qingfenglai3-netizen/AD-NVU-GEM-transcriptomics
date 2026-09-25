"""GSE36980 probe-to-gene mapping and feature extraction for the fixed 18-gene module."""
from pathlib import Path
import gzip, csv, json, sys, re, io, hashlib
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
R=Path(r'<WORKDIR>'); cache=Path(r'<RAW_DATA>\GEO_bulk_screen_cache')
trainfile=Path(r'<ANALYSIS_ROOT>\fuxian_test\day6_GEM_method_enhancement_screen\tmp\bulk_GEM_gene_matrix.csv')
train=pd.read_csv(trainfile);genes=train.columns[2:].tolist()
def textopen(p):
 with open(p,'rb') as f:magic=f.read(2)
 return gzip.open(p,'rt',encoding='utf-8',errors='replace') if magic==b'\x1f\x8b' else open(p,encoding='utf-8',errors='replace')
def read_table(p,marker):
 f=textopen(p)
 for line in f:
  if line.strip()==marker: break
 d=pd.read_csv(f,sep='\t',comment='!',dtype=str);f.close();return d
plat=read_table(cache/'GPL6244.soft.gz','!platform_table_begin')
print('annotation columns',plat.columns.tolist())
print(plat.head(1).to_string(index=False)[:1300])
rows=list(csv.reader((R/'intermediate/GSE36980_headers.txt').read_text(encoding='utf-8').splitlines(),delimiter='\t'))
meta={r[0]:r[1:] for r in rows if r and r[0] in ['!Sample_title','!Sample_geo_accession','!Sample_source_name_ch1']}
dm=pd.DataFrame({'sample_id':meta['!Sample_geo_accession'],'title':meta['!Sample_title'],'source':meta['!Sample_source_name_ch1']})
dm['region']=dm.title.str.extract(r'_(FC|TC|HPC|HI)')[0].replace({'HI':'HPC'})
dm['group']=np.where(dm.title.str.startswith('AD_'),'AD','CN')
dm=dm[dm.region.eq('TC')].copy()
dm.to_csv(R/'results/external_GSE36980_manifest.csv',index=False)
expr=read_table(cache/'GSE36980_series_matrix.txt.gz','!series_matrix_table_begin').set_index('ID_REF').apply(pd.to_numeric)
symcol=next(c for c in plat if c.lower() in ['gene_assignment','gene symbol','gene_symbol'])
mapping=[]
for _,row in plat.iterrows():
 raw=str(row[symcol]);parts=raw.split(' /// ')
 symbols=set()
 for part in parts:
  ss=part.split(' // ')
  if len(ss)>1:symbols.add(ss[1].strip().upper())
 if len(symbols)==1 and next(iter(symbols)) in genes:
  mapping.append({'probe':str(row['ID']),'gene':next(iter(symbols)),'raw_annotation':raw})
mp=pd.DataFrame(mapping);mp=mp[mp.probe.isin(expr.index)]
present=[g for g in genes if g in set(mp.gene)]
pd.DataFrame({'gene':genes,'available_GSE36980':[g in present for g in genes]}).to_csv(R/'results/external_feature_coverage.csv',index=False)
mp.to_csv(R/'results/external_probe_gene_map.csv',index=False)
# Gene value: mean of all uniquely assigned probes.
X=pd.DataFrame({g:expr.loc[mp.loc[mp.gene.eq(g),'probe']].mean(axis=0) for g in present})
out=dm.join(X,on='sample_id')
out.to_csv(R/'results/external_GSE36980_features.csv',index=False)
train.to_csv(R/'results/discovery_GSE132903_features.csv',index=False)
spec={'discovery':'GSE132903','external':'GSE36980','primary_region':'TC','secondary_regions':[],'frozen_features':genes,'missing_external_features':[g for g in genes if g not in present], 'mapping':'only uniquely assigned gene symbols; arithmetic mean of all mapped probes; no diagnosis-based probe choice', 'preprocessing':'published log2 expression; discovery-fitted StandardScaler only; no test-fitted scaling or batch correction', 'model_selection':'report all 72 fixed original pipelines; primary all-gene L2 logistic; no test-based selection', 'bootstrap':'2000 stratified resamples of external subjects; percentile 95% intervals', 'caveat':'external datasets had previously been inspected for individual-gene associations; retrospective external validation, not a pristine prospective holdout', 'training_file_sha256':hashlib.sha256(trainfile.read_bytes()).hexdigest()}
(R/'intermediate/external_validation_frozen_spec.json').write_text(json.dumps(spec,indent=2),encoding='utf-8')
print('Feature coverage',len(present),'/',len(genes),'missing',spec['missing_external_features']); print(dm.groupby(['region','group']).size().to_string());print('values',X.min().min(),X.max().max())
