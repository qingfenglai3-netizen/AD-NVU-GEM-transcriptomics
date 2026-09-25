"""Analysis plan and download of GSE118553."""
from pathlib import Path
import json,urllib.request,datetime,hashlib,gzip,csv
R=Path(__file__).resolve().parents[2];O=R/'intermediate/GSE118553';O.mkdir(exist_ok=True)
spec={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'training':'GSE132903 label-blind gene-wise probe means,195 subjects','external':'GSE118553','region':'Temporal_Cortex only','contrast':'AD vs Control; AsymAD excluded a priori','unit':'individual donor, average repeated technical samples within donor and region before prediction; require consistent diagnosis','features':'same fixed18 genes, no external selection','probe_mapping':'same GPL10558 annotation and all corresponding probe means, without labels; all18 required; no imputation','preprocessing':'published log2 values; training-only StandardScaler; no joint or external normalization/recalibration','models':'same72 fixed pipelines sklearn1.5.2; fixed primary all18 Logistic_L2(liblinear,C1); no outcome-based best-model choice','metric':'donor-level AUC; primary percentile95%CI from2000 stratified bootstrap resamples seed20260910','stop_rules':'do not fit if donor labels/identity cannot be resolved or fewer than10 AD or10controls or any fixed feature missing'}
f=O/'frozen_spec.json'
if not f.exists():f.write_text(json.dumps(spec,indent=2),encoding='utf-8')
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118553/matrix/GSE118553_series_matrix.txt.gz';p=O/'GSE118553_series_matrix.txt.gz'
if not p.exists():urllib.request.urlretrieve(url,p)
with gzip.open(p,'rt',encoding='utf-8') as fh:
 lines=[]
 for line in fh:
  if line.startswith('!series_matrix_table_begin'):break
  lines.append(line)
(O/'headers.txt').write_text(''.join(lines),encoding='utf-8')
(O/'source.json').write_text(json.dumps({'url':url,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size},indent=2))
for line in lines:
 if line.startswith(('!Series_summary','!Series_overall_design','!Series_pubmed_id')):print(line[:3000])
rows=list(csv.reader(lines,delimiter='\t'));print([(r[0],r[1:4]) for r in rows if r and r[0].startswith('!Sample_characteristics')])
