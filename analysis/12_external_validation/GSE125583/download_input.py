from pathlib import Path
import hashlib
from urllib.request import urlretrieve
p=Path(__file__).resolve().parent/'GSE125583_RAW.tar'
if not p.exists():
    urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125583/suppl/GSE125583_RAW.tar', p)
assert hashlib.sha256(p.read_bytes()).hexdigest() == 'f13f8b565230403684a9ff0658b38ac0dd49d682df949b5cd9fceb6a1d88e3ec', 'Input checksum mismatch'
print('GSE125583 expression archive verified.')
