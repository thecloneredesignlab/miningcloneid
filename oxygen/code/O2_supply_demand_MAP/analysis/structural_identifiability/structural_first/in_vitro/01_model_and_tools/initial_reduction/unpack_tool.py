from pathlib import Path
import hashlib,sys,tarfile
s=Path(sys.argv[1]);scratch=Path(sys.argv[2]);archive=s/'toolchain/verified_depot.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()==(s/'toolchain/verified_depot.sha256').read_text().strip()
with tarfile.open(archive) as tar:
 for x in tar.getmembers():assert not Path(x.name).is_absolute() and '..' not in Path(x.name).parts
 tar.extractall(scratch)
print('Verified frozen depot extracted to node-local scratch',flush=True)
