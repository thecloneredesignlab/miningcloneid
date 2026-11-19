"""Install a checksum-verified official Julia binary in this study only."""
from pathlib import Path
import hashlib,sys,tarfile,urllib.request
s=Path(sys.argv[1]);t=s/'toolchain';t.mkdir(exist_ok=True)
name='julia-1.10.13-linux-x86_64.tar.gz'
checksum=urllib.request.urlopen('https://julialang-s3.julialang.org/bin/checksums/julia-1.10.13.sha256',timeout=60).read().decode()
(t/'official_julia_checksums.txt').write_text(checksum)
expected=[line.split()[0] for line in checksum.splitlines() if line.split()[-1].lstrip('*')==name]
assert len(expected)==1,'Official checksum not unique'
archive=t/name
if not archive.exists():
 with urllib.request.urlopen('https://julialang-s3.julialang.org/bin/linux/x64/1.10/'+name,timeout=60) as source,archive.open('wb') as out:
  while True:
   block=source.read(1048576)
   if not block:break
   out.write(block)
digest=hashlib.sha256(archive.read_bytes()).hexdigest();assert digest==expected[0],'Julia checksum mismatch'
(t/'julia_archive.sha256').write_text(digest+'  '+name+'\n')
if not (t/'julia-1.10.13/bin/julia').exists():
 with tarfile.open(archive) as pack:
  for member in pack.getmembers():
   assert not Path(member.name).is_absolute() and '..' not in Path(member.name).parts
  pack.extractall(t)
(t/'environment').mkdir(exist_ok=True)
print('Verified Julia 1.10.13',digest,flush=True)
