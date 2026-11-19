"""Standard-library metadata and finite authorization contracts."""
import csv
import hashlib
import os
from pathlib import Path
import tempfile

CAP=720000
RETRY=28800
SETTLEMENT=14400
SIF='/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif'
SIF_SHA='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'
MECHANISMS=('lam_max','p_mis_base','p_misseg','k_o_mis','p_wgd','buffer_smax','buffer_beta','buffer_n_exp','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')
NUISANCES=('init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N','sigma_kary','sigma_growth')
def read(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(path,rows,fields=None):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);rows=list(rows)
    fd,tmp=tempfile.mkstemp(prefix='.pending_',dir=p.parent)
    with os.fdopen(fd,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()
def kv(path):return {r['key']:r['value'] for r in read(path)}
def roles_valid(rows):
    names=[r['parameter'] for r in rows]
    if len(names)!=20 or len(set(names))!=20 or set(names)!=set(MECHANISMS+NUISANCES):raise ValueError('Original20-role identity mismatch')
    if sum(r['scientific_role']=='mechanism' for r in rows)!=14:raise ValueError('Mechanism count mismatch')
def verify(study):
    for r in read(Path(study)/'source_manifest.tsv'):
        p=Path(study)/r['path']
        if p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:raise ValueError('Immutable source mismatch: '+r['path'])
