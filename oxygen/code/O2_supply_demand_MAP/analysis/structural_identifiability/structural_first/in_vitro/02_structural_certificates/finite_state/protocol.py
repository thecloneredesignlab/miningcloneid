"""Stage B metadata, immutable provenance and aggregate resource contracts."""
import csv
import hashlib
import os
from pathlib import Path
import tempfile

CAP=230400
RETRY_RESERVE=14400
SETTLEMENT_RESERVE=14400
MECHANISMS=('lam_max','p_mis_base','p_misseg','k_o_mis','p_wgd','buffer_smax',
 'buffer_beta','buffer_n_exp','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')
INITIALS=('init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N')
NUISANCES=INITIALS+('sigma_kary','sigma_growth')
SIF='/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif'
SIF_SHA='dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()
def read(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(path,rows,fields=None):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);rows=list(rows)
    fd,tmp=tempfile.mkstemp(prefix='.pending_',dir=p.parent)
    with os.fdopen(fd,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
        f.flush();os.fsync(f.fileno())
    os.replace(tmp,p)
def kv(path):return {r['key']:r['value'] for r in read(path)}
def roles_valid(rows):
    names=[r['parameter'] for r in rows]
    if len(names)!=20 or len(set(names))!=20 or set(names)!=set(MECHANISMS+NUISANCES):raise ValueError('Original 14+6 inventory mismatch')
    if sum(r['scientific_role']=='mechanism' for r in rows)!=14:raise ValueError('Mechanism count mismatch')
def verify_archive(path,expected):
    path=Path(path)
    assert sha(path/'collection_manifest.tsv')==expected==(path/'COMPLETE').read_text().splitlines()[0]
    rows=read(path/'collection_manifest.tsv')
    for r in rows:
        p=path/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],str(p)
    return len(rows)
def reservation(settled,live,site,used_retry=0):
    available=CAP-settled-live-SETTLEMENT_RESERVE-(0 if used_retry else RETRY_RESERVE)
    seconds=min(available,site)//60*60
    if seconds<=0:raise ValueError('No aggregate-budget-safe reservation remains')
    return seconds
def check_payload(rows,task,release):
    if not rows:raise ValueError('Incomplete result')
    identities=[]
    for r in rows:
        if r['task']!=task or r['release_sha256']!=release:raise ValueError('Stale or cross-task result')
        identities.append(r['item'])
    if len(identities)!=len(set(identities)):raise ValueError('Duplicate result identity')

def successful_worker(gates,application,policy,state,exit_code):
    return bool(gates and application and policy and state=='COMPLETED' and exit_code=='0:0'
                and gates[0].get('status')=='PASS' and application[0].get('status')=='PASS'
                and policy[0].get('status')=='PASS')
