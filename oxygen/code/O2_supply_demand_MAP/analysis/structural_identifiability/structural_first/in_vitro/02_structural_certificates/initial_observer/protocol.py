"""Standard-library metadata contracts; no scientific imports or timers."""
import csv
import hashlib
from pathlib import Path

MECHANISMS = ('lam_max','p_mis_base','p_misseg','k_o_mis','p_wgd',
              'buffer_smax','buffer_beta','buffer_n_exp','alpha_o2','gamma_growth',
              'mu_hp','gamma_mu','O2_crit','n_O')
NUISANCES = ('init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N',
             'sigma_kary','sigma_growth')
PARENTS = {
 '20261116_culture_step04_public_tools_01':'6ea489717a56b25bb1932635b91f8ff8de91a71036bfe0184ea32203088d00ef',
 '20261116_culture_step04_four_controls_01':'0ff28dac7b66fbec51883a52fc20f8944d7b860e58b0882d4880797da1e61e06',
 '20261116_culture_step04_saved_review_01':'222cfface062a37a3a8695a8f51fc2697f34475ce266b8ca05dc4ef4ca03896d',
 '20261116_culture_step04_C04_certificate_01':'6b55fa7b944ef36dc2ffe4eefc302b7df08d2f3f6732018c2e3cb9b2cc28f706'}
CAP = 28800
SETTLEMENT_RESERVE = 3600
RETRY_RESERVE = 3600

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    with Path(path).open() as f: return list(csv.DictReader(f,delimiter='\t'))
def table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    rows=list(rows)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)

def validate_roles(rows):
    names=[r['parameter'] for r in rows]
    if len(names)!=20 or len(set(names))!=20 or set(names)!=set(MECHANISMS+NUISANCES):
        raise ValueError('Original twenty-coordinate inventory mismatch')
    if sum(r['scientific_role']=='mechanism' for r in rows)!=14:
        raise ValueError('Mechanism count mismatch')
    if any(r['structural_result']!='not_analyzed' for r in rows):
        raise ValueError('Original culture classification is outside this package')

def validate_schema(columns,required,cohort_map=None,assay_widths=None,diagnostic=False):
    if len(columns)!=len(set(columns)) or set(columns)!=set(required):
        raise ValueError('Missing, extra or duplicated unknown coordinate')
    if 'sigma_growth' in columns:
        raise ValueError('Variance-only coordinate is outside the mean-column schema')
    if cohort_map:
        for cohort in ('2N','4N'):
            values=[tuple(p) for c,e,p in cohort_map if c==cohort]
            if len(values)!=2 or len(set(values))!=1:
                raise ValueError('Initial-family sharing is not preserved')
        a=next(tuple(p) for c,e,p in cohort_map if c=='2N')
        b=next(tuple(p) for c,e,p in cohort_map if c=='4N')
        if set(a)&set(b): raise ValueError('Distinct cohort initial families collapsed')
    if assay_widths and not diagnostic and len(set(assay_widths))!=1:
        raise ValueError('Shared assay width duplicated between experiments')

def validate_collection(rows,expected_case,expected_release,expected_input):
    if not rows: raise ValueError('Missing matrix records')
    identities=[]
    for r in rows:
        if r['case']!=expected_case or r['release_sha256']!=expected_release or r['input_sha256']!=expected_input:
            raise ValueError('Stale or cross-case record')
        identities.append((r['case'],r['order'],r['path']))
    if len(identities)!=len(set(identities)):
        raise ValueError('Duplicate task/matrix identity')

def allocation_seconds(settled,live=0,retry=False,site_seconds=43200):
    if settled<0 or live<0 or site_seconds<60: raise ValueError('Invalid accounting')
    available=CAP-settled-live-SETTLEMENT_RESERVE-(0 if retry else RETRY_RESERVE)
    seconds=min(available,site_seconds)//60*60
    if seconds<=0: raise ValueError('No budget-safe allocation remains')
    return seconds

def verify_archive(study,expected=None):
    study=Path(study)
    digest=sha(study/'collection_manifest.tsv')
    if (study/'COMPLETE').read_text().splitlines()[0]!=digest or (expected and digest!=expected):
        raise ValueError('Closed manifest mismatch')
    rows=read(study/'collection_manifest.tsv')
    for r in rows:
        p=study/r['path']
        if p.stat().st_size!=int(r['bytes']) or sha(p)!=r['sha256']:
            raise ValueError('Closed file mismatch: '+str(p))
    return rows
