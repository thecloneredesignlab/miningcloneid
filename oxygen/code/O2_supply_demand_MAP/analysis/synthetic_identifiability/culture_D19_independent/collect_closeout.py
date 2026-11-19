"""Collect and independently reconcile a terminal D19 release. No scheduler mutations."""
import csv
import hashlib
from pathlib import Path, PurePosixPath
import re
import shlex
import subprocess
import sys
import tarfile

HOST = '4482173@red.moffitt.org'
REMOTE = '/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261118_culture_D19_independent_01'
PYTHON = '/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3'
TERMINAL = {'COMPLETED', 'FAILED', 'TIMEOUT', 'CANCELLED', 'NODE_FAIL', 'OUT_OF_MEMORY',
            'PREEMPTED', 'BOOT_FAIL', 'DEADLINE', 'REVOKED', 'SPECIAL_EXIT'}
PREDICTION_TERMINALS = {'D19_PREDICTION_RECOVERY_PASS', 'D19_PREDICTION_RECOVERED_LOCAL_UNRESOLVED',
                        'DET_OPTIMIZER_UNRESOLVED'}
ORIGINS = ('D_A','D_B','D_D','D_F')
IMAGE_SHA = 'dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372'
CAPS = dict(global_=73400, local=6000, pattern=2000, joint=800, restore=800, verify=600, handoff=400)
CAPS['global'] = CAPS.pop('global_')

def ssh(command):
    return subprocess.check_output(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=12',
        '-o','ServerAliveInterval=12','-o','ServerAliveCountMax=2',HOST,command])

def in_study(command):
    return 'cd '+shlex.quote(REMOTE)+' && '+command

def table(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def truth(value):
    return str(value).lower() == 'true'

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def safe_path(value):
    path = PurePosixPath(value)
    if path.is_absolute() or '..' in path.parts or str(path) != value or not path.parts:
        raise ValueError('Unexpected archive path: '+value)
    return path

def collect(paths, dest):
    names = sorted(set(paths))
    for name in names:
        safe_path(name)
    if not names:
        return
    quoted = ' '.join(shlex.quote(path) for path in names)
    cmd = in_study('tar -cf - -- '+quoted)
    received = set()
    args = ['ssh','-o','BatchMode=yes','-o','ConnectTimeout=12','-o','ServerAliveInterval=12',
            '-o','ServerAliveCountMax=2',HOST,cmd]
    with subprocess.Popen(args, stdout=subprocess.PIPE) as proc:
        with tarfile.open(fileobj=proc.stdout, mode='r|') as archive:
            for record in archive:
                if record.name not in names or record.name in received or not record.isfile():
                    raise ValueError('Unexpected tar member: '+record.name)
                target = dest/safe_path(record.name)
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.extractfile(record) as source, target.open('wb') as sink:
                    for block in iter(lambda: source.read(1024*1024), b''):
                        sink.write(block)
                received.add(record.name)
        if proc.wait() != 0:
            raise RuntimeError('Tar transfer failed')
    if set(names) != received:
        raise RuntimeError('Incomplete archive transfer')

def verified_file(path, sealed_base, rows):
    rel = path.relative_to(sealed_base).as_posix()
    found = [r for r in rows if r['path']==rel]
    if len(found) != 1:
        raise ValueError('Missing or duplicate sealed path: '+rel)
    r = found[0]
    if path.stat().st_size != int(r['bytes']) or sha(path) != r['sha256']:
        raise ValueError('Path digest does not agree with the sealed receipt: '+rel)

def scheduler(root, jobs):
    ids = [j['job_id'] for j in jobs]
    if any(not value.isdigit() for value in ids) or len(set(ids)) != len(ids):
        raise ValueError('Missing or duplicate native allocation identifier')
    requested = ','.join(ids)
    cmd = ('sacct -nP -j '+requested+
           ' --format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS && '
           'printf "\n--- LIVE ---\n" && squeue -h -u 4482173 -o "%i %q %T"')
    fresh = ssh('bash -lc '+shlex.quote(cmd)).decode()
    (root/'accounting/independent_slurm_closeout.txt').write_text(fresh)
    records, live = fresh.split('--- LIVE ---',1)
    if any(line.split()[0] in ids for line in live.splitlines() if line.strip()):
        raise RuntimeError('A supposedly terminal release still owns live allocations')
    by_id = {}
    for line in records.splitlines():
        fields = line.split('|')
        if len(fields) >= 7:
            by_id[fields[0]] = fields
    for row in jobs:
        jid = row['job_id']
        main = by_id[jid]
        state = main[1].split()[0]
        allocated = max(int(v[3] or 0)*int(v[4] or 0) for k,v in by_id.items() if k==jid or k.startswith(jid+'.'))
        if state not in TERMINAL or state != row['state'] or main[2] != row['exit_code'] or main[6] != 'xxlarge':
            raise RuntimeError('D19 allocation state mismatch: '+jid)
        if int(row['settled_seconds']) != allocated or int(row['reserved_seconds']) != 0:
            raise RuntimeError('D19 settlement mismatch: '+jid)
    return by_id

def select_files(root):
    quality = table(root/'reports/verified_candidates.tsv')
    output = set()
    for row in quality:
        if not re.fullmatch(r'[0-9a-f]{64}', row['key']) or row['origin'] not in ORIGINS:
            raise ValueError('Unknown replay key or origin')
        for field in ('first_process_task','second_process_task'):
            task = row[field]
            if not re.fullmatch(re.escape(row['origin']) + r'_replay[12][ab]',task):
                raise ValueError('Unrecognized replay phase')
            output.add(f"verification/{task}/point_{row['key']}.rds")
            output.add(f'verification/{task}/phase_manifest.tsv')
    for origin in set(r['origin'] for r in quality):
        for suffix in ('verify1','verify2','replay2b'):
            phase = f'verification/{origin}_{suffix}'
            output.add(phase+'/phase_manifest.tsv')
            if suffix.startswith('verify'):
                output.update(phase+'/'+name for name in ('verification_polls.tsv','checked_representative.rds'))
            else:
                output.add(phase+'/state.rds')
    return output

def validate_settlement(root, jobs):
    budget, = table(root/'accounting/final_budget.tsv')
    attempts = table(root/'accounting/candidate_attempts.tsv')
    charged = sum(int(x['attempts']) for x in attempts)
    settled = sum(int(x['settled_seconds']) for x in jobs)
    if settled != int(budget['settled_core_seconds']) or settled > 432000:
        raise RuntimeError('Settled allocation use exceeds or disagrees with the hard cap')
    if len(jobs)!=int(budget['owned_jobs']) or int(budget['active_owned_jobs'])!=0:
        raise RuntimeError('Incorrect owned allocation accounting')
    if sum(j['purpose']=='retry' for j in jobs)>1 or int(budget['scientific_reserve_rounds'])>1:
        raise RuntimeError('Retry or scientific reserve authority used more than once')
    if charged>336400 or charged!=int(budget['attempts']):
        raise RuntimeError('Attempt sum exceeds or disagrees with the hard cap')
    for r in attempts:
        bound = 400 if r['origin']=='engineering' else CAPS[r['bucket']]
        if int(r['cap'])!=bound or int(r['attempts'])>bound:
            raise RuntimeError('Inconsistent origin or engineering attempt accounting')
    return settled, charged

def collect_repair(root, jobs):
    retry, = table(root/'accounting/technical_retry.tsv')
    history = retry['prior_history']
    failed = retry['failed_phase']
    patch = 'repair_launch_22050987'
    recovered = [j for j in jobs if j['task']==retry['task'] and j['attempt']=='2' and j['purpose']=='retry']
    original = [j for j in jobs if j['job_id']==retry['previous_job'] and j['state']=='FAILED']
    if len(recovered)!=1 or len(original)!=1 or int(retry['failed_scientific_calls'])!=1:
        raise RuntimeError('Technical retry evidence does not match the settled ledger')
    files = [f'{failed}/{name}' for name in ('phase_manifest.tsv','checks.tsv','console.log','evaluations.tsv',
                                             'phase_exit_codes.tsv','policy_gate.tsv')]
    files += [f'{history}/{name}' for name in ('REPORT.md','TERMINAL','final_budget.tsv','stop_record.tsv')]
    files += [f'{patch}/{name}' for name in ('manifest.tsv','launch.sh','resume.py','diagnosis.tsv')]
    files += ['controller/replay_env_repair.log',f"logs/{retry['task']}_{recovered[0]['job_id']}.txt"]
    collect(files,root)
    for rel in files:
        path = root/rel
        if path.parent==root/failed and path.name!='phase_manifest.tsv':
            verified_file(path,path.parent,table(path.parent/'phase_manifest.tsv'))
        if path.parent==root/patch and path.name!='manifest.tsv':
            verified_file(path,path.parent,table(path.parent/'manifest.tsv'))
    charged = [r for r in table(root/failed/'evaluations.tsv') if r['event']=='START']
    if len(charged)!=1 or charged[0]['bucket']!='verify':
        raise RuntimeError('Original failed verification call did not remain auditable')

def main(dest):
    root = Path(dest).resolve()
    root.mkdir(parents=True, exist_ok=True)
    terminal = ssh(in_study('cat controller/TERMINAL')).decode().strip()
    inventory = ssh(in_study('find accounting reports controller -maxdepth 1 -type f -printf "%p\\n"')).decode().splitlines()
    base = {'REPORT.md','controller/TERMINAL','source_manifest.tsv','inputs/start_registry.tsv',
            'inputs/public/D19.rds','verification/E/checks.tsv','verification/E/phase_manifest.tsv',
            'verification/E/image_sha256.txt','verification/E/forward_attempts.tsv',
            'verification/G/checks.tsv','verification/G/phase_manifest.tsv','verification/G/policy_gate.tsv'}
    base.update(p for p in inventory if p.endswith(('.tsv','.txt')) or p in ('controller/SELECTION_LOCKED','controller/controller.log'))
    collect(base,root)
    source = table(root/'source_manifest.tsv')
    verified_file(root/'inputs/public/D19.rds',root,source)
    verified_file(root/'inputs/start_registry.tsv',root,source)
    fresh = ssh(in_study(shlex.quote(PYTHON)+' -B -S release/verify_source.py .')).decode()
    if 'IMMUTABLE_SOURCE_PASS' not in fresh:
        raise RuntimeError('Source verification did not pass on HPC')
    jobs = table(root/'accounting/jobs.tsv')
    scheduler(root,jobs)
    for phase in ('E','G'):
        path = root/'verification'/phase
        for file in path.iterdir():
            if file.name!='phase_manifest.tsv':
                verified_file(file,path,table(path/'phase_manifest.tsv'))
    if (root/'verification/E/image_sha256.txt').read_text().strip()!=IMAGE_SHA:
        raise RuntimeError('Production image hash mismatch')
    positive = {r['check']:truth(r['passed']) for r in table(root/'verification/E/checks.tsv')}
    for test in ('positive_control_feasible','positive_scale_exact','positive_D_tolerance','sealed_reference_decisions'):
        if positive.get(test) is not True:
            raise RuntimeError('Isolated D19 positive control lacks a passing receipt')
    settled, attempts = validate_settlement(root,jobs)
    if (root/'accounting/technical_retry.tsv').is_file():
        collect_repair(root,jobs)
    if terminal in PREDICTION_TERMINALS:
        selected = select_files(root)
        collect(selected,root)
        for name in selected:
            path = root/name
            if path.name != 'phase_manifest.tsv':
                verified_file(path,path.parent,table(path.parent/'phase_manifest.tsv'))
        assessor = 'verification/P'
        results = ['assessment_gate.tsv','all20_parameter_recovery.tsv','exact_vector_labels.tsv',
                   'origin_search_summary.tsv','sigma_kary_kernel_diagnostics.tsv',
                   'function_error_summary.tsv','growth_prediction_comparison.tsv',
                   'full_karyotype_prediction_comparison.tsv','phase_manifest.tsv']
        collect([f'{assessor}/{f}' for f in results],root)
        manifest = table(root/assessor/'phase_manifest.tsv')
        for name in results:
            if name != 'phase_manifest.tsv':
                verified_file(root/assessor/name,root/assessor,manifest)
    print('TERMINAL_ARCHIVE_VERIFIED',terminal,settled,attempts)

if __name__=='__main__':
    main(sys.argv[1])
