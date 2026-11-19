"""Copy the finite supplement, validate receipts, and leave full evidence on HPC."""
import csv
from pathlib import Path
import re
import shlex
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
import collect_closeout as transfer

transfer.REMOTE = '/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261118_culture_D19_B_lowest_local_review_01'
PREFLIGHT = ('anchor_admission.tsv','poll_coverage_expected.tsv','launcher_validation.tsv')
POLL = ('poll_coverage.tsv','poll_candidates.tsv','per_sample_tv.tsv','width_direction_diagnostics.tsv',
        'local_resolution_audit.tsv','selected_point.rds')
REPLAY = ('replayed_point.rds','selected_point_replays.tsv')
PHASE_RECEIPTS = ('phase_manifest.tsv','application_result.tsv','policy_gate.tsv','phase_exit_codes.tsv',
                  'checks.tsv','image_sha256.txt','loaded_namespaces.tsv')

def main(dest):
    root=Path(dest).resolve()
    root.mkdir(parents=True,exist_ok=True)
    paths=['accounting/jobs.tsv','accounting/slurm_snapshot.txt','source_manifest.tsv',
           'release/source_and_anchor_manifest.tsv','inputs/public/D19.rds',
           'verification/D_B_anchor/replay_a.rds','verification/D_B_anchor/replay_b.rds']
    transfer.collect(paths,root)
    source=transfer.table(root/'source_manifest.tsv')
    for name in paths[3:]:
        transfer.verified_file(root/name,root,source)
    jobs=transfer.table(root/'accounting/jobs.tsv')
    stages={}
    ids=[]
    for r in jobs:
        stage=r['stage']
        stages[stage]=r
        jid=r['job_id']
        if not jid.isdigit() or jid in ids:
            raise RuntimeError('Invalid supplement job identity')
        ids.append(jid)
    if any((r['state'],r['exit_code'])!=('COMPLETED','0:0') for r in stages.values()):
        raise RuntimeError('Incomplete latest supplement allocation')
    if not {'admit','poll'}<=set(stages) or not set(stages)<={'admit','poll','replayA','replayB'}:
        raise RuntimeError('Invalid stage set')
    if sum(int(r['settled_seconds']) for r in jobs)>7200 or sum(int(r['charged_calls']) for r in jobs)>155:
        raise RuntimeError('Independent budget overflow')
    snapshot=transfer.ssh('bash -lc '+shlex.quote(
        'sacct -nP -j '+','.join(ids)+' --format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,QOS && '
        'printf "\\n--- LIVE ---\\n" && squeue -h -u 4482173 -o "%i"')).decode()
    (root/'accounting/independent_slurm_closeout.txt').write_text(snapshot)
    summary,live=snapshot.split('--- LIVE ---',1)
    if any(line.strip() in ids for line in live.splitlines()):
        raise RuntimeError('Owned supplement job remains active')
    records={v[0]:v for line in summary.splitlines() if len(v:=line.split('|'))>=6}
    for r in jobs:
        job=r['job_id']
        v=records[job]
        charged=max(int(x[3] or 0)*int(x[4] or 0) for k,x in records.items() if k==job or k.startswith(job+'.'))
        if v[1].split()[0]!=r['state'] or v[2]!=r['exit_code'] or v[5]!='xxlarge' or charged!=int(r['settled_seconds']):
            raise RuntimeError('Native scheduler reconciliation failed: '+job)
    files = dict(admit=PREFLIGHT,poll=POLL,replayA=REPLAY,replayB=REPLAY)
    receipt_root={}
    for stage,r in stages.items():
        task=r['task']
        if not re.fullmatch(r'D_B_(admit|poll|replayA|replayB)(_retry)?',task):
            raise RuntimeError('Invalid verified phase label')
        folder='verification/'+task
        receipt_root[stage]=root/folder
        transfer.collect([folder+'/'+f for f in PHASE_RECEIPTS+files[stage]],root)
        manifest=transfer.table(root/folder/'phase_manifest.tsv')
        for name in PHASE_RECEIPTS+files[stage]:
            if name!='phase_manifest.tsv':
                transfer.verified_file(root/folder/name,root/folder,manifest)
        if (root/folder/'image_sha256.txt').read_text().strip()!=transfer.IMAGE_SHA:
            raise RuntimeError('Production image mismatch')
        checks=transfer.table(root/folder/'checks.tsv')
        if not checks or any(r['passed']!='TRUE' for r in checks):
            raise RuntimeError('Scientific test failed in '+stage)
        policy,=transfer.table(root/folder/'policy_gate.tsv')
        if policy['status']!='PASS':
            raise RuntimeError('Scientific policy audit failed in '+stage)
    d,=transfer.table(receipt_root['poll']/'local_resolution_audit.tsv')
    required={'admit','poll','replayA','replayB'} if d['lower_found']=='TRUE' else {'admit','poll'}
    if required!=set(stages): raise RuntimeError('Replay allocation count conflicts with locked scientific result')
    print('SUPPLEMENT_ARCHIVE_VERIFIED',len(jobs),sum(int(r['charged_calls']) for r in jobs))

if __name__=='__main__':
    main(sys.argv[1])
