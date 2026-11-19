"""Freeze the authorized, isolated D19 restart. Metadata/file copies only."""
from pathlib import Path
import os,shutil,subprocess,sys
from common import *

def prepare(repo):
 repo=Path(repo).resolve();s=repo/'oxygen/results/synthetic_identifiability/20261117_culture_D19_recovery_restart_01'
 assert not s.exists(),'Never overwrite an existing study'
 f=repo/'oxygen/results/synthetic_identifiability/20261117_culture_D19_feasibility_01'
 old=repo/'oxygen/results/synthetic_identifiability/20261116_culture_growth_karyotype_recovery_01'
 code=Path(__file__).resolve().parent
 for d in ['inputs/public','inputs/private','inputs/banks','inputs/rng','inputs/production','verification','accounting','controller','logs','reports','release','runtime']:(s/d).mkdir(parents=True)
 expected=dict(A='8e3ebac616f384f4bc20b9c1ce3c3db800d81efc2f13078684330678cefa2115',B='e27665ffc68ae7a5a43dbed163e99fbba969b6d9be819d63189ea9e7582f1e44',C='765a740842c7eab130ae7ace9305fafbc2fd03c2ce9033d0637bd6318e36c4ab',D='f4e29c3532af47af3afc43045a6d3981e34a563e81d4452be259182d60b647b7')
 evidence=[]
 def cp(src,dst,expected_sha=None):
  h=sha(src)
  if expected_sha:assert h==expected_sha,'Source hash conflict '+str(src)
  q=s/dst;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,q)
  evidence.append(dict(source=str(src.relative_to(repo)),destination=dst,bytes=q.stat().st_size,sha256=h))
 for letter,h in expected.items():cp(f/f'verification/F_{letter}_de4/bank.rds',f'inputs/banks/D_{letter}.rds',h)
 cp(old/'verification/I01/cases/D19.rds','inputs/public/D19.rds','0b1787b85af14c024c420c281c59c50ef92ac5a4f10ada926f1cc614ced1cbfe')
 cp(f/'verification/P0_R1/feasibility_contract.rds','inputs/public/contract.rds')
 for name in ['protocol_edges.tsv','observation_days.tsv','exact_bound_exports.tsv']:cp(f/'verification/P0_R1'/name,'inputs/public/'+name)
 cp(f/'inputs/parameter_roles.tsv','inputs/public/parameter_roles.tsv')
 cp(f/'inputs/private/truth.rds','inputs/private/truth.rds')
 cp(f/'inputs/private/reference.rds','inputs/private/reference.rds')
 for p in sorted((f/'inputs/production').iterdir()):cp(p,'inputs/production/'+p.name)
 shutil.copytree(f/'runtime',s/'runtime',dirs_exist_ok=True)
 for p in sorted(code.iterdir()):
  if p.is_file():shutil.copy2(p,s/'release'/p.name)
 # Non-scientific negative-read fixtures exist before deployment; filenames matter.
 for name in ['fit_data.Rds','jobs_2N.Rds','jobs_4N.Rds']:(s/'inputs/private'/name).write_bytes(b'NON_SCIENTIFIC_DENIAL_FIXTURE\n')
 (s/'inputs/private/private_fixture.txt').write_text('No scientific coordinates in this fixture.\n')
 (s/'inputs/private/engineering_fixture.txt').write_text('Privileged engineering output fixture.\n')
 (s/'inputs/public/alias_private.txt').symlink_to('../private/private_fixture.txt')
 (s/'inputs/public/alias_target.rds').symlink_to('D19.rds')
 (s/'verification/D_B_fixture').mkdir();(s/'verification/D_B_fixture/state.rds').write_bytes(b'FOREIGN_ORIGIN_FIXTURE\n')
 (s/'inputs/public/alias_foreign.rds').symlink_to('../../verification/D_B_fixture/state.rds')
 plan=repo/'oxygen/results/analysis_documentation/in_vitro/synthetic_recovery/culture_D19_deterministic_recovery_restart_plan_20261117.md'
 cp(plan,'inputs/private/confirmed_plan.md')
 table(s/'inputs/evidence_registry.tsv',evidence)
 table(s/'accounting/authority.tsv',[dict(scope='D19 restart only; four primaries; one E/F round; no N20, profiles, tumor',core_seconds=CAP,technical_retries=1,scientific_reserve_rounds=1,mechanisms=14,nuisances=6,active=19,known_control='sigma_growth',full_forward_external_cap=512)])
 table(s/'accounting/retry_ledger.tsv',[dict(key='used',value=0)])
 table(s/'accounting/jobs.tsv',[],['job_id','task','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
 table(s/'accounting/purpose_envelopes.tsv',[dict(purpose=k,planned_core_seconds=v) for k,v in [('engineering',28800),('D_A',100800),('D_B',100800),('D_C',100800),('D_D',100800),('reserve_E',108000),('reserve_F',108000),('assessment',28800),('collection',14400),('technical_retry',28800)]])
 # Freeze pre-existing tracked dirty content; no automatic Git operations.
 dirty=subprocess.check_output(['git','diff','--name-only','-z'],cwd=repo).split(b'\0')
 table(s/'accounting/preservation.tsv',[dict(path=os.fsdecode(p),sha256=sha(repo/os.fsdecode(p))) for p in dirty if p and (repo/os.fsdecode(p)).is_file()])
 (s/'accounting/preexisting_status.txt').write_text(subprocess.check_output(['git','status','--short','--branch'],cwd=repo,text=True))
 table(s/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=43200,checked='2026-11-17; refresh remotely before submission')])
 paths=[p for sub in ['release','runtime','inputs'] for p in sorted((s/sub).rglob('*')) if p.is_file() and not p.is_symlink()]
 table(s/'source_manifest.tsv',[dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p)) for p in paths])
 verify(s);print(s)
if __name__=='__main__':prepare(sys.argv[1])
