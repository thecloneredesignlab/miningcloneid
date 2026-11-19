"""Metadata-only final collection; never imports a simulator or optimizer."""
from collections import Counter
from pathlib import Path
import shutil
import sys
from common import *

def verify_manifest(folder,name):
 rows=read(folder/name)
 for r in rows:
  p=folder/r['path'];assert p.is_file() and sha(p)==r['sha256']
  if 'bytes' in r:assert p.stat().st_size==int(r['bytes'])
 return len(rows)

def main():
 project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();assert not (study/'CLOSED').exists()
 verify(study);verify_manifest(study,'repair_manifest.tsv');roles=read(study/'inputs/parameter_roles.tsv');roles_valid(roles)
 tasks=['I01_R1','P01','G01','R01','R02','Q01','Q02','D01','D02','B01'];jobs=read(study/'accounting/jobs.tsv')
 assert len(jobs)==11 and {r['task'] for r in jobs}==set(tasks)|{'I01'}
 failed=next(r for r in jobs if r['task']=='I01');assert failed['state']=='FAILED' and failed['exit_code']=='1:0'
 verify_manifest(study/'verification/I01','phase_manifest.tsv')
 assert not (study/'verification/I01/forward_ledger.tsv').exists() and not (study/'verification/I01/observation_attempts.tsv').exists()
 accepted_rows=[]
 for task in tasks:
  j=next(r for r in jobs if r['task']==task);out=study/'verification'/task
  assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
  assert read(out/'application_result.tsv')[0]['status']=='PASS'
  assert read(out/'policy_gate.tsv')[0]['status']=='PASS'
  assert read(out/'gate.tsv')[0]['status'].startswith('PASS')
  assert all(r['passed']=='TRUE' for r in read(out/'checks.tsv')) if task!='G01' else True
  n=verify_manifest(out,'phase_manifest.tsv')
  accepted_rows.append(dict(task=task,job_id=j['job_id'],slurm_state=j['state'],exit_code=j['exit_code'],application='PASS',policy='PASS',files=n,charged_core_seconds=j['actual_core_seconds']))
 spent=sum(int(j['actual_core_seconds']) for j in jobs);assert spent<=14400 and not any(int(j['reserved_core_seconds']) for j in jobs)
 assert int(kv(study/'accounting/retry_ledger.tsv')['used'])<=1
 old=project/'oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01'
 assert sha(old/'collection_manifest.tsv')=='48f6549ccd2094b2794b6b579b6caedb12f0a90cf5c89597433fe19b4dcc11cc'
 assert verify_manifest(old,'collection_manifest.tsv')==807
 prior=[]
 for r in read(study/'preservation_baseline.tsv'):
  p=project/r['project_path'];assert p.is_file() and sha(p)==r['sha256'],r['project_path'];prior.append(dict(**r,status='UNCHANGED'))
 report=study/'reports';table(report/'phase_acceptance.tsv',accepted_rows);table(report/'preservation_checks.tsv',prior)
 forward=[dict(task=t,**r) for t in ('R01','R02') for r in read(study/'verification'/t/'forward_ledger.tsv')]
 attempts=[dict(task=t,**r) for t in ('Q01','Q02','D01','D02') for r in read(study/'verification'/t/'observation_attempts.tsv')]
 fixtures=read(study/'verification/Q01/fixture_attempts.tsv')
 assert len(forward)<=4 and all(r['status']=='returned' for r in forward)
 assert len(attempts)<=1027 and all(r['status']=='returned' for r in attempts)
 assert len(fixtures)<=12 and all(r['status']=='returned' for r in fixtures)
 qa=[r for r in attempts if r['task'] in ('Q01','Q02')];assert len(qa)==1024 and len({r['id'] for r in qa})==1024
 qa_checks=read(study/'verification/Q02/generator_QA.tsv');assert len(qa_checks)==1824 and all(r['passed']=='TRUE' for r in qa_checks)
 first=read(study/'verification/D01/dataset_receipt.tsv')[0];second=read(study/'verification/D02/dataset_receipt.tsv')[0]
 assert first['payload_sha']==second['payload_sha'] and first['growth']=='114' and first['kary_cells']=='220'
 assert sha(study/'verification/D01/synthetic_dataset.rds')==sha(study/'verification/D02/synthetic_dataset.rds')
 geometry=read(study/'verification/G01/local_observation_geometry.tsv');assert len(geometry)==19
 assert next(r for r in geometry if r['parameter']=='sigma_kary')['status']=='NUMERICALLY_FLAT_AT_REGISTERED_PROBES'
 public=study/'public';public.mkdir()
 for src,dst in [('verification/D01/evaluation_bundle.rds','evaluation_bundle.rds'),('verification/D01/synthetic_dataset.rds','synthetic_dataset.rds'),('verification/D01/synthetic_growth.tsv','synthetic_growth.tsv'),('verification/D01/synthetic_karyotype_cells.tsv','synthetic_karyotype_cells.tsv'),('verification/D01/synthetic_observation_manifest.tsv','synthetic_observation_manifest.tsv'),('verification/I01_R1/public_bounds.tsv','public_bounds.tsv')]:shutil.copy2(study/src,public/dst)
 (public/'README.md').write_text('''# Conditional artificial dataset and future evaluation interface

One fixed-input production-compatible artificial dataset:114 endpoint-Gaussian growth responses and12 integer karyotype samples/220 cells. There are no synthetic raw counts, flow or live/dead observations. Historical counts/history and data-derived bounds are disclosed external algorithmic information; this is not a physical OLS experiment. Unknown primary coordinates remain14 mechanisms plus6 nuisances. Sigma_kary is numerically flat at the registered seed144 probes; kernel/distribution targets must be reported separately from raw width. This bundle contains no generating parameter object or latent trajectories. A scoped native-guard test denied reading the private truth object. No fit/recovery claim or fitting authorization is included. Do not read private generation artifacts when building a future blind fitting bundle. Provenance/constraints are in the parent report.
''')
 table(public/'manifest.tsv',[dict(path=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(public.iterdir()) if p.is_file() and p.name!='manifest.tsv'])
 table(report/'forward_attempts.tsv',forward);table(report/'observation_attempts.tsv',attempts)
 for src,dst in [('verification/I01_R1/conditioning_manifest.tsv','conditioning_manifest.tsv'),('verification/I01_R1/nuisance_configurations.tsv','nuisance_configurations.tsv'),('verification/G01/local_observation_geometry.tsv','growth_karyotype_geometry.tsv'),('verification/G01/singular_values.tsv','mask_singular_values.tsv'),('verification/Q02/generator_QA.tsv','generator_QA.tsv')]:shutil.copy2(study/src,report/dst)
 table(report/'parameter_roles.tsv',[dict(**r,future_primary='unknown; retained diagnostic if unsupported',current_fit_status='NOT_STARTED',sigma_policy='kernel/distribution target at flat seed144' if r['parameter']=='sigma_kary' else 'positive variance law; not mean rank' if r['parameter']=='sigma_growth' else 'retain original activity and bounds') for r in roles])
 table(report/'release_gates.tsv',[dict(scope=s,status=v) for s,v in [('conditional_generator','CONDITIONAL_GENERATOR_VALIDATED_SINGLE_DATASET_FROZEN'),('sparse_discrete_identification','NOT_PROVED'),('physical_experiment_design','HELD'),('parameter_recovery','NOT_TESTED'),('fitting_profiles_tumor','NOT_STARTED_OR_AUTHORIZED')]])
 table(study/'accounting/resource_summary.tsv',[dict(key=k,value=v) for k,v in dict(cap_core_seconds=14400,settled_core_seconds=spent,retries_used=kv(study/'accounting/retry_ledger.tsv')['used'],forward_attempts=len(forward),fixture_passes=len(fixtures),QA_realizations=1024,total_observation_attempts=len(attempts),distinct_primary_datasets=1,fits=0,profiles=0,tumor=0).items()])
 counts=Counter(r['status'] for r in geometry);sv=read(study/'verification/G01/singular_values.tsv');ranks={r['scaling']:(r['nominal_rank_1e6'],r['nominal_rank_1e8'],r['nominal_rank_1e10']) for r in sv}
 replays=[read(study/'verification'/t/'replay_checks.tsv')[0] for t in ('R01','R02')]
 (study/'REPORT.md').write_text(f'''# Conditional growth+karyotype generator execution

Status: CONDITIONAL_GENERATOR_VALIDATED_SINGLE_DATASET_FROZEN.
Sparse discrete structural identification: NOT_PROVED. Physical-study design: HELD.
Parameter recovery: NOT_TESTED. Fits/profiles/tumor: NOT_STARTED.

## Evidence and implementation

The confirmed conditional case retains14 mechanisms and6 nuisances,133 live
states22:154, production v2, ploidy_related death, original bounds/zero mappings,
root families, lineage inheritance, count-target selection, independent parent
boundary choice and deterministic propagation. Counts/history remain external
algorithmic inputs. Artificial Gaussian endpoint growth is not raw-count OLS
or its true conditional law given historical counts. Flow and death do not enter
the new score, generation, denominators or public response fields. Remaining
per-passage and per-sample averaging/weights are unchanged. All20 explicit
candidate coordinates are required; sigma_growth is not read from fixed nuisance
truth in the unknown-nuisance interface. Known-nuisance configurations are
registered controls only. None was fitted.

Rehashed closed Step5 evidence (807 files) unchanged. Its three ideal19-rank
certificates,12 unresolved old direct runs and broader HOLD retain their scope.
The new1698-row19-column saved-mask numerical geometry is a fresh projection,
with karyotype-only floor flags, not copied full-mask ranks. Status counts:
{dict(counts)}. Nominal ranks1e-6/1e-8/1e-10: {ranks}. These diagnostics do not
prove structural identification or parameter compensation. No structural case
or saved perturbation was repropagated.

Two new full two-background/history replays passed original saved-state, growth,
karyotype, selection, parent-boundary and independently scored new-mask checks.
Receipts: {replays}. The original biological functions and C++ are unchanged;
private compilation bypasses only the old timed-lock bootstrap. No optimizer,
old controller, old mutable counter or deadline was reused.

## Observation validation and frozen data

Registered RNG: R4.4.2, L'Ecuyer-CMRG, Inversion normals, Rejection sampling,
master seed20261116;1024 QA streams and one distinct primary stream. Per-sample
substreams are canonical and separate. Sampling uses the normalized finite
integer kernel before the1e-12 score floor. No invalid vector was repaired with
uniform probabilities or an unbounded rounded Gaussian. No process noise,
random bottleneck, oxygen observation or synthetic raw count was introduced.

All1824 predeclared simultaneous marginal QA bands passed on1024 observation-
only realizations from the same latent bundle, with total assumed-law false-
alarm allowance0.01. This is generator QA, not a fitted replication study,
parameter interval, physical-model validation or power calculation. No QA seed
was selected or redrawn for a favorable result. Exactly one primary artificial
dataset was frozen:114 growth responses and12 samples/220 integer chromosome
observations. Negative growth draws retained:{first['negative_growth']}.
Independent-process regeneration gave the same payload and RDS file hash.
Payload SHA:{first['payload_sha']}.

Sigma_kary remains the original small seed144 value. Its four registered local
probes remain numerically flat; a declared width1 interface control responds.
The control did not replace the generating truth. More finite cell noise does
not recover information absent from those evaluated kernels. Future reporting
must separate raw-width and kernel/distribution recovery and retain the unknown
coordinate honestly. Sigma_growth is a positive variance-law target, outside
the deterministic mean map; it was not fitted on exact means.

The public evaluation bundle has conditional design, public bounds/fixed
configuration and synthetic responses, with no generating active parameter
object or latent states. Active defaults are independent neutral values and
complete injection is mandatory. The native read guard denied the registered
private-truth canary read in the scoped future-evaluator test. This is a tested
interface boundary, not an unrestricted hostile-code security claim. Historical
bounds/counts remain disclosed information. Karyotype physical reassignment,
count-error/SOP, flow and live/dead contracts remain unresolved for broader work.

## Resources and policy

Ten accepted xxlarge allocations completed0:0; the original I01 interface failure
and its charged allocation remain preserved. The single authorized repair reused
the original endpoint-roundoff helper before any propagation or generation; application, source/hash and scientific
worker/child dependency/write audits passed independently. Charged{spent}
allocated core-seconds ({spent/3600:.6f} core-hours) of4 authorized. Retry ledger
used:{kv(study/'accounting/retry_ledger.tsv')['used']}/1. Forward attempts:{len(forward)}/4;
standalone fixture passes:{len(fixtures)}/12; observation attempts:{len(attempts)}/1027.
No application concurrency cap or discretionary research deadline was imposed;
Slurm reservations derive only from aggregate budget and unavoidable site limits.
No own allocation remains running/queued. Raw MaxRSS is scheduler sampled.

English files; Python/R/C++/shell only, no Java/JavaScript/JSON workflow/files.
Production/model/sampling calculations ran inside the pinned SIF on xxlarge;
saved numerical SVD used the recorded guarded host runtime on xxlarge. No fit,
profile, tumor analysis, author contact, old-task cancellation/restart or real-
data global-optimum search occurred. Old evidence and unrelated work were
preserved. Public bundle: public/. Private generator artifacts remain outside
that bundle. Unused budget authorizes no new work.

A later deterministic/noisy recovery proposal needs separate authorization,
optimizer validation and measured fitting budget. Generator success cannot be
reported as20-coordinate recovery or a closed structural/practical gate.
''')
 (report/'next_stage_requirements.md').write_text('''# Next-stage boundary

No fitting is released. Review the conditional law and new-mask diagnostics, the sigma_kary flatness/function estimands, all20 unknown roles and public bundle isolation. Before any deterministic/noisy fit, separately authorize the estimator interface, blind-start streams, finite evaluation/search ceilings and a measured fitting budget. Do not demand noiseless free-sigma_growth recovery, select truth by favorable recovery, replace known-nuisance controls for Nall or apply uncalibrated likelihood-ratio thresholds. Physical OLS/count/SOP, sample provenance, flow and live/dead contracts remain held. No new global-optimum search is prerequisite.
''')
 code=Path(__file__).parent
 table(report/'final_implementation_manifest.tsv',[dict(project_path=str(p.relative_to(project)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(code.iterdir()) if p.is_file()])
 assert not [p for p in study.rglob('*') if p.is_file() and (p.suffix.lower() in ('.json','.jsonl','.js','.mjs','.java','.jar','.class') or '__pycache__' in p.parts)]
 table(study/'collection_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(study.rglob('*')) if p.is_file() and p.name not in ('collection_manifest.tsv','CLOSED')])
 (study/'CLOSED').write_text(sha(study/'collection_manifest.tsv')+'\nCONDITIONAL_GENERATOR_VALIDATED_SINGLE_DATASET_FROZEN\nFITTING_NOT_STARTED\n')
 print('CONDITIONAL_PACKAGE_CLOSED',spent,'core-seconds; one dataset; zero fits')
if __name__=='__main__':main()
