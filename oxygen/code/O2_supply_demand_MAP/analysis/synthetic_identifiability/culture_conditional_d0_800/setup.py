"""Register a new immutable release locally; does not submit or simulate."""
from pathlib import Path
import shutil,sys
import workflow as w

def setup(project):
    project=Path(project).resolve();run=project/'oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02'
    run.mkdir(exist_ok=False)
    for name in ('accounting','checks','data','design','dispatch','engineering','fits','inputs','logs','manifests','reports','runtime','state','tables','truth','collection'): (run/name).mkdir()
    remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability'
    w.typed(run/'study.tsv',dict(project=remote,study_id=run.name,scope='D0_production_conditioned_expected',truth='seed144',free_mechanisms=14,known_nuisance=6,physical_sampling_release=False,old_run_closed='20261113_culture_d0_recovery_01',old_run_actual_core_seconds=23))
    w.typed(run/'authorization.tsv',dict(total_core_hours=800,qos='xxlarge',partition='red',application_concurrency_cap='none',requeue=False,blind_primary_origins=4,reserve_rounds=1,reserve_packages=2,infrastructure_retries=2,noisy_datasets=0,profiles=0,new_tumor_fits=0,truth_selection='User-prespecified seed144 before generated outcomes',language='English',formats='TSV CSV RDS Markdown PNG PDF',runtime='R Python shell native C++; no Java JavaScript or JSON'))
    controls=dict(NP=256,global_min_generations=300,global_max_generations=1000,global_worker_max_callbacks=(256*1001+6)//7+512,global_soft_seconds=41400,local_soft_seconds=20700,local_max_callbacks=20000,independent_score_tolerance=1e-8,attainability_delta=1e-6,max_growth_standardized_error=.01,rms_growth_standardized_error=.001,max_kary_tv=.001,max_flow_tv=.001,directional_improvement_tolerance=1e-7,parameter_transformed_resolution=.01)
    controls.update({'generation_seed_'+o:2026111311+11*i for i,o in enumerate('ABCDEF')});w.typed(run/'tables/search_controls.tsv',controls)
    source=Path(__file__).parent;release=run/'releases/r01';shutil.copytree(source,release,ignore=shutil.ignore_patterns('__pycache__'))
    files=sorted(p for p in release.rglob('*') if p.is_file())
    w.table(release/'manifest.tsv',[dict(path=str(p.relative_to(release)),sha256=w.sha(p)) for p in files],['path','sha256'])
    wp0=project/'oxygen/results/structural_identifiability/20261109_seed144_wp0'
    external=[]
    for directory in ('source/code/oxygen/code/O2_supply_demand_MAP/model','source/code/oxygen/code/O2_supply_demand_MAP/util'):
        external.extend(p for p in (wp0/directory).rglob('*') if p.is_file() and not any(x.startswith('.rcpp_cache') for x in p.parts) and p.suffix in ('.R','.cpp','.h','.hpp'))
    external += [wp0/'source/seed144'/n for n in ('fit_result.rds','parameter_table_input.csv')]
    external += [wp0/'source/inputs/oxygen'/p for p in ('ploidyOxygen/data/fit_objects/fit_data.Rds','ploidyOxygen/data/fit_objects/jobs_2N.Rds','ploidyOxygen/data/fit_objects/jobs_4N.Rds','data/g0g1_ploidy_density_grid.csv')]
    external += [project/'oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson/implementation'/n for n in ('helpers.R','diagnostics.R')]
    external += [project/'scripts/agentRrunner.sh']
    rows=[dict(path=str(p.relative_to(project)),sha256=w.sha(p)) for p in sorted(set(external))]
    w.table(run/'manifests/external_manifest.tsv',rows,['path','sha256'])
    plan='''# Approved conditional D0 engineering and blind recovery

The independent ceiling is 800 allocated core-hours: engineering24, positive6,
primary420, reserve210, report12, infrastructure128. No borrowing across phase
ceilings, no requeue, no application concurrency cap, red/xxlarge only.
The closed earlier128-hour package and its23 consumed core-seconds are historical
evidence; neither its failed-start allowance nor its old stopping policy is reset.
This new package follows the separately confirmed800-hour plan.

Freeze the full-precision validated seed144 source, effective roundoff mapping,
all14 mechanisms, all6 known generating nuisance, fixed runtime configuration,
lineages, initial absolute scale, assigned oxygen and conditional count histories.
No real-data global-optimum recovery is prerequisite. Boundary truth is retained.

Expected growth uses the production selected-endpoint log ratio/day, which differs
from the experiment's OLS estimator. Karyotype uses all133 probability weights;
flow uses normalized mass on200 grid points. No finite draws, pseudo-cells, FCS,
noise, or death scoring. Actual sample sizes remain design metadata. Candidate
parameters recompute count-conditioned selected days and propagation. Unknown
physical growth/stopping rules, flow processing and live/dead mapping remain held.

Blinding: the fitting bundle contains known nuisance and neutral active defaults,
not the generating active coordinates. Four independent origin RNG streams use
a registered50:50 mixture of full transformed-domain uniform proposals and the
previous engineering growth-feasibility proposal. The latter replaces lam_max
with Uniform(1,1.5), alpha_o2(.5,.75), gamma_growth(1.5,2), mu_hp(.005,.01),
gamma_mu(1.2,1.6), buffer_smax(.8,1); every fourth such draw uses p_misseg
Uniform(0,.003) and log10(p_mis_base) Uniform(-6,-4). Other coordinates retain
full-domain draws. These proposals affect starts only, never bounds or priors.
Screen first feasible in draw order, not by objective or truth distance, at most
2000 candidate calls or3300 soft seconds. Require8 feasible anchors. Construct
the256-member population from registered draws with those first8 anchors, and
verify every actual DEoptim initial input against it. No truth/old-fit injection.

DEoptim uses NP256, minimum early-stop patience300, maximum1000 generations,
8 allocated CPUs with7 objective workers, maximum12 hours; soft cutoff41400s.
Local SBPLX then bounded pattern search uses at most20000 callbacks/6 hours,
soft cutoff20700s. Independent checks replay source vectors and test both signs
of every coordinate at meshes.001/.0001 against the same selected vector.
Positive control is truth-started, separately labelled and never independent.

Same-vector score difference <=1e-8; raw signed Q-Qtruth may not be below-1e-8.
Attainability delta <=1e-6, growth standardized RMS <=.001 and max <=.01,
per-sample karyotype and flow-mass TV <=.001. Final directional improvement
must not exceed1e-7. At least2 independent blind origins must meet all gates.
Parameter errors and transformed-domain1% engineering resolution are reported
separately; nonunique parameter recovery alone is not a bridge failure or a
finite structural proof. Lowest verified score selects estimates, never truth
distance. Preserve original native status, budget truncation, all failed attempts,
exact starts, boundary activity and separate best-found versus checked vectors.

Run first complete blind chain before releasing the other3. One reserve round
contains at most2 packages; no indefinite retries. Infrastructure failures get
at most2 explicit reviewed resubmissions within128 hours; ambiguous submission
requires receipt reconciliation. Source/interface mismatch stops release.
After4 primary and at most2 reserve origins, stop with PASS or optimizer-unresolved.
Do not generate noisy datasets, profiles or tumor fits. All outputs are English.
No Java, JavaScript, JSON/JSONL/NDJSON, forbidden parsers or hidden table payloads.
'''
    w.text(run/'APPROVED_PLAN.md',plan);w.text(run/'STATUS.md','# Conditional D0 execution\n\nRegistered; engineering not yet released.\n')
    print(run)
if __name__=='__main__':setup(sys.argv[1])
