"""Stage A metadata/source compatibility audit only. No model/tool execution."""
from pathlib import Path
import argparse,csv,datetime,hashlib,os,shutil,subprocess,sys

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with Path(p).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tab(p,rows,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def kv(p,values):tab(p,[dict(key=k,value=str(v).lower() if isinstance(v,bool) else str(v)) for k,v in values.items()],['key','value'])

def main(project,study,benchmark):
 assert not study.exists(),'New namespace required; no overwrite of closed evidence'
 for d in ['inputs','reports','accounting','releases/r01']: (study/d).mkdir(parents=True)
 owned=[Path(__file__).resolve()]
 plans=['culture_structural_step04_tool_first_revision_plan_20261114.md','culture_structural_tools_literature_review_20261114.md','meeting_aligned_synthetic_identifiability_master_plan_20261112.md']
 for name in plans:
  source=project/Path({'culture_structural_step04_tool_first_revision_plan_20261114.md': 'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_tool_first_revision_plan_20261114.md', 'culture_structural_tools_literature_review_20261114.md': 'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_tools_literature_review_20261114.md', 'meeting_aligned_synthetic_identifiability_master_plan_20261112.md': 'oxygen/results/analysis_documentation/frozen_inputs/in_vitro/design_and_contracts/meeting_aligned_synthetic_identifiability_master_plan_20261112.md'}[name]);shutil.copy2(source,study/'inputs'/name);owned.append(source)
 role_path=project/'oxygen/results/structural_identifiability/20261113_culture_symbolic_01/inputs/parameter_roles.tsv';roles=read(role_path);assert len(roles)==20 and sum(r['scientific_role']=='mechanism' for r in roles)==14
 shutil.copy2(role_path,study/'inputs/parameter_inventory.tsv');owned.append(role_path)
 old=project/'oxygen/results/structural_identifiability/20261114_culture_reduction_revision_01'
 previous=read(old/'inputs/external_sources.tsv');proof=[]
 for row in previous:
  p=project/row['path'];actual=sha(p);assert actual==row['sha256'],'Frozen source drift: '+row['path']
  proof.append(dict(path=row['path'],expected_sha256=row['sha256'],actual_sha256=actual,pass_gate=True));owned.append(p)
 tab(study/'reports/frozen_source_verification.tsv',proof)
 sources={str(p.relative_to(project)):dict(path=str(p.relative_to(project)),sha256=sha(p),bytes=p.stat().st_size) for p in owned if p.is_relative_to(project)}
 tab(study/'inputs/source_manifest.tsv',list(sources.values()))
 commit=subprocess.check_output(['git','--git-dir='+str(benchmark),'rev-parse','HEAD'],text=True).strip()
 listing=subprocess.check_output(['git','--git-dir='+str(benchmark),'ls-tree','-r','--name-only','HEAD'],text=True).splitlines()
 selected=[p for p in listing if p.startswith(('EAR/','Combos/')) or p=='README.md']
 (study/'reports/benchmark_tree_listing.txt').write_text('\n'.join(selected)+'\n')
 assert not any(Path(p).name.lower().startswith(('findcombos','identifiabilityanalysis')) for p in listing)
 caller=subprocess.check_output(['git','--git-dir='+str(benchmark),'show','HEAD:EAR/HIV2.txt'],text=True)
 assert 'Needs["IdentifiabilityAnalysis`"]' in caller
 kv(study/'reports/benchmark_repository_audit.tsv',dict(origin='https://github.com/Xabo-RB/Benchmarking_files',commit=commit,scope='Git tree and selected text read only; no checkout, package import or execution',EAR_files=sum(p.startswith('EAR/') for p in selected),COMBOS_files=sum(p.startswith('Combos/') for p in selected),core_source_found=False,root_license_found=any(p.lower() in ('license','license.txt','copying','copying.txt') for p in listing),external_package_call_present=True,notebook_files_executed=0,forbidden_files_checked_out=0))
 searches=[
 ('S01','local workspace','oxygen/code, docs, scripts, container; bounded depth 7','No matching core package filenames'),
 ('S02','local supplied downloads','Downloads; bounded depth 2','No matching core package filenames'),
 ('S03','HPC read-only metadata','Project oxygen/code, tools, vendor and Docker; bounded depth 7','No matching core package filenames; no allocation or runtime call'),
 ('S04','public search','FindCombos.maxima; COMBOS Maxima source/download','Paper/component references; no verified core archive'),
 ('S05','public search','IdentifiabilityAnalysis.m; IdentifiabilityAnalysis GitHub/kernel/version','Provider and calling examples; no verified application core'),
 ('S06','provider page','https://www.fcc.chalmers.se/software/other-software/identifiabilityanalysis/','Provider directs package requests to author contact; no request sent'),
 ('S07','primary paper','https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0110261','Published code/license acquisition is author request; no request sent'),
 ('S08','published service/component URL','https://biocyb1.cs.ucla.edu/combos/ and /maxima/FindCombos.maxima','Not accessible via browsing tool; does not establish global unavailability'),
 ('S09','benchmark repository',commit,'Input/calling files only for the reviewed tools; no matching application cores'),
 ('S10','name disambiguation','arcos-combos ComBoS search result','Unrelated distributed-computing simulator; excluded'),
 ('S11','local application directories','Existing native Wolfram/Mathematica application roots; bounded depth 4','Library/Wolfram existed; no matching package file or directory'),
 ('S12','HPC application directories','Existing user Wolfram/Mathematica application roots, Downloads and project; bounded depth 4','Only project candidate root existed; no matching application package; no runtime launched')]
 tab(study/'reports/source_search_log.tsv',[dict(search_id=i,surface=s,scope=q,observed_result=o,checked_utc=utc()) for i,s,q,o in searches])
 contacts={'COMBOS':'joed@ucla.edu (published route; current responsiveness untested)','IdentifiabilityAnalysis':'identifiabilityanalysis@fcc.chalmers.se (current provider instruction)'}
 tab(study/'reports/source_access.tsv',[dict(tool=t,core_obtained=False,license_text_obtained=False,source_audit_pass=False,acquisition_route=contacts[t],author_contact_authorized=False,author_contact_attempted=False,status='access_hold',reason='No verified core plus license obtained through permitted public/local routes') for t in contacts])
 tab(study/'reports/dependency_audit.tsv',[
  dict(tool='COMBOS_web',basis='Published implementation, not executed source',policy_status='excluded',detail='Published frontend uses JavaScript and JSON examples; not downloaded or executed'),
  dict(tool='COMBOS_local_core',basis='Candidate architecture only',policy_status='unverified_source_blocked',detail='Native SBCL/ECL Maxima candidate; ABCL/JVM excluded; no acquired core imports to inspect'),
  dict(tool='IdentifiabilityAnalysis_kernel',basis='Published package and kernel documentation only',policy_status='unverified_source_blocked',detail='Direct native kernel candidate; package imports/JLink/serialization and license not inspected')])
 tab(study/'reports/runtime_gate.tsv',[dict(tool=t,source_license_gate=False,static_package_gate='not_run_access_hold',runtime_gate='not_run_access_hold',xxlarge_eligible=False,process_compliance='not_certified',package_execution_count=0) for t in contacts])
 tab(study/'reports/process_audit.tsv',[dict(scope='upstream analysis applications',executed=False,process_tree='not_applicable',reason='No eligible core/license; no application process was started',HPC_allocations=0)])
 controls=[]
 for tool in contacts:
  for case in ['identifiable_decay','nonidentifiable_parameter_sum','rational_parameterized_initial_map','unsupported_original_map_control']:
   controls.append(dict(tool=tool,control=case,status='not_run_access_hold',job_id='none',hours_reserved=0,observed_result='none',reason='Source/license gate did not pass; no replacement toy or unsupported project task submitted'))
 tab(study/'reports/control_results.tsv',controls)
 base='oxygen/results/structural_identifiability/20261109_seed144_wp0/source/code/oxygen/code/O2_supply_demand_MAP/'
 R=base+'model/model_O2_supply_demand_MAP.R';C=base+'model/model_O2_supply_demand_MAP.cpp';I=base+'util/o2_supply_demand_map_invitro_lineage_simulation_utils.R';O=base+'util/o2_supply_demand_map_invitro_objective_utils.R'
 specifications={
 'lam_max':('rational_prefactor','lam_max/(1+alpha_o2*h*(N/44)^gamma_growth)',R,724,'Joint growth map remains nonrational; no parameter isolated'),
 'p_misseg':('rational_in_effective_death','p_mis_base+p_misseg*mu/(k_o_mis+mu)',R,817,'Exact zero removes induction; k_o_mis inactive in this null'),
 'k_o_mis':('rational_in_effective_death','mu/(k_o_mis+mu)',R,817,'Inactive at p_misseg=0; not an identified fixed value'),
 'buffer_smax':('prefactor_of_nonrational_function','buffer_smax*exp(-buffer_beta*(44/N)^buffer_n_exp)',C,533,'All three buffering coordinates remain coupled'),
 'buffer_beta':('nonrational_exponential_parameter','exp(-buffer_beta*(44/N)^buffer_n_exp)',C,533,'Exact zero removes parent-ploidy dependence, not all error-size dependence'),
 'buffer_n_exp':('nonrational_unknown_power','(44/N)^buffer_n_exp inside exponential',C,533,'Inactive at buffer_beta=0; do not demand recovery there'),
 'p_wgd':('rational_division_probability','one doubled lineage with probability p_wgd',C,945,'Retain full generator and mother removal; no structural label'),
 'alpha_o2':('rational_prefactor_with_nonrational_coupling','alpha_o2*h*(N/44)^gamma_growth',R,724,'Fixed oxygen does not establish rational original parameter dependence'),
 'gamma_growth':('nonrational_unknown_power','(N/44)^gamma_growth',R,724,'Known integer N does not make dependence on unknown gamma rational'),
 'mu_hp':('rational_prefactor_with_nonrational_coupling','mu_hp*h*(1+above_dip^gamma_mu)',R,765,'Preserve stress-factor coupling across oxygens'),
 'gamma_mu':('nonrational_unknown_power','above_dip^gamma_mu',R,765,'Production diploid_NULL switch is fixed; exponent remains unknown'),
 'O2_crit':('nonrational_joint_Hill_map','O2_crit^n_O/(O2_crit^n_O+O^n_O)',R,724,'Positive O2_crit; joint unknown exponent cannot be silently fixed'),
 'n_O':('nonrational_unknown_power','O2_crit^n_O/(O2_crit^n_O+O^n_O)',R,724,'At O=0 and n_O>0, h=1; this special input is not the multi-input proof'),
 'p_mis_base':('rational_probability_entry','p_mis_base plus induced increment',R,817,'Retain baseline at p_misseg=0; generic estimate flags do not override map'),
 'sigma_growth':('outside_deterministic_mean_map','Gaussian growth-score variance',base+'util/o2_supply_demand_map_invitro_objective_utils.R',387,'Requires probabilistic observation analysis; no mean-output identification claim'),
 'sigma_kary':('nonrational_unknown_observation_width','normalized Gaussian observation-kernel columns',O,408,'Unknown in S7; absent from true-composition outputs, not identified by omission'),
 'init_mean_2N':('nonrational_initial_manifold','finite-grid normalized Gaussian weights',I,52,'Known in S1-S5, unknown jointly in S6-S7; shared by cohort across ideal oxygen cases'),
 'init_sd_2N':('nonrational_initial_manifold','finite-grid normalized Gaussian weights',I,52,'Positive SD; not interchangeable with independent state weights'),
 'init_mean_4N':('nonrational_initial_manifold','finite-grid normalized Gaussian weights',I,52,'Separate cohort family, same known/unknown case distinction'),
 'init_sd_4N':('nonrational_initial_manifold','finite-grid normalized Gaussian weights',I,52,'Preserve normalization and finite support; positive SD')}
 assert set(specifications)=={r['parameter'] for r in roles}
 role_rows=[]
 for r in roles:
  kind,expression,path,line,note=specifications[r['parameter']]
  role_rows.append(dict(parameter=r['parameter'],scientific_role=r['scientific_role'],units=r['units'],lower=r['lower'],upper=r['upper'],transform=r['transform'],expression_class=kind,source_expression=expression,source_file=path,source_line=line,source_sha256=sha(project/path),handling=note,structural_result='not_analyzed',Stage_A_estimation='none'))
 tab(study/'reports/parameter_and_nuisance_audit.tsv',role_rows)
 features=[
 ('F01','state_grid','133 live states, N=22:154','Large size; no complexity benchmark performed','none'),
 ('F02','joint_dynamics','Unknown powers in growth, death and Hill response; exponential buffering','Not a rational function of the original joint unknowns','R1 model assumptions; R2 Sec2 p942'),
 ('F03','finite_binomial_kernel','Integer chromosome/error counts; integer powers in p and survival','Rational coefficient layer possible; original coupled coefficient manifold still required','Frozen C++ generator'),
 ('F04','initial_family','Four unknown Gaussian mean/SD coordinates with normalization','Nonrational initial map; rational initial-map support does not cover it directly','R2 Sec2 and Sec3.1; frozen initial helper'),
 ('F05','observations','Unknown sigma_kary in Gaussian kernel; derived fixed flow width','Nonrational unknown observer map; no extra flow nuisance','Frozen observation helper'),
 ('F06','variance','sigma_growth affects variance/score, not deterministic mean outputs','Outside this mean-output question; no identification by fixing/omission','Frozen objective'),
 ('F07','branches','Crowding/clipping and parameter-dependent pruning masks','Need valid analytic region and threshold inequalities; no global rationality claim','Frozen source'),
 ('F08','propagation','D: production DT=0.05; C: separately declared ODE','Reviewed ODE methods do not establish discrete-model identifiability','R1/R2 ODE scope; closed DT evidence'),
 ('F09','passage_protocol','Count-selected end day and inherited composition','Do not replace hybrid/sparse design by an ideal ODE observation claim','Audited schedule and lineage helper'),
 ('F10','boundary_domain','Finite 22:154 drop support; separate biological/numerical sinks','Domain-conditional conclusions; old positive-tail holds retained','Closed reduction audit'),
 ('F11','exact_reparameterization','eta=Phi(theta), initial and observer constraints retained','No accepted exact supported transformation established in Stage A','Compatibility audit only'),
 ('F12','prior_reduction','5/10/20-state adequacy gates failed','Do not use an unvalidated coarse rational surrogate to label full mechanisms','Closed 20261114 results')]
 tab(study/'reports/rationality_ledger.tsv',[dict(feature_id=i,feature=f,source_behavior=b,compatibility_consequence=c,evidence=e,status='static_audit_complete_not_structural_analysis') for i,f,b,c,e in features])
 cases=[('S1','0.2','L','known initial composition','14 mechanisms'),('S2','0.2','L and true composition','known initial composition','14 mechanisms'),('S3','0.2;1;20.5','L','known initial composition','14 mechanisms'),('S4','0.2;1;20.5','L and true composition','known initial composition','14 mechanisms'),('S5','0;0.2;1;20.5','L and true composition','known initial composition','14 mechanisms'),('S6','0;0.2;1;20.5','L and true composition','unknown four-coordinate Gaussian family','14 mechanisms plus four initial nuisances'),('S7','0;0.2;1;20.5','L plus expected karyotype and flow','unknown four-coordinate Gaussian family','14 mechanisms plus four initial nuisances plus sigma_kary')]
 case_rows=[]
 for tool in contacts:
  for cid,ox,obs,initial,unknown in cases:
   for form in ['C','D']:
    case_rows.append(dict(tool=tool,case_id=cid,model_form=form,oxygen_pct=ox,outputs=obs,initial_role=initial,unknowns=unknown,source_license_gate='hold',direct_mapping='not_supported_by_published_original_parameter_contract',additional_obstacle='Nonrational dynamics'+('; nonrational initial manifold' if cid in ('S6','S7') else '')+('; nonrational observer' if cid=='S7' else '')+('; no discrete-time guarantee in reviewed tool' if form=='D' else ''),accepted_exact_lift=False,compute_released=False,structural_result='not_analyzed'))
 tab(study/'reports/per_case_backend_status.tsv',case_rows)
 tab(study/'reports/initial_observation_constraints.tsv',[dict(component='initial_four',constraint='Normalized finite-grid Gaussian; four coordinates shared by cohort as declared',known_case='S1-S5',unknown_case='S6-S7',replacement_rejected='133 unrelated free initial weights'),dict(component='sigma_kary',constraint='Column-normalized Gaussian chromosome-observation kernel',known_case='Not a true-composition parameter in S1-S6',unknown_case='S7',replacement_rejected='True composition supplied in place of measured operator'),dict(component='sigma_growth',constraint='Variance-only role in existing growth score',known_case='No mean-output identification claim',unknown_case='Requires separate probability model',replacement_rejected='A fitted scale declared identified because fixed'),dict(component='flow_width',constraint='Configured/derived production setting with grid-step normalization',known_case='Fixed setting',unknown_case='No seventh fitted nuisance',replacement_rejected='An invented fitted flow coordinate'),dict(component='initial_abundance',constraint='Known input in proposed ideal cases; not a measured-scale claim',known_case='All proposed cases conditional on known count',unknown_case='Not analyzed here',replacement_rejected='Transfer to uncertain absolute inoculum scale')])
 gates=[dict(gate='A0_source_and_legal_access',audit_complete=True,execution_gate_pass=False,outcome='access_hold',detail='Neither application core/license obtained; no author contact'),dict(gate='A1_runtime_compliance',audit_complete=True,execution_gate_pass=False,outcome='not_run_access_hold',detail='Gate disposition recorded; no process/import compliance certificate claimed'),dict(gate='A2_exact_model_class_map',audit_complete=True,execution_gate_pass=False,outcome='mapping_hold',detail='All14+6 covered; unchanged joint map is outside the reviewed rational contracts; no accepted lift'),dict(gate='A3_reuse_decision',audit_complete=True,execution_gate_pass=False,outcome='complete_with_HOLD',detail='No backend/project run released; source and exact mapping prerequisites both unresolved')]
 tab(study/'reports/stage_A_acceptance.tsv',gates)
 tab(study/'accounting/budget_ledger.tsv',[dict(phase=p,authorized_core_hours=h,actual_allocated_core_seconds=0,live_reserved_core_seconds=0,remaining_core_hours=h) for p,h in [('runtime_compliance',4),('control_jobs',8),('mapping_verification_report',4)]])
 tab(study/'accounting/jobs.tsv',[],['task','job_id','qos','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
 tab(study/'accounting/call_ledger.tsv',[],['kind','call_id','status','reason'])
 kv(study/'reports/final_status.tsv',dict(stage_A_execution_complete=True,scientific_status='HOLD_access_and_exact_mapping',eligible_backends=0,applications_executed=0,source_license_gates_passed=0,known_example_tests_run=0,HPC_jobs_submitted=0,actual_allocated_core_hours=0,live_reserved_core_hours=0,authorized_core_hour_cap=16,author_contacts=0,model_replays=0,model_simulations=0,synthetic_datasets=0,fits=0,profiles=0,tumor_fits=0,custom_solver_implemented=False,structural_classification_complete=False,Stage_B_released=False,Step_5_released=False,closed_utc=utc()))
 kv(study/'inputs/authorization.tsv',dict(user_confirmation_date='2026-11-16 America/New_York',scope='Step 4 Stage A only',budget='16 independent allocated core-hours',qos='xxlarge if eligible',application_concurrency_cap='none',author_contact=False,java=False,javascript=False,json=False,Stage_B=False,Step_5=False,old_controllers_and_tasks='Do not modify'))
 release=study/'releases/r01';shutil.copy2(Path(__file__),release/'audit_stage_a.py')
 kv(study/'reports/audit_runtime.tsv',dict(python=sys.version.split()[0],role='Local standard-library metadata/hash/table writer only; not a scientific solver',upstream_code_executed=False,HPC_access='Read-only bounded filename inventory; no scheduler or runtime invocation'))
 print('AUDIT_TABLES_READY',len(proof),'frozen source files',len(role_rows),'roles',len(case_rows),'backend/case/form rows',commit)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('project',type=Path);p.add_argument('study',type=Path);p.add_argument('benchmark',type=Path);a=p.parse_args();main(a.project.resolve(),a.study.resolve(),a.benchmark.resolve())
