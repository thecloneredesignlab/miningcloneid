"""Close the registered reduction hold; metadata only, no model/tool execution."""
from pathlib import Path
import argparse,csv,re,shutil,subprocess,tarfile
import workflow as w
def kv(p):return {r['key']:r['value'] for r in w.read(p)}
def write_kv(p,data):w.table(p,[dict(key=k,value=str(v).lower() if isinstance(v,bool) else str(v)) for k,v in data.items()],['key','value'])
def close(s,release):
 assert not (s/'COMPLETE').exists(),'Already closed'
 settings=kv(s/'study.tsv');project=Path(settings['project'])
 for r in w.read(s/'inputs/manifest.tsv'):assert w.sha(s/'inputs'/r['path'])==r['sha256']
 for r in w.read(s/'inputs/external_sources.tsv'):assert w.sha(project/r['path'])==r['sha256']
 rows=w.reconcile(s);assert all(r['state'] in w.TERMINAL for r in rows),'Owned allocation still active'
 native=s/'results/reduction_native_rate_exact';result=kv(native/'result.tsv')
 assert result['source_and_kernel_integrity']=='true'
 assert result['production_to_reduction_approximation_pass']=='false' and result['coarse_projection_approximation_pass']=='false'
 assert result['fine_tail_routing_pass']=='false'
 assert not any(r['phase']=='structural' for r in rows),'Unexpected structural dispatch'
 assert (s/'results/tool_postfreeze/CONSTRUCTOR_PASS').exists()
 for r in w.read(native/'payload_manifest.tsv'):assert w.sha(native/r['path'])==r['sha256']
 assert (native/'COMPLETE').read_text().strip()==w.sha(native/'payload_manifest.tsv')
 calls=w.read(s/'state/native_calls.tsv');assert len(calls)==24 and len(calls)<=48
 case_rows=w.read(native/'reduction_validation.tsv');account=w.read(native/'accounting_and_native_checks.tsv');assert len(account)==24
 summary=[]
 for comparison in sorted({x['comparison'] for x in case_rows}):
  group=[x for x in case_rows if x['comparison']==comparison]
  summary.append(dict(comparison=comparison,cases=len(group),engineering_passes=sum(x['engineering_resolution_pass']=='TRUE' for x in group),max_class_tv=max(float(x['max_class_tv']) for x in group),max_abs_log_live_error=max(float(x['max_log_live_error']) for x in group)))
 w.table(s/'reports/reduction_summary.tsv',summary,list(summary[0]))
 w.table(s/'reports/reduction_validation.tsv',case_rows,list(case_rows[0]))
 cases=w.read(release/'case_manifest_r04.tsv');held=[]
 for c in cases:
  held.append(dict(case=c['case'],model=c['model'],outputs=c['kind'],oxygen_levels=c['oxygens'],initial_unknown=c['unknown'],initial_type=c['initial'],variant=c['variant'],status='Not released: reduction adequacy gate failed',parameter_classification='Unresolved; no solver or exact-jet certificate executed',scope='Candidate reduced model only'))
 w.table(s/'reports/structural_identifiability_results.tsv',held,list(held[0]))
 additional=w.read(s/'inputs/additional_cases.tsv');w.table(s/'reports/additional_case_status.tsv',[dict(**r,status='Conditional analytical discussion only; hybrid/project gate not passed') for r in additional],list(additional[0])+['status'])
 mappings=[('lam_max','bmax','Division meaning retained; production response shapes still differ'),('p_mis_base','p0','Per-chromosome baseline role retained'),('p_misseg','a','Death-linked increment retained with simplified death response'),('k_o_mis','k','Death-risk half-saturation role retained'),('p_wgd','w','One retained doubled lineage per WGD division'),('buffer_smax','s2;s4','Three original buffering coordinates represented by two anchors'),('buffer_beta','s2;s4','No unique original-coordinate inverse established'),('buffer_n_exp','s2;s4','No unique original-coordinate inverse established'),('alpha_o2','alpha','Growth exponent fixed to one'),('gamma_growth','Fixed at one','Original exponent not analyzed'),('mu_hp','dH','Death exponent fixed to one'),('gamma_mu','Fixed at one','Original exponent not analyzed'),('O2_crit','K_O','Resource Hill exponent fixed to one; not carrying capacity'),('n_O','Fixed at one','Original exponent not analyzed')]
 w.table(s/'reports/parameter_reduction_map.tsv',[dict(full_parameter=n,candidate_reduced_quantity=q,reduction_assumption=a,full_model_identifiability='Not established',step05_release=False) for n,q,a in mappings],['full_parameter','candidate_reduced_quantity','reduction_assumption','full_model_identifiability','step05_release'])
 nuisance=[]
 for n in ['init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N']:
  nuisance.append(dict(parameter=n,structural_role='Known or independent unknown coarse initial vectors; not the production Gaussian parameterization',status='Original nuisance identifiability not established'))
 nuisance.extend([dict(parameter='sigma_growth',structural_role='Absent from ideal mean outputs; requires a probabilistic variance model',status='Not analyzed'),dict(parameter='sigma_kary',structural_role='Ideal true composition omits the observation convolution',status='Not analyzed')])
 w.table(s/'reports/nuisance_role_status.tsv',nuisance,list(nuisance[0]))
 transforms=[dict(id='T01',model='R5D2-C and R5D2-L',assumptions='One known positive constant oxygen; positive interior',transformation='K_O->K_O_new; alpha->alpha*h/h_new; dH->dH*h/h_new; all other coordinates unchanged',preserved='Every b,d,p,s,w and therefore full state/output dynamics',status='Exact conditional analytical equivalence; not transferred to full simulator'),dict(id='T02',model='R5D2-L',assumptions='Only normalized composition; unknown initial absolute abundance',transformation='X->c*X for each independent experiment',preserved='Normalized live fractions',status='Exact conditional analytical equivalence; not a symmetry of primary known-capacity model'),dict(id='T03',model='R5D2-C with hypothetical unknown capacity',assumptions='Capacity is an additional unknown; fraction-only outputs',transformation='X->c*X; Ccap->c*Ccap',preserved='Normalized fractions and crowding factor',status='Conditional setting-sensitivity argument; capacity was not freed or fitted')]
 w.table(s/'reports/structural_transformations.tsv',transforms,list(transforms[0]))
 (s/'reports/identifiable_combinations.md').write_text('''# Candidate combinations and analytical limits\n\nFor a single prescribed positive oxygen, alpha*h(O;K_O) and dH*h(O;K_O) are invariant under the explicit T01 transformation. This proves ambiguity of the separate response-scale coordinates in the candidate rational model; it does not establish that the two invariants are recoverable from the existing outputs.\n\nIn the linear, generic minimal-order-five, growth-only model, characteristic-polynomial coefficients are ideal-output combinations. No project Hankel/local certificate was executed after the reduction gate failed, so the registered finite-jet route remains unexecuted. This statement is not transferred to the primary known-capacity nonlinear model.\n\nThe hypothetical event-augmented inverse in proofs/global_event_augmented_inverse.md is a conditional sufficient-output argument for the ten-coordinate reduction, pending scientific review. It is not a result for actual culture observations or all fourteen original mechanisms. No primary recovery estimand or Step05 design is released.\n''')
 for f in ['model_addendum_r04.md','case_manifest_r04.tsv']:
  shutil.copy2(release/f,s/'symbolic_models'/f)
 for f in ['exact_model.py','julia_case.jl']:shutil.copy2(release/f,s/'symbolic_models'/f)
 shutil.copy2(s/'inputs/model_definition.md',s/'symbolic_models/model_definition_R5D2_L_historical.md')
 ids=','.join(r['job_id'] for r in rows)
 raw=subprocess.run([str(w.SLURM/'sacct'),'-n','-P','-j',ids,'--format=JobIDRaw,State%30,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,TotalCPU,QOS'],check=True,stdout=subprocess.PIPE,universal_newlines=True).stdout
 (s/'accounting/final_sacct_raw.txt').write_text(raw)
 total=sum(int(r['actual']) for r in rows)
 write_kv(s/'accounting/final_accounting.tsv',dict(allocated_core_seconds=total,allocated_core_hours=total/3600,authorized_core_hours=200,operational_phase_cap_core_hours=168,active_owned_parents=0,owned_parents=len(rows),infrastructure_retries=sum(r['phase']=='infrastructure' for r in rows),native_segments=len(calls),reused_native_segments=23,new_segments_after_reference_correction=1,symbolic_allocations=0,synthetic_datasets=0,recovery_fits=0,profiles=0,tumor_fits=0,step05_released=False))
 for p in (s/'runtime').rglob('*'):
  if p.is_file() and p.name in ['model_source_manifest.tsv','image_sha256.txt','launch_exit_code.txt']:
   target=s/'reports/runtime_provenance'/p.relative_to(s/'runtime');target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 for p in (s/'results/tool_setup_retry1').glob('*'):
  if p.is_file():target=s/'reports/tool_provenance'/p.name;target.parent.mkdir(exist_ok=True);shutil.copy2(p,target)
 for name in ['Project.toml','Manifest.toml']:
  target=s/'reports/tool_provenance'/name;shutil.copy2(s/'toolchain/environment'/name,target)
 forbidden={'.json','.jsonl','.ndjson','.js','.java','.jar','.class'}
 removed_archive=[]
 for p in (s/'toolchain').glob('julia-*.tar.gz'):
  with tarfile.open(p) as archive:bad=[m.name for m in archive.getmembers() if m.isfile() and Path(m.name).suffix.lower() in forbidden]
  if bad:removed_archive.append(dict(path=str(p.relative_to(s)),sha256=w.sha(p),prohibited_unused_members=len(bad),reason='Official checksum retained; discard downloaded bundle containing unused prohibited fixtures'));p.unlink()
 w.table(s/'reports/discarded_installation_downloads.tsv',removed_archive,['path','sha256','prohibited_unused_members','reason'])
 (s/'reports/structural_limitations.md').write_text('''# Structural limitations and stop decision\n\nThe production runtime/native propagation is valid; the proposed reduction is not an accepted quantitative bridge. The registered response rationalization, state aggregation and fine upper-tail routing all show material limitations. Prespecified fixed weighted lifting and shifted boundaries did not resolve the full gate. No expensive project symbolic elimination or finite-jet classification was released. Tool toy examples and analytical candidate-model equivalences do not satisfy the project structural gate.\n\nAll fourteen original mechanisms remain inventoried. Three original response exponents were fixed by rationalization; three buffering coordinates became two anchors. Four initial-distribution and two observation-error coordinates have not been structurally identified. Known capacity can break a linear abundance symmetry, so nonlinear and linear models cannot exchange classifications. Euler discretization is distinct from continuous structural dynamics.\n\nControlled switching, target-count passage selection and the nonlinear hybrid reset case are unresolved. Actual growth OLS/count dependence, physical stopping rules, precise flow processing/event counts and live/dead classification/time mapping remain open. Ideal class fractions and retained material are not existing assay operators.\n\nStep05 is not released. A revised Step04 reduction/tool package requires review before more model computations: isolate coarsening using production response values, reconsider within-class closure/tail representation, then validate a stated response-function representation. Do not obtain acceptance by dropping failed references/oxygen levels, removing mechanism parameters, loosening thresholds or spending the unused allowance without a new decision.\n''')
 max_native=max(float(r['native_euler_scaled_gap']) for r in account);max_route=max(float(r['upper_tail_route_per_final_live']) for r in account)
 lines=['# Revised culture Step04: stopped at the reduction gate','', 'Closed '+w.utc()+'.','', '**Tool/environment and native-propagation gates passed. Reduction adequacy failed. The main project structural case batch and Step05 were not released.**','', 'This is a registered finite stop, not a structural nonidentifiability verdict for the full simulator. Existing controllers/tasks and old budgets were not changed.','', '## Verified numerical evidence','', '- Julia1.10.13 and StructuralIdentifiability.jl0.5.36 passed identifiable, nonidentifiable and identifiable-combination examples and the fresh explicit-ODE-constructor check.','- The production SIF is unchanged, SHA256 dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372.','- All24 fixed production segments pass the unchanged1e-12 propagation threshold; maximum scaled gap {:.16g}.'.format(max_native),'- Twenty-three saved native segments were hash-verified and reused; only the missing segment was added. The cumulative native-call ledger has24 slots out of48.','- The primary rational model preserves known logistic crowding and constant cell-volume bookkeeping. The linear controls are separate.','', '| Comparison | Cases passing | Maximum class TV | Maximum absolute log-live error |','|---|---:|---:|---:|']
 for r in summary:lines.append('| {} | {}/{} | {:.9g} | {:.9g} |'.format(r['comparison'],r['engineering_passes'],r['cases'],r['max_class_tv'],r['max_abs_log_live_error']))
 lines+=['','Targets: TV<=0.05 and absolute log-live error<=0.1. Fine upper-tail routing/final-live maximum {:.9g}, above the1e-6 target. The routed mass remains explicitly counted; its inadequacy limits the comparator.'.format(max_route),'','## Structural scope','', 'No project solver/local-certificate classification was executed. The twenty-two registered exact-jet cases remain unresolved/not released. Conditional analytical transformations and an ideal event-output inverse are supplied separately; neither transfers to the original fourteen mechanisms. Source/exponent/anchor assumptions and nuisance limitations remain explicit.','', '## Allocation closure','', '| Task | Job | State | Exit | Allocated core-hours |','|---|---:|---|---|---:|']
 for r in rows:lines.append('| {} | {} | {} | {} | {:.9f} |'.format(r['task'],r['job_id'],r['state'],r['exit_code'],int(r['actual'])/3600))
 lines+=['','Actual allocated use: {:.9f}/200 core-hours; all{} owned parents terminal. The one infrastructure retry is consumed. Failed parents remain counted. No recovery fits, profiles, noisy datasets or tumor fits were started.'.format(total/3600,len(rows)),'', 'The source-guard and propagation-reference failures remain preserved. The former exposed production crowding before any forward loop. The latter exposed approximately1e-15 rate-accumulation differences in a weak dead stock; native consumed rates resolved propagation without changing parameters or relaxing the state threshold.','', '## Delivery and next decision','', 'Scientific payloads, frozen releases, assumptions, failure logs, manifests and typed tables are deliverables. Large software/native caches remain on HPC with software versions/archive hashes and runtime provenance. A separate local collection receipt determines collection acceptance. The root COMPLETE marker means execution/manifest closure, not scientific passage of Step04.','', 'Review a revised Step04 reduction before Step05. This run releases no further computation.']
 (s/'REPORT.md').write_text('\n'.join(lines)+'\n')
 write_kv(s/'reports/step04_acceptance.tsv',dict(tool_gate=True,source_native_gate=True,reduction_adequacy_gate=False,main_structural_case_gate=False,conditional_analytical_evidence=True,full_fourteen_structural_status='Unresolved',step05_release=False,execution_closed=True))
 write_kv(s/'reports/STOP.tsv',dict(reason='reduction_adequacy_gate_failed',utc=w.utc(),dependent_symbolic_dispatches=0,step05_release=False,requires='Review a revised Step04 reduction'))
 manifest=[]
 for p in sorted(s.rglob('*')):
  rel=p.relative_to(s)
  if not p.is_file() or any(part in ('runtime','toolchain','state','collection','__pycache__') for part in rel.parts) or (len(rel.parts)==1 and p.name in ('delivery_manifest.tsv','delivery_verification.tsv','COMPLETE')):continue
  assert p.suffix.lower() not in forbidden,'Prohibited scientific artifact: '+str(rel)
  if p.suffix in ('.py','.R','.jl','.sh','.md','.tsv','.csv','.txt','.toml'):
   content=p.read_text();assert not re.search('[\u3400-\u9fff]',content),'Non-English scientific text: '+str(rel)
  manifest.append(dict(path=str(rel),sha256=w.sha(p),bytes=p.stat().st_size))
 w.table(s/'delivery_manifest.tsv',manifest,['path','sha256','bytes'])
 for row in manifest:assert w.sha(s/row['path'])==row['sha256']
 write_kv(s/'delivery_verification.tsv',dict(pass_gate=True,files=len(manifest),bytes=sum(r['bytes'] for r in manifest),manifest_sha256=w.sha(s/'delivery_manifest.tsv'),scientific_gate_passed=False,all_owned_parents_terminal=True))
 (s/'COMPLETE').write_text(w.sha(s/'delivery_manifest.tsv')+'\nEXECUTION_CLOSED_AT_REDUCTION_HOLD\nSCIENTIFIC_GATE_PASSED=false\n')
 print('CLOSED_REDUCTION_HOLD',len(rows),total,'core_seconds',len(manifest),'scientific files',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);p.add_argument('release',type=Path);a=p.parse_args();close(a.study,a.release)
