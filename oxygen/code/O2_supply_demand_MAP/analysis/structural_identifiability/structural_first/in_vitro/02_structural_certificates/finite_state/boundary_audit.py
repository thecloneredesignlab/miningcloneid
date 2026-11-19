"""Finite compositional boundary proofs and a source-based discrete scope audit."""
from pathlib import Path
import sys
from runtime import boot
from protocol import read, table, sha, MECHANISMS, NUISANCES, roles_valid, kv
from algebra import decode, independent_minor_check

study=Path(sys.argv[1]); task=sys.argv[2]
assert task=='B6'
sp,flint,np,out=boot(study,task)
roles_valid(read(study/'inputs/parameter_roles.tsv'))
accepted=kv(study/'accounting/operator_acceptance.tsv')['phase']
assert read(study/'verification'/accepted/'operator_gate.tsv')[0]['status']=='PASS'
assert read(study/'verification/B1/validation_gate.tsv')[0]['status']=='PASS'
checks=[]
def check(name,condition,detail=''):
    checks.append(dict(check=name,passed=bool(condition),detail=detail))
    table(out/'boundary_checks.tsv',checks)
    assert condition,name
def zero(name,expr):
    residual=sp.simplify(expr)
    check(name,residual==0,str(residual))

p={name:sp.Symbol(name,positive=True) for name in MECHANISMS}
u,O=sp.symbols('parent_ratio oxygen',positive=True)
H=p['O2_crit']**p['n_O']/(p['O2_crit']**p['n_O']+O**p['n_O'])
mu=p['mu_hp']*H*u**p['gamma_mu']
prob=p['p_mis_base']+p['p_misseg']*mu/(p['k_o_mis']+mu)
logs=sp.log(p['buffer_smax'])-p['buffer_beta']*(1/u)**p['buffer_n_exp']
zero('exact_induction_zero_retains_basal_probability',prob.subs(p['p_misseg'],0)-p['p_mis_base'])
zero('exact_induction_zero_removes_k',sp.diff(prob.subs(p['p_misseg'],0),p['k_o_mis']))
zero('exact_survival_zero_retains_error_size_dependence',logs.subs(p['buffer_beta'],0)-sp.log(p['buffer_smax']))
zero('exact_survival_zero_removes_exponent',sp.diff(logs.subs(p['buffer_beta'],0),p['buffer_n_exp']))
zero('zero_oxygen_positive_shape_Hill',H.subs(O,0)-1)
zero('zero_oxygen_inactive_critical',sp.diff(H.subs(O,0),p['O2_crit']))
zero('zero_oxygen_inactive_exponent',sp.diff(H.subs(O,0),p['n_O']))
m=sp.Symbol('error_size',integer=True,positive=True)
zero('nonzero_error_survival_keeps_smax',sp.exp(m*logs.subs(p['buffer_beta'],0))-p['buffer_smax']**m)
zero('zero_error_survival_is_one',sp.exp(0*logs)-1)
check('retained_original_inventory',len(MECHANISMS)==14 and len(NUISANCES)==6)

source=study/'inputs/production'
cpp=(source/'model_O2_supply_demand_MAP.cpp').read_text()
lineage=(source/'o2_supply_demand_map_invitro_lineage_simulation_utils.R').read_text()
requirements={
 'source_logistic_current_live':'1.0 - (Ntot_live_now / K_use)',
 'source_division_scalar':'DT_use * crowd_mult * tx_mult',
 'source_direct_death_outside_scalar':'DT_use * mu_i * src_live',
 'source_buffer_division_linked':'scalar * db_rate_i * src_live',
 'source_live_update':'v_live[i] + scalar * growth[i] - death_flow_hypoxia[i]',
 'source_hypoxia_clearance':'dead_h_prev + death_flow_hypoxia[i] - DT_use * k_clear_use * dead_h_prev',
 'source_buffer_clearance':'dead_b_prev + death_flow_buffer[i] - DT_use * k_clear_use * dead_b_prev',
 'source_nonnegative_clamp':'!std::isfinite(next_v) || next_v < 0.0',
}
for name,text in requirements.items():check(name,text in cpp,text)
lineage_requirements={
 'target_first_index_tie':'order(abs(live_cells[candidate_idx] - target_live_cells_use), candidate_idx)',
 'positive_days_first':'positive_day_idx <- which(is.finite(obs_days_use) & obs_days_use > 0)',
 'missing_target_last_candidate':'candidate_idx[[length(candidate_idx)]]',
 'inoculum_eligible_fallback':'live_state_totals[eligible_idx] - reseed_live_cells',
 'inoculum_failure':'required_inoculum_not_reached',
 'inherited_parent_fraction':'as.numeric(parent_boundary_selection$reseeded_state)',
 'lineage_infeasible_condition':'invitro_protocol_infeasible',
 'recorded_final_count_conditioning':'target_live_cells_use <- as.numeric(seg$final_cells)',
 'fallback_next_inoculum':'target_live_cells_use <- as.numeric(next_seg$initial_cells)',
}
for name,text in lineage_requirements.items():check(name,text in lineage,text)
table(out/'source_protocol_audit.tsv',[
 dict(file=p.name,sha256=sha(p),scope='static source correspondence; no lineage or propagation call')
 for p in (source/'model_O2_supply_demand_MAP.cpp',source/'o2_supply_demand_map_invitro_lineage_simulation_utils.R')])

# A compositional proof: the kernel depends on lambda, mu, p, logS, WGD and
# fixed support/masks; the initial and observation interfaces do not depend on
# the coordinates removed below. No objective/penalty is included in this map.
table(out/'boundary_activity.tsv',[
 dict(boundary='p_misseg=0',inactive='k_o_mis',retained='p_mis_base; WGD; growth; death; survival',
      proof='p=p_mis_base; every binomial weight and pruning predicate is k-independent; all later forward maps compose these values',
      scope='exact zero; original nonnegative domain; no penalty/prior claim'),
 dict(boundary='buffer_beta=0',inactive='buffer_n_exp',retained='buffer_smax and error-size exponent m',
      proof='m>0 survival=smax^m, m=0 survival=1; no other forward dependence on buffer_n_exp',
      scope='exact zero; original nonnegative domain; no penalty/prior claim'),
 dict(boundary='O=0 only',inactive='O2_crit;n_O',retained='alpha_o2;mu_hp;other response coordinates',
      proof='H=1 for positive O2_crit and n_O; every rate, mask and downstream forward map is shape-independent',
      scope='all segments at zero oxygen; not a positive/multiple-oxygen result')])

table(out/'continuous_discrete_scope.tsv',[
 dict(item='continuous_positive_local_certificates',status='CONDITIONAL_ONLY',
      evidence='case-specific Arb minors and independently reconstructed Neumann certificates',
      limitation='133 live states; fixed drop support; strict fixed pruning masks; ideal continuous observations; 0<L<K'),
 dict(item='single_O_discrete_real_arithmetic_symmetry',status='PROVED_BY_COMPOSITION_AND_INDUCTION',
      evidence='B1 exact rate invariance plus B2 event identity: identical initial state and lambda/mu/p/logS imply identical G, masks and indexed stock inflows; equal states give equal updates and clamps at every step',
      limitation='same one constant O across all segments; identical known DT=0.05, inputs, seeds, assay settings, targets and reseeding counts; no oxygen feedback'),
 dict(item='passage_selection_induction',status='PROVED_UNDER_IDENTICAL_INPUTS',
      evidence='identical daily states/counts imply equal nearest-target ordering, first-index ties, eligible fallback, failures, selected fractions and reseeded states; induction continues to the next segment',
      limitation='target observed final counts remain conditioning information; inherited fractions are not reset to a new Gaussian'),
 dict(item='multiple_O_production_lineages',status='NOT_PROVED',
      evidence='single-O compensating transform uses H(O)/H_new(O)',
      limitation='one transform is not asserted to preserve several distinct oxygen levels'),
 dict(item='positive_discrete_identifiability',status='NOT_ANALYZED',
      evidence='no discrete sensitivity or hybrid solver; no trajectory replay',
      limitation='continuous full-rank jets do not prove discrete/passage identifiability'),
 dict(item='literal_floating_point_equality',status='NOT_CLAIMED',
      evidence='real-arithmetic source composition only',
      limitation='cache keys, numerical rounding, underflow and saturation/clamp boundaries can differ'),
 dict(item='actual_growth_flow_live_dead_readiness',status='HELD_EXPERIMENTAL_CONTRACT',
      evidence='reuse Steps1-3 and Stage A audit without replacing actual records by ideal outputs',
      limitation='observed OLS growth versus scored endpoint log-ratio; physical flow calibration; eligible live/dead denominators remain separate unresolved design/experimental issues')])

# Verify saved exact dyadic enclosures without another determinant evaluation.
# Both the primary and the independent event-route reconstruction must pass.
flint.ctx.prec=256
reviews=[]
for wid in ('W1','W2','W3'):
    for case,expected in (('B01',16),('B04',18)):
        directory=study/'verification'/(case+'_'+wid)
        result=read(directory/'case_result.tsv')[0]
        check('saved_certificate_rank_'+case+'_'+wid,int(result['rank_lower_bound'])==expected)
        successes=[r for r in read(directory/'certificate_attempts.tsv') if r['status']=='CERTIFIED']
        check('one_saved_certificate_'+case+'_'+wid,len(successes)==1)
        a=successes[0];prefix=a['candidate']+'_'+a['precision']
        determinant=decode(read(directory/(prefix+'_determinant.tsv'))[0])
        check('saved_determinant_excludes_zero_'+case+'_'+wid,determinant.is_finite() and not determinant.contains(0))
        for suffix in ('matrix','reference'):
            path=directory/(prefix+'_'+suffix+'.tsv');entries=read(path)
            check('saved_matrix_shape_'+case+'_'+wid+'_'+suffix,len(entries)==expected**2)
            seen=set();matrix=flint.arb_mat(expected,expected)
            for row in entries:
                i,j=int(row['row']),int(row['column'])
                check_name='saved_entry_identity_'+case+'_'+wid+'_'+suffix
                assert 0<=i<expected and 0<=j<expected and (i,j) not in seen,check_name
                seen.add((i,j));matrix[i,j]=decode(row)
            good,bound=independent_minor_check(matrix)
            check('saved_neumann_'+case+'_'+wid+'_'+suffix,good,str(bound))
            reviews.append(dict(case=case,witness=wid,route=suffix,precision=256,
                sha256=sha(path),rank=expected,status='PASS',neumann_bound=str(bound),
                new_determinant_evaluations=0))
table(out/'saved_certificate_review.tsv',reviews)
table(out/'boundary_gate.tsv',[dict(status='PASS',checks=len(checks),native_calls=0,
    forward_trajectory_calls=0,original_mechanisms=14,original_nuisances=6,
    discrete_positive_certificate=False,Step_5=False,fits=0)])
table(out/'application_result.tsv',[dict(status='PASS',task=task,phase='B6')])
print('FINITE_BOUNDARY_AND_DISCRETE_SCOPE_AUDIT_PASSED',len(checks),flush=True)
