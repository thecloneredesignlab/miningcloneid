"""Reread saved enclosures and close finite source-based algebraic obligations."""
from pathlib import Path
import sys
from common import read,table,sha,roles_valid,MECHANISMS,NUISANCES,matrix_records
from runtime import boot

def main():
    study=Path(sys.argv[1]);task=sys.argv[2]
    sp,flint,np,out=boot(study,task)
    from algebra import decode,encode,independent_minor_check
    import re
    assert task=='P02'
    roles_valid(read(study/'inputs/parameter_roles.tsv'))
    flint.ctx.prec=256;checks=[];reviews=[];summaries=[]
    def check(name,condition,detail=''):
        checks.append(dict(check=name,passed=bool(condition),detail=detail));table(out/'review_checks.tsv',checks)
        assert condition,name
    for case,n in (('B01',16),('B04',18)):
        for witness in ('W1','W2','W3'):
            d=study/'inputs/evidence/verification'/(case+'_'+witness)
            result=read(d/'case_result.tsv')[0]
            check('original_rank_'+case+witness,int(result['rank_lower_bound'])==n)
            a=[r for r in read(d/'certificate_attempts.tsv') if r['status']=='CERTIFIED'];check('one_positive_'+case+witness,len(a)==1)
            stem=a[0]['candidate']+'_'+a[0]['precision']
            candidate=next(r for r in read(d/'candidate_registry.tsv') if r['candidate']==a[0]['candidate'])
            ids=[int(v) for v in candidate['rows'].split(',')];names=candidate['columns'].split(',')
            registry={int(r['index']):r for r in read(d/'row_registry.tsv')}
            cols={r['parameter'] for r in read(d/'column_registry.tsv')}
            check('candidate_identities_'+case+witness,len(ids)==n and len(set(ids))==n and len(names)==n and len(set(names))==n and set(ids)<=set(registry) and set(names)<=cols)
            expected=set(MECHANISMS+NUISANCES[:4])-({'O2_crit','n_O'} if case=='B01' else set())
            check('original_coordinate_columns_'+case+witness,set(names)==expected)
            detrow=read(d/(stem+'_determinant.tsv'))[0]
            check('saved_determinant_radius_'+case+witness,int(detrow['rad_man'])>=0)
            determinant=decode(detrow)
            check('saved_determinant_nonzero_'+case+witness,determinant.is_finite() and not determinant.contains(0))
            matrices=[]
            for route in ('matrix','reference'):
                path=d/(stem+'_'+route+'.tsv');records=read(path);matrix_records(records,n)
                M=flint.arb_mat(n,n)
                for row in records:
                    v=decode(row);assert v.is_finite(),'Nonfinite saved matrix entry'
                    M[int(row['row']),int(row['column'])]=v
                check('finite_entries_'+case+witness+route,True)
                verified,bound=independent_minor_check(M)
                check('Neumann_'+case+witness+route,verified and bound.is_finite() and bound<1,str(bound))
                reviews.append(dict(case=case,witness=witness,route=route,rank=n,precision=256,status='PASS',
                    matrix_sha256=sha(path),neumann_upper=str(bound),**encode(bound),new_determinant_evaluations=0))
                matrices.append(M)
            check('independent_routes_overlap_'+case+witness,all(matrices[0][i,j].overlaps(matrices[1][i,j]) for i in range(n) for j in range(n)))
            summaries.append(dict(case=case,witness=witness,status='PASS',rank=n,saved_determinant_excludes_zero=True,matrices_reviewed=2,
                proof_scope='original continuous finite-support fixed-mask local certificate; derivative/source proofs reused, not regenerated'))
    table(out/'saved_certificate_review.tsv',reviews);table(out/'saved_certificate_inventory.tsv',summaries)
    table(out/'saved_certificate_gate.tsv',[dict(status='PASS',certificates=6,matrix_reviews=12,precision=256,new_project_cases=0,new_determinants=0)])

    # Exact scalar identities plus a manually source-reviewed compositional route.
    p={name:sp.Symbol(name,positive=True) for name in MECHANISMS}
    O,u,Cnew,nnew=sp.symbols('O u Cnew nnew',positive=True)
    H=p['O2_crit']**p['n_O']/(p['O2_crit']**p['n_O']+O**p['n_O'])
    Hnew=Cnew**nnew/(Cnew**nnew+O**nnew)
    lam=p['lam_max']/(1+p['alpha_o2']*H*u**p['gamma_growth'])
    mu=p['mu_hp']*H*u**p['gamma_mu']
    prob=p['p_mis_base']+p['p_misseg']*mu/(p['k_o_mis']+mu)
    logs=sp.log(p['buffer_smax'])-p['buffer_beta']*(1/u)**p['buffer_n_exp']
    obligations=[]
    def zero(name,expr,scope):
        residual=sp.simplify(expr);check(name,residual==0,str(residual))
        obligations.append(dict(obligation=name,residual=str(residual),status='PASS_EXACT',scope=scope))
    zero('induction_zero_retains_base',prob.subs(p['p_misseg'],0)-p['p_mis_base'],'p_misseg=0')
    zero('induction_zero_inactive_k',sp.diff(prob.subs(p['p_misseg'],0),p['k_o_mis']),'p_misseg=0')
    zero('survival_zero_retains_smax',logs.subs(p['buffer_beta'],0)-sp.log(p['buffer_smax']),'buffer_beta=0')
    zero('survival_zero_inactive_nexp',sp.diff(logs.subs(p['buffer_beta'],0),p['buffer_n_exp']),'buffer_beta=0')
    m=sp.Symbol('m',integer=True,positive=True)
    zero('survival_error_size_retained',sp.exp(m*logs.subs(p['buffer_beta'],0))-p['buffer_smax']**m,'m>0; buffer_beta=0')
    zero('zero_error_survival',sp.exp(0*logs)-1,'m=0')
    zero('zero_O_Hill',H.subs(O,0)-1,'O=0; positive O2_crit,n_O')
    for name in ('O2_crit','n_O'):zero('zero_O_inactive_'+name,sp.diff(H.subs(O,0),p[name]),'O=0 only')
    changed={p['O2_crit']:Cnew,p['n_O']:nnew,p['alpha_o2']:p['alpha_o2']*H/Hnew,p['mu_hp']:p['mu_hp']*H/Hnew}
    for name,v in (('lambda',lam),('mu',mu),('probability',prob),('log_survival',logs)):
        zero('single_O_family_'+name,v.subs(changed,simultaneous=True)-v,'same positive O; nearby admissible transformed coordinates')
    dt,a,Ostate,Otarget=sp.symbols('dt relaxation Ostate Otarget')
    zero('constant_O_relaxation',(Ostate+a*(Otarget-Ostate)).subs({Ostate:O,Otarget:O})-O,'feedback disabled; equal initial and target O')
    # Equal arguments of a deterministic clamp/select/reset function give equal
    # outputs, including ties and infeasibility; no numerical trajectory is used.
    x,gx,death,b,c,Hstock,Bstock,k=sp.symbols('x Gx death b c Hstock Bstock k')
    z=sp.Symbol('same_value')
    for name,expression in (('live_step',x+dt*c*gx-dt*death*x),
        ('direct_death_stock',Hstock+dt*death*x-dt*k*Hstock),
        ('division_loss_stock',Bstock+dt*c*b*x-dt*k*Bstock)):
        zero('equal_arguments_'+name,expression.subs(x,z)-expression.subs(x,z),'inductive identical-state/operator premise; componentwise clamp follows by congruence')
    table(out/'boundary_proof_obligations.tsv',obligations)

    cpp=study/'inputs/production/model_O2_supply_demand_MAP.cpp'
    lineage=study/'inputs/production/o2_supply_demand_map_invitro_lineage_simulation_utils.R'
    cpptext=cpp.read_text();rtext=lineage.read_text()
    mapping=[
      ('hypoxia_weight_cpp','H=C^n/(C^n+O^n); positive declared C/n; inactive invalid-input fallbacks excluded'),
      ('lambda_eff_runtime_cpp','O2_growth TRUE dispatch to the audited oxygen-damped lambda; no hidden k or survival exponent'),
      ('mu_eff_soft_cpp','actual ploidy_related branch; mu=mu_hp*H*(N/44)^gamma_mu'),
      ('death_rate_for_N_cpp','direct-death dispatcher returns the same audited mu_eff_soft_cpp with actual canonical mode'),
      ('resolve_pmis_for_death','only semantic use of k_o_mis is in p_base+p_amp*mu/(k+mu); output clamp is deterministic'),
      ('binom_prob_int','pi depends on N,m,p only'),
      ('buffering_survival_modifier','m=0 returns1; m>0 uses (smax*exp(-beta*(44/N)^nexp))^m'),
      ('o2simps_pr_delta_internal','keep pi>=eps; two signed daughters; m=0 added twice; survival affects daughter weights'),
      ('finalize_pr_delta_internal','pre-boundary loss2-surviving daughter sum includes viability and pruned surviving mass'),
      ('append_block_with_boundary','drop out-of-support daughters and route their weights into buffer'),
      ('cpp_o2simps_build_G_for_o2_triplet','ordinary lambda*(1-w), diagonal -lambda, single WGD lambda*w; direct death remains separate'),
      ('g_cache_signature_cpp','signature contains semantic model coordinates; changed signature clears cache; cache correctness assumed without hash collision'),
      ('cpp_o2simps_simulate_one','fixed oxygen initialization; Euler/crowding/death/clearance/clamps; no simulator call is made here')]
    source_map=[]
    for name,argument in mapping:
        matches=list(re.finditer(r'^(?:inline\s+)?(?:double|void|List|std::size_t)\s+'+re.escape(name)+r'\s*\(',cpptext,re.M))
        check('source_definition_'+name,len(matches)==1)
        line=cpptext.count('\n',0,matches[0].start())+1
        source_map.append(dict(file=cpp.name,source_sha256=sha(cpp),function=name,line=line,reviewed_dependency_argument=argument,proof_route='source-reviewed composition plus exact response identities; token presence alone is not proof'))
    for name,argument in (
      ('ivt_set_segment_o2','feedback FALSE; segment target overwrites S0 and o2_min'),
      ('ivt_gaussian_initial_state','root mean/SD unchanged by each claimed parameter family'),
      ('ivt_run_segment_fixed_o2','parent live override dominates Gaussian; live normalized to known inoculum; no retained-state overrides are supplied'),
      ('ivt_select_segment_observation','positive-day candidate subset; nearest observed target/first-index tie; missing target last candidate'),
      ('ivt_extract_passage_end_state','eligible fallback smallest excess then index; insufficient supply returns failure; fractions preserved under downsampling'),
      ('ivt_run_lineage','parent-child extraction; protocol failure propagation; inherited fraction fallback; recorded final count then next-inoculum conditioning')):
        pos=rtext.find(name+' <-');check('source_definition_'+name,pos>=0)
        source_map.append(dict(file=lineage.name,source_sha256=sha(lineage),function=name,line=rtext.count('\n',0,pos)+1,reviewed_dependency_argument=argument,proof_route='source-reviewed control-flow argument; pure selectors additionally tested in SIF P03'))
    table(out/'source_dependency_map.tsv',source_map)
    # Verify crucial actual switches from their frozen cfg/run inventory.
    settings=read(study/'inputs/fixed_settings.tsv')
    required={'DT':'0.05','o2_burden_feedback':'FALSE','O2_growth':'TRUE','Crowding':'TRUE','crowding':'logistic','K':'1e+12','ploidy_O2_death':'ploidy_related'}
    for key,value in required.items():
        matches=[r for r in settings if r['setting']==key and r['scope']=='archived_cfg']
        check('actual_switch_'+key,len(matches)==1 and matches[0]['value']==value)
    table(out/'boundary_activity.tsv',[
        dict(boundary='p_misseg=0',inactive='k_o_mis',status='PROVED_FORWARD_DEPENDENCE',retained='p_mis_base; other mechanisms',scope='valid nonnegative domain and finite real arithmetic; not a penalty/prior or boundary rank statement'),
        dict(boundary='buffer_beta=0',inactive='buffer_n_exp',status='PROVED_FORWARD_DEPENDENCE',retained='smax^m for m>0 and1 for m=0',scope='valid survival domain; not a boundary rank statement'),
        dict(boundary='O=0 only',inactive='O2_crit;n_O',status='PROVED_FORWARD_DEPENDENCE',retained='alpha_o2;mu_hp; other mechanisms',scope='positive shape coordinates; all relevant inputs at zero O')])
    (out/'boundary_proofs.md').write_text('''# Boundary compositional proofs

The exact identities in boundary_proof_obligations.tsv and the source-reviewed
dependency map jointly support the following finite statements. Source string
matches identify locations; they do not establish the statements alone.

At p_misseg=0 the saturating increment is zero for every mother N. Therefore p
is independent of k_o_mis. Every binomial pi, keep predicate, signed daughter
weight and live/drop event is unchanged when only k changes. Lambda, mu, WGD,
survival and initial/observer maps have no other semantic dependence on k.
Retained inflows and subsequent deterministic steps/resets are compositions of
these same arguments. Basal missegregation and all other mechanisms remain.

At buffer_beta=0 the per-copy survival is smax, so m>0 survival is smax^m;
m=0 remains exactly1. No other active response depends on buffer_n_exp. Thus
all live events, nonviability/pruning/boundary accounts and downstream maps
are exponent-independent. Error-size dependence is retained through m.

At O=0 with positive C and n, H=1. Growth, death and induced probability no
longer depend on C,n. Initial and assay settings are unchanged. The same event
composition applies. This assertion is specific to zero-input experiments.

For the event account, let A be the sum of kept surviving signed daughters
before support drop and S their in-support sum. Pre-boundary buffer is
lambda*(1-w)*(2-A). Support drop adds lambda*(1-w)*(A-S), plus off-grid WGD
lambda*w. Their sum is lambda*(1-w)*(2-S)+offgrid_WGD. Biological mortality,
pruned surviving probability and lower/upper support losses remain separate
concepts; the buffer is not pure biological death. This uses accepted B2 event
evidence; no kernel or event table was regenerated.

Assumptions are the original admissible finite parameter domain, N=22:154,
fixed drop support, positive denominator/widths and the stated production
branches. No generic interior rank is transferred to a null boundary. The
map excludes priors/penalties. All14 mechanisms and6 nuisances remain in the
inventory; variance-only sigma_growth is not promoted to a mean parameter.
''')
    (out/'discrete_induction.md').write_text('''# Constant-oxygen discrete/passage induction

This is a source-reviewed real-arithmetic argument, conditional on the pure
selection/reset contract tests in P03. It is not an executed trajectory or a
positive discrete identifiability certificate.

For one positive constant O, change C,n to nearby admissible C',n' and set
alpha'=alpha*H(O;C,n)/H(O;C',n'), mu_hp'=mu_hp*H(O;C,n)/H(O;C',n'). Exact
response identities prove equal lambda, mu, p and survival for every N. Initial
families, observer settings and all external inputs remain unchanged.

Feedback is FALSE. Segment setup makes S0 equal the target; the archived cap21
contains the declared O values. Initial O_state and O_target equal O; relaxation
O_state+a*(O_target-O_state) preserves O. Treatment is absent. Thus each step
uses equal rates, masks, G and per-parent buffer rates in both parameterizations.

Let x_N,H_N,B_N be indexed live/hypoxia/buffer states; c=max(0,1-sum(x)/K).
The source update is componentwise clamp of x+DT*c*Gx-DT*mu*x,
H+DT*mu*x-DT*k_clear*H, and B+DT*c*b*x-DT*k_clear*B. Direct death is outside
the division crowding multiplier; buffer inflow is division-linked. Equal
indexed states and operators yield equal unclamped arguments, hence equal
clamps, observations and any deterministic extinction/observation bookkeeping.
Induction proves equality at every real-arithmetic step for the same DT=0.05.

The same daily arrays and external targets give the same positive-day candidate
set, nearest-target ordering and first-index tie break. With no target the last
candidate is selected. If target-selected supply is below the requested next
inoculum, the actual extraction function chooses an eligible state with smallest
excess supply then first index, or returns the same protocol failure. The same
selected fractions are multiplied by the same known inoculum. Parent overrides
the root Gaussian in inherited segments. The wrapper supplies no retained-state
overrides, so both retained compartments reset to their zero defaults at segment
initialization; this is not an assumption that dead stocks are inherited.
The identical-selection/rescaling/reset premise restarts the induction in the
next segment. Source fallback to prior fractions and invalid-root/config fallback
is deterministic under the unchanged inputs; claims use admissible root widths.

Caches are treated as correct memoization for a run restricted to the same
constant O, with a fresh or valid same-O cache, no signature collision and no
alias from a previous distinct O in the same cache bin. Source signatures cover
model inputs and clear the cache on signature change. These assumptions restrict
the mathematical claim: literal hash collisions, floating-point rounding,
underflow and cache/bin bit patterns are not proved equivalent. No multi-O
history, oxygen feedback or positive discrete identifiability result follows.

Observed final counts (or the next-inoculum fallback) used to select a day are
conditioning inputs even if growth is not scored. This is not independent
prediction validation. Observed OLS growth versus scored endpoint log-ratio,
physical flow calibration and eligible live/dead denominators remain unresolved
experimental/design issues. Step5 remains unreleased.
''')
    table(out/'continuous_discrete_scope.tsv',[
        dict(claim='single_O_real_arithmetic_D_and_passage_invariance',status='CONDITIONAL_COMPOSITION_REVIEWED_PENDING_P03',limitation='constant O; no feedback; valid same-O cache; same external inputs; P03 selector contracts required'),
        dict(claim='positive_discrete_identifiability',status='NOT_ANALYZED',limitation='continuous certificate does not transfer'),
        dict(claim='literal_floating_point_or_multioxygen_history',status='NOT_CLAIMED',limitation='rounding/underflow/cache aliases excluded'),
        dict(claim='actual_sparse_data_or_variance_nuisance_identifiability',status='NOT_ESTABLISHED',limitation='experimental contracts and practical recovery outside scope')])
    table(out/'boundary_gate.tsv',[dict(status='PASS',boundary_arguments=3,exact_identities=len(obligations),source_routes=len(source_map),passage_contract_confirmation='P03_required')])
    table(out/'application_result.tsv',[dict(status='PASS',task=task,phase='saved_evidence_and_boundary_review',new_project_cases=0)])
    print('SAVED_CERTIFICATES_AND_BOUNDARY_COMPOSITION_PASS',len(reviews),len(obligations),flush=True)

if __name__=='__main__':main()
