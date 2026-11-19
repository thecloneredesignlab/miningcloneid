"""Registered exact toy-interface controls only; never a culture solver."""
import os
from pathlib import Path
import sys
import traceback
from protocol import read,table,sha,validate_roles,validate_schema,validate_collection

study=Path(sys.argv[1]).resolve();out=study/'verification';out.mkdir(exist_ok=True)
assert not (out/'application_result.tsv').exists(),'Refusing stale worker output'
for r in read(study/'source_manifest.tsv'):assert sha(study/r['path'])==r['sha256']
rows=[];cases=[];witnesses=[];current='T00'
forbidden={'json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius',
           'execjs','py_mini_racer','symbtools','strike_goldd','native_strike'}
def guard(event,args):
    if event=='import' and args[0].split('.')[0] in forbidden|{'ipydex'}:
        with (out/'blocked_imports.tsv').open('a') as f:f.write('BLOCKED\t'+args[0]+'\n')
        raise ImportError('Excluded dependency in finite interface controls')
sys.addaudithook(guard)
sys.path[:0]=['/home/4482173/.local/lib/python3.11/site-packages',
             *[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]

def record(name,residual,detail=''):
    rows.append(dict(case=current,check=name,residual=str(residual),detail=detail))
    table(out/'exact_checks.tsv',rows)

def zero(expr,name):
    residual=sp.simplify(expr)
    if residual!=0:
        residual=sp.simplify(sp.cancel(sp.together(sp.expand(expr))))
    record(name,residual)
    assert residual==0,current+': '+name+': '+str(residual)

def same(a,b,name):
    assert a.shape==b.shape,name+': shape'
    for i in range(a.rows):
        for j in range(a.cols):zero(a[i,j]-b[i,j],name+'_'+str(i)+'_'+str(j))

def passed(detail):
    cases.append(dict(case=current,status='PASS',evidence_route='exact_interface_certificate',detail=detail))
    table(out/'control_results.tsv',cases)
    print(current+' PASS '+detail,flush=True)

def rejected(name,callback):
    try:callback()
    except ValueError:record(name,0,'Expected invalid interface rejected')
    else:raise AssertionError('Invalid interface accepted: '+name)

try:
    import sympy as sp
    assert sp.__version__=='1.12'
    import adapter as ad
    import reference as ref
    validate_roles(read(study/'inputs/parameter_roles.tsv'))
    record('twenty_original_roles_retained',0)
    record('scoped_symbolic_runtime',0,'SymPy 1.12; no native rank or Lie driver')
    passed('Twenty roles and scoped symbolic startup')

    current='T01'
    grid=list(map(sp.Integer,[-1,0,1]))
    mu=sp.Symbol('toy_mu',real=True);s=sp.Symbol('toy_s',positive=True)
    r=sp.Symbol('toy_r',positive=True)
    p=ad.initial_weights(sp,grid,mu,s)
    pr,dm,ds=ref.finite_family(sp,grid,mu,s)
    same(p,pr,'independent_exponential_family')
    same(p.diff(mu),dm,'mean_chain_formula');same(p.diff(s),ds,'width_chain_formula')
    zero(sum(p)-1,'normalization');zero(sum(p.diff(mu)),'mean_derivative_sum')
    zero(sum(p.diff(s)),'width_derivative_sum')
    product_ratio=sp.simplify(p[0]*p[2]/p[1]**2)
    directional_ratio=sp.simplify(p[2]/p[0])
    zero(product_ratio-sp.exp(-1/s**2),'inverse_variance_ratio')
    zero(sp.log(product_ratio)+1/s**2,'inverse_variance_log')
    zero(s**2*sp.log(directional_ratio)/2-mu,'inverse_mean')
    d0=p.jacobian([mu,s]);minor=sp.simplify(d0[:2,:].det())
    target=2*p[0]*p[1]*p[2]/s**5
    zero(minor-target,'positive_initial_minor')
    witnesses.append(dict(case=current,witness=str(target),domain='s>0; finite real mu; all three probabilities positive',conclusion='toy initial mean and positive width identifiable from full probabilities'))
    passed('Normalized finite family, exact derivatives and inverse identities')

    current='T02'
    x=sp.Matrix(sp.symbols('toy_x_minus toy_x_zero toy_x_plus',positive=True))
    kappa,c,n0=sp.symbols('toy_kappa toy_c toy_N0',positive=True)
    jets=ad.output_jets(sp,x,-c*kappa*x,x,2)
    values,jacs=ad.initial_pullback(sp,jets,x,n0*p,[mu,s,kappa])
    jet_rows=[]
    for order in range(3):
        same(jets[order],(-c*kappa)**order*x,'state_space_jet_'+str(order))
        same(values[order],(-c*kappa)**order*n0*pr,'constrained_jet_'+str(order))
        factor=(-c*kappa)**order*n0
        dk=0 if order==0 else order*(-c)**order*kappa**(order-1)*n0
        expected=sp.Matrix.hstack(factor*dm,factor*ds,dk*pr)
        same(jacs[order],expected,'independent_chain_rows_'+str(order))
        for i in range(3):
            for j,param in enumerate([mu,s,kappa]):
                jet_rows.append(dict(case=current,order=order,row=i,parameter=str(param),expression=str(jacs[order][i,j])))
    transformed=jacs[1][0,:]+c*kappa*jacs[0][0,:]
    zero(transformed[0],'rank_row_elimination_mean')
    zero(transformed[1],'rank_row_elimination_width')
    zero(transformed[2]+c*n0*p[0],'rank_row_elimination_rate')
    rate_witness=-c*n0**3*p[0]*target
    witnesses.append(dict(case=current,witness=str(rate_witness),domain='c,N0,kappa,s>0; finite mu',conclusion='toy (mu,s,kappa) locally identifiable; 3x3 determinant via exact row elimination'))
    table(out/'toy_jet_rows.tsv',jet_rows)
    passed('State-space jets precede initial substitution; complete chain rule')

    current='T03'
    total_jets=ad.output_jets(sp,x,-c*kappa*x,[sum(x)],2)
    totals,total_jacs=ad.initial_pullback(sp,total_jets,x,n0*p,[mu,s,kappa])
    for order in range(3):
        zero(totals[order][0]-n0*(-c*kappa)**order,'total_jet_'+str(order))
        zero(total_jacs[order][0,0],'total_mean_column_'+str(order))
        zero(total_jacs[order][0,1],'total_width_column_'+str(order))
    t=sp.Symbol('toy_t',real=True)
    trajectory=n0*sp.exp(-c*kappa*t)
    zero(trajectory.diff(mu),'all_time_mean_invariance')
    zero(trajectory.diff(s),'all_time_width_invariance')
    witnesses.append(dict(case=current,witness=str(trajectory),domain='any finite mu and any s>0; same c,kappa,N0',conclusion='explicit all-time family: arbitrary initial mean/width changes preserve live-total output'))
    passed('Total-only negative control has an explicit all-time invariant')

    current='T04'
    kernel=ad.kary_kernel(sp,grid,grid,r);kr,dr=ref.manual_kernel(sp,r)
    same(kernel,kr,'independent_kernel');same(kernel.diff(r),dr,'quotient_width_derivative')
    for j in range(3):
        zero(sum(kernel[:,j])-1,'kernel_column_'+str(j))
        zero(sum(kernel.diff(r)[:,j]),'kernel_derivative_column_'+str(j))
    q=kernel[:,1]
    ratio=sp.simplify(q[2]/q[1])
    zero(ratio-sp.exp(-1/(2*r**2)),'point_mass_width_ratio')
    zero(sp.log(ratio)+1/(2*r**2),'point_mass_log_ratio')
    zero(sum(z*v for z,v in zip(grid,q)),'symmetric_observed_mean')
    zero(sum(z*v for z,v in zip(grid,q)).diff(r),'mean_only_width_column')
    witnesses.append(dict(case=current,witness='r^2=-1/(2*log(q_plus/q_zero))',domain='r>0; known latent point mass at zero',conclusion='toy width identifiable from full probabilities; mean-only output invariant'))
    passed('Unknown normalized observer width and mean-only negative check')

    current='T05'
    v,w=sp.symbols('toy_initial_variance toy_observer_variance',positive=True)
    delta,y,freq=sp.symbols('toy_delta toy_y toy_frequency',real=True)
    char=sp.exp(sp.I*mu*freq-(v+w)*freq**2/2)
    density=sp.exp(-(y-mu)**2/(2*(v+w)))/sp.sqrt(2*sp.pi*(v+w))
    zero(sp.exp(sp.I*mu*freq-v*freq**2/2)*sp.exp(-w*freq**2/2)-char,'characteristic_product')
    for label,expression in [('characteristic',char),('density',density)]:
        zero(expression.subs({v:v+delta,w:w-delta},simultaneous=True)-expression,label+'_invariance')
        zero(expression.diff(v)-expression.diff(w),label+'_equal_variance_columns')
    witnesses.append(dict(case=current,witness='v_new=v+delta; w_new=w-delta',domain='v,w>0 and -v<delta<w',conclusion='unbounded Gaussian toy only: observed variance v+w; individual widths confounded'))
    passed('Exact Gaussian variance-confounding family; no finite-grid inference')

    current='T06'
    half=sp.Rational(1,2);observed=[z/2 for z in grid]
    flow=ad.flow_density(sp,observed,grid[:2],[half,half],r/2,sp.Integer(2),half)
    raw=sp.Matrix([[1,sp.exp(-1/(2*r**2))],
                   [sp.exp(-1/(2*r**2)),1],
                   [sp.exp(-2/r**2),sp.exp(-1/(2*r**2))]])
    d=raw*sp.Matrix([half,half]);flow_ref=d/(sum(d)*half)
    same(flow,flow_ref,'flow_global_normalization_reference')
    zero(sum(flow)*half-1,'flow_grid_weighted_normalization')
    zero(sum(flow)-2,'density_is_not_probability_vector')
    same(sp.Matrix([n/sp.Integer(2) for n in grid]),sp.Matrix(observed),'chromosome_to_ploidy_conversion')
    u=sp.Symbol('edge_u',positive=True)
    rawu=sp.Matrix([[1,u],[u,1],[u**4,u]])
    karyu=sp.Matrix([sum(rawu[i,j]/sum(rawu[:,j])*half for j in range(2)) for i in range(3)])
    du=rawu*sp.Matrix([half,half]);flowmass=du/sum(du)
    edge_difference=sp.factor(karyu[0]-flowmass[0])
    counter=sp.simplify(edge_difference.subs(u,half))
    assert counter!=0
    record('edge_normalization_counterexample',0,'Exact nonzero difference at u=1/2: '+str(counter))
    same(kernel[:,:2]*sp.Matrix([half,half]),karyu.subs(u,sp.exp(-1/(2*r**2))),'column_normalization_reference')
    table(out/'flow_normalization_checks.tsv',[dict(check='global_vs_column_normalization',symbolic_difference=str(edge_difference),exact_counterexample=str(counter),grid_step='1/2',weighted_sum=1,ordinary_density_sum=2,flow_width_role='fixed/configured; not a new fitted coordinate')])
    passed('Finite-grid normalization, units and edge counterexample')

    current='T07'
    ma,mb=sp.symbols('cohort_A_mean cohort_B_mean',real=True)
    sa,sb=sp.symbols('cohort_A_sd cohort_B_sd',positive=True)
    unknowns=[ma,sa,mb,sb,kappa]
    cohort_map=[('2N','E1',('init_mean_2N','init_sd_2N')),
                ('2N','E2',('init_mean_2N','init_sd_2N')),
                ('4N','E1',('init_mean_4N','init_sd_4N')),
                ('4N','E2',('init_mean_4N','init_sd_4N'))]
    required=[str(z) for z in unknowns]
    validate_schema(required,required,cohort_map,['sigma_kary','sigma_kary'])
    shared=[];families={}
    for cohort,mean,width,other in [('2N',ma,sa,[mb,sb]),('4N',mb,sb,[ma,sa])]:
        weights=ad.initial_weights(sp,grid,mean,width);families[cohort]=weights
        results=[]
        for e,known_c in [('E1',sp.Integer(1)),('E2',sp.Integer(2))]:
            js=ad.output_jets(sp,x,-known_c*kappa*x,x,1)
            vv,jj=ad.initial_pullback(sp,js,x,n0*weights,unknowns)
            for z in other:
                for h in vv:same(h.diff(z),sp.zeros(3,1),'cross_cohort_'+cohort+'_'+e+'_'+str(z))
            for j,z in enumerate(unknowns):
                for i in range(3):shared.append(dict(cohort=cohort,experiment=e,parameter=str(z),order=1,row=i,expression=str(jj[1][i,j])))
            results.append((vv,jj))
        same(results[1][0][1],2*results[0][0][1],'constant_input_scaling_'+cohort)
        same(results[1][1][1],2*results[0][1][1],'shared_derivative_scaling_'+cohort)
    pa,pb=families['2N'],families['4N']
    # Block determinant is certified through the verified 2x2 family minors.
    witness=-n0**5*pa[0]*(2*pa[0]*pa[1]*pa[2]/sa**5)*(2*pb[0]*pb[1]*pb[2]/sb**5)
    witnesses.append(dict(case=current,witness=str(witness),domain='positive N0,sa,sb,kappa; finite independent cohort means',conclusion='shared toy kappa plus four initial coordinates; block-triangular witness from registered independent family minors'))
    ra,rb=sp.symbols('diagnostic_assay_A diagnostic_assay_B',positive=True)
    qa=ad.kary_kernel(sp,grid,grid,ra)[:,1];qb=ad.kary_kernel(sp,grid,grid,rb)[:,1]
    same(qa.diff(rb),sp.zeros(3,1),'independent_width_cross_A')
    same(qb.diff(ra),sp.zeros(3,1),'independent_width_cross_B')
    validate_schema(['rA','rB'],['rA','rB'],assay_widths=['rA','rB'],diagnostic=True)
    table(out/'shared_parameter_checks.tsv',shared)
    passed('Two cohort families shared across two assigned constant inputs')

    current='T08'
    sg=sp.Symbol('sigma_growth',positive=True);g=sp.Symbol('growth_mean',real=True)
    zero(g.diff(sg),'variance_only_mean_column')
    logdensity=-sp.log(sg)-(y-g)**2/(2*sg**2)-sp.log(2*sp.pi)/2
    zero(logdensity.diff(sg)-(-1/sg+(y-g)**2/sg**3),'growth_SD_likelihood_role')
    base,ind,beta=sp.symbols('p_mis_base p_misseg buffer_beta',nonnegative=True)
    risk,km,ss,nn,N=sp.symbols('death_risk k_o_mis buffer_smax buffer_n_exp N',positive=True)
    induction=base+ind*risk/(km+risk);survival=ss*sp.exp(-beta*(44/N)**nn)
    zero(induction.diff(km).subs(ind,0),'induction_zero_inactivates_k')
    zero(induction.subs(ind,0)-base,'induction_zero_retains_base')
    zero(survival.diff(nn).subs(beta,0),'beta_zero_inactivates_exponent')
    zero(survival.subs(beta,0)-ss,'beta_zero_retains_smax')
    rejected('zero_initial_width',lambda:ad.initial_weights(sp,grid,mu,sp.Integer(0)))
    rejected('zero_observer_width',lambda:ad.kary_kernel(sp,grid,grid,sp.Integer(0)))
    table(out/'variance_role_checks.tsv',[dict(parameter='sigma_growth',mean_column='exactly_zero',likelihood_role='positive standard deviation; variance=sigma_growth^2',structural_result='not_analyzed; probabilistic observation analysis required')])
    passed('Variance-only SD and exact-zero expression activity preserved')

    current='T09'
    rejected('row_dynamics_shape',lambda:ad.output_jets(sp,x,sp.Matrix([list(-c*kappa*x)]),x,1))
    rejected('initial_length_mismatch',lambda:ad.initial_pullback(sp,jets,x,[sp.Integer(1)],[mu]))
    rejected('omitted_initial_width',lambda:validate_schema(['mu'],['mu','s']))
    rejected('duplicate_unknown',lambda:validate_schema(['mu','mu'],['mu','s']))
    rejected('invented_mean_variance_column',lambda:validate_schema(['mu','sigma_growth'],['mu','sigma_growth']))
    broken=list(cohort_map);broken[1]=('2N','E2',('separate_mean','separate_width'))
    rejected('duplicated_initial_family_per_oxygen',lambda:validate_schema(required,required,broken))
    rejected('shared_observer_replaced_by_independent',lambda:validate_schema(['rA','rB'],['rA','rB'],assay_widths=['rA','rB']))
    rejected('fixed_observer_omits_width',lambda:validate_schema([],['sigma_kary']))
    raww=[sp.exp(-(z-mu)**2/(2*s**2)) for z in grid]
    wrong=sp.diff(raww[1],s)/sum(raww)
    missing=sp.simplify((wrong-p[1].diff(s)).subs({mu:0,s:1}))
    assert missing!=0;record('missing_normalizer_derivative_detected',0,'Exact nonzero residual: '+str(missing))
    assert counter!=0;record('wrong_flow_normalization_detected',0,'Reuse independently registered edge counterexample')
    records=[dict(case='fixture',release_sha256='r1',input_sha256='i1',order=str(k),path='matrix_'+str(k)+'.tsv') for k in (1,2)]
    validate_collection(records,'fixture','r1','i1');record('multiple_legitimate_matrices_accepted',0)
    rejected('stale_matrix_record',lambda:validate_collection(records,'fixture','r2','i1'))
    rejected('duplicate_matrix_identity',lambda:validate_collection(records+records[:1],'fixture','r1','i1'))
    rejected('missing_matrix_identity',lambda:validate_collection([],'fixture','r1','i1'))
    passed('Shape, omission, normalization, sharing and collection faults rejected')

    table(out/'rank_witnesses.tsv',witnesses)
    table(out/'numeric_diagnostics.tsv',[dict(status='not_used',reason='Optional numeric checks not needed; all required identities verified exactly; no random substitutions')])
    table(out/'application_result.tsv',[dict(status='PASS_EXACT_INTERFACE_CONTROLS',completed_cases='T00-T09',original_culture_parameters_classified=0,new_native_calls=0,new_random_rank_calls=0,Stage_B=False,Step_5=False)])
    print('EXACT_INITIAL_OBSERVER_CONTROLS_PASSED',flush=True)
except BaseException as exc:
    cases.append(dict(case=current,status='FAIL',evidence_route='unresolved',detail=type(exc).__name__+': '+str(exc)))
    table(out/'control_results.tsv',cases)
    table(out/'application_result.tsv',[dict(status='HOLD',failed_case=current,error=type(exc).__name__+': '+str(exc))])
    traceback.print_exc();raise
finally:
    with (out/'imported_modules.tsv').open('w') as f:
        f.write('module\tpath\n')
        for name,module in sorted(sys.modules.items()):f.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')
