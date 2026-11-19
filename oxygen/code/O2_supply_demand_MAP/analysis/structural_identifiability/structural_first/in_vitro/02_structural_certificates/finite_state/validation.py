"""Finite arithmetic, fixed-model derivative and structural symmetry controls."""
from pathlib import Path
import sys
from runtime import boot
from protocol import table,MECHANISMS,INITIALS,check_payload,successful_worker

study=Path(sys.argv[1]);task=sys.argv[2]
sp,flint,np,out=boot(study,task)
from algebra import rates,exact,initial,taylor,encode,decode,nonzero,independent_minor_check
arb,arb_mat,fmpq,fmpq_mat=flint.arb,flint.arb_mat,flint.fmpq,flint.fmpq_mat
checks=[]
def check(name,condition,detail=''):
    checks.append(dict(check=name,passed=bool(condition),detail=detail));table(out/'validation_checks.tsv',checks)
    assert condition,name+': '+detail
def zero(name,expr):
    v=sp.simplify(sp.expand_log(expr,force=False))
    if v!=0:v=sp.simplify(sp.cancel(sp.together(v)))
    check(name,v==0,str(v))

flint.ctx.prec=128
M=arb_mat([[1,2],[3,4]])
check('exact_nonsingular_matrix',M.det()==-2)
ok,bound=independent_minor_check(M);check('independent_neumann_certificate',ok,str(bound))
check('exact_singular_matrix',arb_mat([[1,2],[2,4]]).det().is_zero())
delta=fmpq(1,2**180)
near=arb_mat([[1,1],[1,arb(1+delta)]])
check('low_precision_is_inconclusive',near.det().contains(0))
flint.ctx.prec=256
near=arb_mat([[1,1],[1,arb(1+delta)]])
check('higher_precision_separates_zero',nonzero(near.det()) and near.det().contains(delta))
ok,bound=independent_minor_check(near);check('near_singular_independent_verification',ok,str(bound))
check('invalid_balls_rejected',not nonzero(arb('nan')) and not nonzero(arb('inf')))
x=exact('1/3').exp();check('ball_roundtrip_outward',decode(encode(x)).contains(x))
check('elementary_interval_identity',((x.log()-exact('1/3')).contains(0)))
check('rational_input_is_not_binary_float',exact('1/10').contains(fmpq(1,10)))
R,rank=fmpq_mat([[1,2,0],[0,0,1]]).rref()
check('midpoint_rref_nomination_interface',rank==2 and R[0,0]==1 and R[1,2]==1)
print('ARITHMETIC_CONTROLS_PASSED',flush=True)

# The 14 original variables remain symbolic, including unknown exponents.
p={n:sp.Symbol(n,positive=True) for n in MECHANISMS}
u,O=sp.symbols('ratio oxygen',positive=True)
values,grads=rates(p,u,O,log=sp.log)
for label,value,gradient in zip(('lambda','mu','p','logS'),values,grads):
    for name in MECHANISMS:zero('response_'+label+'_'+name,sp.diff(value,p[name])-gradient[name])
print('ORIGINAL_RESPONSE_DERIVATIVES_PASSED',flush=True)

# Direct state-space Lie polynomials independently validate the fixed recurrence.
a,eta=sp.symbols('toy_rate toy_initial_fraction',positive=True)
x1,x2=sp.symbols('toy_x1 toy_x2')
xv=sp.Matrix([x1,x2]);G=a*sp.Matrix([[-1,sp.Rational(1,3)],[2,-sp.Rational(1,3)]])
mu=[sp.Rational(1,5),sp.Rational(1,7)];K=sp.Integer(10)
g0=sp.Matrix([eta,1-eta]);D0=g0.jacobian([a,eta])
A,D=taylor(G,[G.diff(a)],mu,[[0,0]],g0,D0,4,K,known_initial_total=1,zero=sp.zeros,copy=sp.Matrix)
f=(1-(x1+x2)/K)*(G*xv)-sp.diag(*mu)*xv
H=xv
for k in range(5):
    direct=H.subs({x1:eta,x2:1-eta},simultaneous=True)/sp.factorial(k)
    jac=direct.jacobian([a,eta])
    for i in range(2):
        check('toy_Taylor_value_'+str(k)+'_'+str(i),sp.Poly(sp.expand(A[k][i]-direct[i]),a,eta).is_zero)
        for j in range(2):check('toy_Taylor_derivative_'+str(k)+'_'+str(i)+'_'+str(j),sp.Poly(sp.expand(D[k][i,j]-jac[i,j]),a,eta).is_zero)
    if k<4:H=(H.jacobian(xv)*f).applyfunc(sp.expand)
at_capacity=f.subs({x1:K/2,x2:K/2})
zero('direct_death_not_crowding_scaled_1',at_capacity[0]+mu[0]*K/2)
zero('direct_death_not_crowding_scaled_2',at_capacity[1]+mu[1]*K/2)
print('FIXED_MODEL_TAYLOR_RECURRENCE_PASSED',flush=True)

# Gaussian factorization is checked by a general common-factor identity.
w0,w1,w2,c=sp.symbols('w0 w1 w2 common',positive=True)
zero('Gaussian_common_factor_cancellation',c*w0/(c*w0+c*w1+c*w2)-w0/(w0+w1+w2))
pp={'init_mean_2N':exact('1/3'),'init_sd_2N':exact('3/2')}
aa,dd=initial(pp,'2N',list(INITIALS[:2]),grid=[-1,0,1],N0=1)
check('Gaussian_ball_normalization',(sum(aa[i,0] for i in range(3))-1).contains(0))
for j in range(2):check('Gaussian_derivative_sum_'+str(j),sum(dd[i,j] for i in range(3)).contains(0))
minor=arb_mat([[dd[i,j] for j in range(2)] for i in range(2)])
expected=2*aa[0,0]*aa[1,0]*aa[2,0]/pp['init_sd_2N']**5
check('factored_Gaussian_minor',nonzero(minor.det()) and (minor.det()-expected).contains(0))
ok,bound=independent_minor_check(minor);check('Gaussian_independent_minor_verifier',ok,str(bound))

# Explicit single-oxygen family is an original-parameter calculation.
Cnew,nnew=sp.symbols('new_O2_crit new_n_O',positive=True)
h=p['O2_crit']**p['n_O']/(p['O2_crit']**p['n_O']+O**p['n_O'])
hnew=Cnew**nnew/(Cnew**nnew+O**nnew)
changed=dict(p,O2_crit=Cnew,n_O=nnew,alpha_o2=p['alpha_o2']*h/hnew,mu_hp=p['mu_hp']*h/hnew)
newvalues,_=rates(changed,u,O,log=sp.log)
for label,v,nv in zip(('lambda','mu','p','logS'),values,newvalues):zero('single_O_invariance_'+label,nv-v)
table(out/'single_input_invariance.tsv',[dict(component=k,residual=0,domain='positive O,C,n and interior transformed alpha/mu; original initial/assay variables unchanged') for k in ('lambda','mu','p','logS','mask_as_function_of_p','G_and_retained_rates_by_event_identity')])
zero('induction_null_inactive_k',sp.diff(values[2],p['k_o_mis']).subs(p['p_misseg'],0))
zero('survival_null_inactive_exponent',sp.diff(values[3],p['buffer_n_exp']).subs(p['buffer_beta'],0))
zero('zero_oxygen_Hill_limit',h.subs(O,0)-1)
S,B,W=sp.symbols('retained_daughters division_buffer WGD_probability')
birth_sum=(1-W)*S+W-1;buffer=(1-W)*(2-S)
zero('complete_event_account_identity',birth_sum+buffer-(1-W))

# Existing real-series quotient arithmetic must preserve normalization gradients.
flint.ctx.cap=5
num=flint.arb_series([1,2,3],prec=3);den=flint.arb_series([2,1,1],prec=3)
q=num/den;product=q*den-num
for i in range(3):check('flow_series_quotient_'+str(i),product[i].contains(0))
derivative=(flint.arb_series([1,0,2],prec=3)*den-num*flint.arb_series([0,1,1],prec=3))/(den*den)
check('flow_series_derivative_exists',all(derivative[i].is_finite() for i in range(3)))

for name,rows,task_,rel in [('missing',[],'x','r'),('stale',[dict(task='x',release_sha256='old',item='m')],'x','r'),
 ('duplicate',[dict(task='x',release_sha256='r',item='m')]*2,'x','r')]:
    try:check_payload(rows,task_,rel)
    except ValueError:check('collector_reject_'+name,True)
    else:raise AssertionError(name)
complete=[dict(status='PASS')]
check('complete_worker_record_accepted',successful_worker(complete,complete,complete,'COMPLETED','0:0'))
check('missing_OOM_finally_is_not_success',not successful_worker(complete,[],complete,'COMPLETED','0:0'))
check('missing_policy_is_not_success',not successful_worker(complete,complete,[],'COMPLETED','0:0'))
check('OOM_scheduler_state_is_not_success',not successful_worker(complete,complete,complete,'OUT_OF_MEMORY','0:9'))
table(out/'validation_gate.tsv',[dict(status='PASS',exact_checks=len(checks),native_calls=0,original_variables=14,original_nuisances=6,new_general_solver=False)])
table(out/'application_result.tsv',[dict(status='PASS',task=task,phase='B1',project_structural_cases=0)])
print('STAGE_B_ARITHMETIC_AND_DERIVATIVE_GATE_PASSED',len(checks),flush=True)
