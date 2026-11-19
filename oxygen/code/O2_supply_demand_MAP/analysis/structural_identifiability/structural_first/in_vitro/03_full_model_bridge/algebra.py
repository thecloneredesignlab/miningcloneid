"""Fixed culture response derivatives, Taylor recurrence and ball certificates."""
from fractions import Fraction
from flint import arb,arb_mat,fmpq,fmpq_mat
from protocol import MECHANISMS,INITIALS,read

def exact(value):
    q=Fraction(str(value));return arb(fmpq(q.numerator,q.denominator))

def point(study,witness):
    p={}
    for row in read(study/'inputs/witnesses.tsv'):
        if row['witness']!=witness:continue
        lo,hi,q=exact(row['lower']),exact(row['upper']),exact(row['fraction'])
        rule=row['rule']
        value=exact(row['value']) if rule=='explicit' else lo*((hi/lo)**q) if rule=='geometric' else ((1+lo)*(((1+hi)/(1+lo))**q)-1) if rule=='log1p' else lo+(hi-lo)*q
        assert value.is_finite() and value>lo and value<hi,row['parameter']
        p[row['parameter']]=value
    assert len(p)==20
    return p

def rates(p,u,O,log=lambda v:v.log()):
    """Production ploidy_related, with independently checked manual gradients."""
    z=p['lam_max']*0
    H=p['O2_crit']**p['n_O']/(p['O2_crit']**p['n_O']+O**p['n_O'])
    HC=p['n_O']*H*(1-H)/p['O2_crit']
    Hn=H*(1-H)*log(p['O2_crit']/O)
    ug=u**p['gamma_growth'];um=u**p['gamma_mu']
    den=1+p['alpha_o2']*H*ug
    lam=p['lam_max']/den;mu=p['mu_hp']*H*um
    prob=p['p_mis_base']+p['p_misseg']*mu/(p['k_o_mis']+mu)
    v=(1/u)**p['buffer_n_exp']
    logs=log(p['buffer_smax'])-p['buffer_beta']*v
    dl={n:z for n in MECHANISMS};dm=dict(dl);dp=dict(dl);ds=dict(dl)
    dl['lam_max']=1/den
    dl['alpha_o2']=-p['lam_max']*H*ug/den**2
    dl['gamma_growth']=-p['lam_max']*p['alpha_o2']*H*ug*log(u)/den**2
    dl['O2_crit']=-p['lam_max']*p['alpha_o2']*ug*HC/den**2
    dl['n_O']=-p['lam_max']*p['alpha_o2']*ug*Hn/den**2
    dm['mu_hp']=H*um;dm['gamma_mu']=mu*log(u)
    dm['O2_crit']=p['mu_hp']*um*HC;dm['n_O']=p['mu_hp']*um*Hn
    for n in MECHANISMS:dp[n]=p['p_misseg']*p['k_o_mis']*dm[n]/(p['k_o_mis']+mu)**2
    dp['p_mis_base']+=1;dp['p_misseg']+=mu/(p['k_o_mis']+mu)
    dp['k_o_mis']-=p['p_misseg']*mu/(p['k_o_mis']+mu)**2
    ds['buffer_smax']=1/p['buffer_smax'];ds['buffer_beta']=-v
    ds['buffer_n_exp']=-p['buffer_beta']*v*log(1/u)
    return (lam,mu,prob,logs),[dl,dm,dp,ds]

def initial(p,cohort,columns,grid=range(22,155),N0=1000000):
    grid=list(grid);mean=p['init_mean_'+cohort];s=p['init_sd_'+cohort]
    # The common exp(-mean^2/(2*s^2)) factor cancels exactly.
    w=[(arb(n)*mean/s**2-arb(n*n)/(2*s**2)).exp() for n in grid]
    total=sum(w);probs=[v/total for v in w]
    ez=sum(arb(n)*v for n,v in zip(grid,probs))
    e2=sum((arb(n)-mean)**2*v for n,v in zip(grid,probs))
    a=arb_mat([[arb(N0)*v] for v in probs]);D=arb_mat(len(grid),len(columns))
    for i,(n,v) in enumerate(zip(grid,probs)):
        if 'init_mean_'+cohort in columns:D[i,columns.index('init_mean_'+cohort)]=N0*v*(n-ez)/s**2
        if 'init_sd_'+cohort in columns:D[i,columns.index('init_sd_'+cohort)]=N0*v*((n-mean)**2-e2)/s**3
    return a,D

def dim(m):return (m.nrows(),m.ncols()) if hasattr(m,'nrows') else (m.rows,m.cols)

def taylor(G,dG,mu,dmu,a0,D0,order,K,known_initial_total=None,zero=arb_mat,copy=arb_mat):
    """Only the registered quadratic culture recurrence; not a generic solver."""
    n,d=dim(D0);theta=len(dG)
    A=[a0];D=[D0];GA=[];GD=[];Q=[];Ls=[];dLs=[]
    for k in range(order):
        ga=G*A[k];gd=G*D[k];q=zero(n,d)
        for j in range(theta):
            v=dG[j]*A[k]
            for i in range(n):q[i,j]=v[i,0]
        GA.append(ga);GD.append(gd);Q.append(q)
        L=known_initial_total if k==0 and known_initial_total is not None else sum(A[k][i,0] for i in range(n))
        dL=[0 for j in range(d)] if k==0 and known_initial_total is not None else [sum(D[k][i,j] for i in range(n)) for j in range(d)]
        Ls.append(L);dLs.append(dL)
        av=copy(ga);dv=copy(gd+q)
        for i in range(n):
            av[i,0]-=mu[i]*A[k][i,0]
            for j in range(d):
                dv[i,j]-=mu[i]*D[k][i,j]
                if j<theta:dv[i,j]-=dmu[j][i]*A[k][i,0]
        for ell in range(k+1):
            b=k-ell
            for i in range(n):
                av[i,0]-=Ls[ell]*GA[b][i,0]/K
                for j in range(d):dv[i,j]-=(dLs[ell][j]*GA[b][i,0]+Ls[ell]*(GD[b][i,j]+Q[b][i,j]))/K
        A.append(av/(k+1));D.append(dv/(k+1))
    return A,D

def encode(value):
    assert value.is_finite()
    m,e=value.mid().man_exp();r,re=value.rad().man_exp()
    return dict(mid_man=str(m),mid_exp=str(e),rad_man=str(r),rad_exp=str(re))
def decode(row):
    return arb((int(row['mid_man']),int(row['mid_exp'])),(int(row['rad_man']),int(row['rad_exp'])))
def nonzero(value):return value.is_finite() and not value.contains(0)

def independent_minor_check(M):
    """Exact rational midpoint preconditioner plus a rigorous Neumann bound."""
    n,m=dim(M);assert n==m
    mid=fmpq_mat([[M[i,j].mid().fmpq() for j in range(n)] for i in range(n)])
    try:R=mid.inv()
    except ZeroDivisionError:return False,arb('inf')
    E=arb_mat([[int(i==j) for j in range(n)] for i in range(n)])-arb_mat(R)*M
    bounds=[sum(abs(E[i,j]) for j in range(n)) for i in range(n)]
    finite=all(v.is_finite() for v in bounds)
    valid=finite and all(v<1 for v in bounds)
    upper=max(v.upper() for v in bounds) if finite else arb('inf')
    return valid,upper
