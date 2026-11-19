"""Exact-support, fixed-mask 133-state culture operator and declared outputs."""
import hashlib
from flint import arb,arb_mat,arb_series
from algebra import rates,exact,initial,taylor
from protocol import MECHANISMS,INITIALS

GRID=tuple(range(22,155));N=len(GRID)

class UnresolvedBranch(ValueError):pass

def operator(p,O):
    G=arb_mat(N,N);dG=[arb_mat(N,N) for j in MECHANISMS]
    mus=[];dmus=[[] for j in MECHANISMS];values=[];ledger=[];masks=[]
    eps=exact('1e-8');w=p['p_wgd']
    for j,n in enumerate(GRID):
        (lam,mu,prob,logs),grads=rates(p,arb(n)/44,O)
        assert prob>0 and prob<1 and lam>0 and mu>0
        dl,dm,dp,ds=grads;mus.append(mu)
        for h,name in enumerate(MECHANISMS):dmus[h].append(dm[name])
        ordinary=lam*(1-w);whole=lam*w
        dordinary=[dl[name]*(1-w)-lam*int(name=='p_wgd') for name in MECHANISMS]
        dwhole=[dl[name]*w+lam*int(name=='p_wgd') for name in MECHANISMS]
        pi=(1-prob)**n;kept=[];S_in=arb(0);bio=arb(0);prune=arb(0);low=arb(0);high=arb(0)
        margin=None
        for m in range(n+1):
            if m:pi=pi*(n-m+1)/m*prob/(1-prob)
            if pi>eps:keep=True
            elif pi<eps:keep=False
            else:raise UnresolvedBranch('Probability threshold unresolved: N='+str(n)+';m='+str(m))
            kept.append('1' if keep else '0')
            gap=abs(pi-eps).abs_lower()
            if margin is None or gap<margin:margin=gap
            survival=(m*logs).exp() if m else arb(1)
            viable=pi*survival
            bio+=2*pi*(1-survival)
            if not keep:prune+=2*viable;continue
            factor=m/prob-(n-m)/(1-prob)
            dviable=[viable*(factor*dp[name]+m*ds[name]) for name in MECHANISMS]
            for destination in (n+m,n-m):
                if destination<22:low+=viable
                elif destination>154:high+=viable
                else:
                    i=destination-22;S_in+=viable
                    G[i,j]+=ordinary*viable
                    for h in range(14):dG[h][i,j]+=dordinary[h]*viable+ordinary*dviable[h]
        G[j,j]-=lam
        for h,name in enumerate(MECHANISMS):dG[h][j,j]-=dl[name]
        if 2*n<=154:
            G[2*n-22,j]+=whole
            for h in range(14):dG[h][2*n-22,j]+=dwhole[h]
        boundary_wgd=whole if 2*n>154 else arb(0)
        b=ordinary*(2-S_in)+boundary_wgd
        bio*=ordinary;prune*=ordinary;low*=ordinary;high=ordinary*high+boundary_wgd
        account=sum(G[i,j] for i in range(N))+b-lam*(1-w)
        assert account.contains(0) and (b-bio-prune-low-high).contains(0)
        ledger.append(dict(N=n,division=lam,direct_death=mu,buffer=b,biological=bio,probability_pruning=prune,lower_drop=low,upper_drop=high,account_residual=account))
        values.append(dict(N=n,lam=lam,mu=mu,probability=prob,survival=logs.exp(),buffer=b))
        masks.append(dict(N=n,mask=''.join(kept),minimum_margin=margin))
    signature=hashlib.sha256('\n'.join(str(r['N'])+':'+r['mask'] for r in masks).encode()).hexdigest()
    return dict(G=G,dG=dG,mu=mus,dmu=dmus,values=values,ledger=ledger,masks=masks,mask_sha256=signature)

def kary_kernel(width):
    K=arb_mat(N,N);dK=arb_mat(N,N)
    for j,n in enumerate(GRID):
        raw=[(-(arb(i-n)**2)/(2*width**2)).exp() for i in GRID]
        total=sum(raw);probs=[v/total for v in raw]
        moment=sum(v*(i-n)**2 for i,v in zip(GRID,probs))
        for i,(m,v) in enumerate(zip(GRID,probs)):
            K[i,j]=v;dK[i,j]=v*((m-n)**2-moment)/width**3
    return K,dK

def flow_kernel(flow_rows,step):
    grid=[exact(r['ploidy']) for r in flow_rows]
    width=2*step
    if width<exact('1/22'):width=exact('1/22')
    return arb_mat([[(-(v-arb(n)/22)**2/(2*width**2)).exp() for n in GRID] for v in grid])

def columns_for(case):
    result=list(MECHANISMS)
    if case['unknown_initial']=='True':result+=list(INITIALS)
    if case['unknown_observer']=='True':result+=['sigma_kary']
    return result

def output_rows(p,case,cohort,O,op,flow_rows=None,flow_step=None):
    columns=columns_for(case);order=int(case['max_order']);d=len(columns)
    a0,D0=initial(p,cohort,columns)
    A,D=taylor(op['G'],op['dG'],op['mu'],op['dmu'],a0,D0,order,arb(10**12),known_initial_total=arb(10**6))
    output=case['output'];labels=[];rows=[]
    if output=='full':
        for k in range(order+1):
            for i,n in enumerate(GRID):
                labels.append(dict(cohort=cohort,oxygen=str(O),output='live_state',N=n,order=k))
                rows.append([D[k][i,j] for j in range(d)])
    else:
        for k in range(order+1):
            labels.append(dict(cohort=cohort,oxygen=str(O),output='live_total',N='',order=k))
            rows.append([arb(0) for j in range(d)] if k==0 else [sum(D[k][i,j] for i in range(N)) for j in range(d)])
        if output in ('kary','kary_flow'):
            K,dK=kary_kernel(p['sigma_kary']);observer=columns.index('sigma_kary')
            for k in range(order+1):
                J=K*D[k];extra=dK*A[k]
                for i in range(N):J[i,observer]+=extra[i,0]
                for i,n in enumerate(GRID[:-1]):
                    labels.append(dict(cohort=cohort,oxygen=str(O),output='L_times_kary_probability',N=n,order=k))
                    rows.append([J[i,j] for j in range(d)])
        if output=='kary_flow':
            assert flow_rows and flow_step is not None
            F=flow_kernel(flow_rows,flow_step);C=[F*a for a in A];CJ=[F*da for da in D]
            nf=len(flow_rows);Z=[flow_step*sum(v[i,0] for i in range(nf)) for v in C]
            DZ=[[flow_step*sum(v[i,j] for i in range(nf)) for j in range(d)] for v in CJ]
            denominator=arb_series(Z,prec=order+1);assert Z[0]>0
            for i in range(nf-1):
                numerator=arb_series([v[i,0] for v in C],prec=order+1)
                derivatives=[]
                for j in range(d):
                    dn=arb_series([v[i,j] for v in CJ],prec=order+1)
                    dz=arb_series([v[j] for v in DZ],prec=order+1)
                    derivatives.append((dn*denominator-numerator*dz)/(denominator*denominator))
                for k in range(order+1):
                    labels.append(dict(cohort=cohort,oxygen=str(O),output='flow_density',N=flow_rows[i]['ploidy'],order=k))
                    rows.append([v[k] for v in derivatives])
    return columns,labels,rows
