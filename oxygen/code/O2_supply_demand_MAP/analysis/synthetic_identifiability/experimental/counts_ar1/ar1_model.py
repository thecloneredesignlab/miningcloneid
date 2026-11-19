"""Joint Gaussian AR(1) likelihood for the 221 audited post-inoculation counts."""
import ctypes as C
from hashlib import sha256
import os
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize, nnls

HERE = Path(__file__).resolve().parent
COUNTS = HERE.with_name("20261118_experimental_counts_practical_01")
sys.path.insert(0,str(COUNTS))
os.environ["CULTURE_COUNT_LIB"]=os.environ.get("CULTURE_AR1_LIB",
    "/private/tmp/culture_count_ar1_practical_01.dylib")
from count_model import CountKernel
from practical import Configuration, FD, FI, MECHANISMS, read, row_x, write
from inference import THRESHOLD, likelihood_ratio, pool

TARGETS=("lam_max","buffer_smax","buffer_beta","buffer_n_exp","n_O",
         "p_misseg","k_o_mis","gamma_growth","mu_hp","O2_crit")
RHO_MAX=.99
STATIONARY=.025
MIN_BUFFER=1e-9


class CalendarAR1(CountKernel):
    def __init__(self):
        super().__init__()
        self.lib.ar1_count_loss.argtypes=[C.c_void_p,FD,C.c_int,FI,FI,FI,FD,FD,C.c_void_p]
        self.lib.ar1_count_loss.restype=C.c_int
        g=self.group[self.capture]
        self.groups=g
        self.days=np.array([int(r["calendar_day"]) for r in self.observations if r["fitted_count"]])
        self.previous=np.full(self.n,-1,dtype=np.int32)
        self.transition=np.flatnonzero(g[1:]==g[:-1])+1
        self.previous[self.transition]=self.transition-1
        self.gaps=np.zeros(self.n,dtype=np.int32)
        self.gaps[self.transition]=self.days[self.transition]-self.days[self.transition-1]
        assert self.n==221 and len(np.unique(g))==114 and len(self.transition)==107
        assert np.array_equal(np.bincount(g).nonzero()[0],np.unique(g))
        assert all(np.count_nonzero(self.previous[g==z]<0)==1 for z in np.unique(g))
        assert np.all(self.gaps[self.transition]>0)

    def metrics(self,p,target,gradient=False):
        p=np.ascontiguousarray(p)
        target=np.ascontiguousarray(target)
        out=np.empty(4)
        grad=np.empty(47) if gradient else None
        ok=self.lib.ar1_count_loss(self.handle,p,self.n,self.indices,self.previous,self.gaps,
                                   target,out,None if grad is None else grad.ctypes.data)
        if not ok:
            return np.array([1e13,-10.,np.nan,np.nan]), np.zeros(47) if gradient else None
        return out,grad

    def constrained(self,p,target,conf,gradient=True):
        out,g=self.metrics(p,target,gradient)
        return out[0],g,out[1]

    def supply(self,p,gaussian=False,pick=None):
        margin,g=super().supply(np.ascontiguousarray(p[:46]),gaussian,pick)
        return margin,np.pad(g,((0,0),(0,1)))

    def predict_logs(self,p):
        return super().log_predict(np.ascontiguousarray(p[:46]))[0][self.capture]

    def noise_covariance(self,rho):
        gap=abs(self.days[:,None]-self.days[None,:])
        return np.where(self.groups[:,None]==self.groups[None,:],rho**gap,0.)


class Config:
    def __init__(self,kernel):
        self.base=Configuration(kernel)
        self.names=[*self.base.names,"rho"]
        self.active=np.r_[self.base.active,46]
        self.x_anchor=np.r_[self.base.x_anchor,.806/RHO_MAX]

    def unmap(self,x):
        p,jac=self.base.unmap(x[:-1])
        return np.r_[p,x[-1]*RHO_MAX],np.r_[jac,RHO_MAX]

    def gradient(self,g,transform_jacobian):
        return np.r_[self.base.gradient(g[:-1],transform_jacobian[:-1]),g[-1]*RHO_MAX]

    def load_old(self,r,rho=.806):
        return np.r_[row_x(self.base,r),rho/RHO_MAX]


MODEL=None


def model():
    global MODEL
    if MODEL is None:
        MODEL=CalendarAR1()
    return MODEL


def solve(k,c,y,start,limit=90,fixed=None,divisor=12.,margin_buffer=MIN_BUFFER):
    """Scale the objective, generate donor constraints and validate every live supply."""
    x=np.clip(np.asarray(start),0,1).copy()
    fixed={} if fixed is None else dict(fixed)
    for i,v in fixed.items():
        x[i]=v
    active=np.array([i for i in range(len(x)) if i not in fixed])
    search=x[active].copy()
    y=np.ascontiguousarray(y)
    begin=time.monotonic()
    calls,iterations=0,0
    best=None

    def full(u):
        z=x.copy()
        z[active]=u
        p,j=c.unmap(z)
        return p,j,z

    def objective(u):
        nonlocal calls,best
        p,j,z=full(u)
        calls+=1
        f,g,s=k.constrained(p,y,c)
        if s>=0 and (best is None or f<best[0]):
            best=(f,z.copy(),s)
        return f/divisor,c.gradient(g,j)[active]/divisor

    m,_=k.supply(c.unmap(x)[0])
    smallest=np.argsort(m)[:4]
    selected=set(smallest[m[smallest]<.06])
    for call in range(1,7):
        take=np.array(sorted(selected),dtype=np.int32)

        def gfun(u):
            m,_=k.supply(full(u)[0])
            return m[take]-margin_buffer

        def gjac(u):
            p,j,_=full(u)
            _,g=k.supply(p,pick=take)
            return np.array([c.gradient(v,j)[active] for v in g])

        res=minimize(objective,search,jac=True,method="SLSQP",bounds=[(0.,1.)]*len(active),
                     constraints=[dict(type="ineq",fun=gfun,jac=gjac)] if len(take) else (),
                     options=dict(maxiter=limit,ftol=2e-9))
        search=np.clip(res.x,0.,1.)
        iterations+=res.nit
        p,j,z=full(search)
        m,_=k.supply(p)
        new=set(np.flatnonzero(m<0))-selected
        if not new:
            break
        selected.update(sorted(new,key=lambda i:m[i])[:4])
    f,g,s=k.constrained(p,y,c)
    if best is not None and (s<0 or best[0]<f-1e-6):
        z=best[1]
        p,j=c.unmap(z)
        f,g,s=k.constrained(p,y,c)
        m,_=k.supply(p)
    binding=np.flatnonzero(m<=5e-7)
    _,jg=k.supply(p,pick=binding)
    normals=[c.gradient(v,j)[active] for v in jg]
    for i,pos in enumerate(active):
        if z[pos]<5e-7 or z[pos]>1.-5e-7:
            v=np.zeros(len(active))
            v[i]=1. if z[pos]<.5 else -1.
            normals.append(v)
    g=c.gradient(g,j)[active]
    norm=nnls(np.array(normals).T,g,maxiter=1000)[1] if normals else np.linalg.norm(g)
    return dict(x=z,loss=f,supply=s,gradient=float(norm),iterations=iterations,calls=calls,
                success=bool(res.success),termination=str(res.message),rounds=call,
                seconds=time.monotonic()-begin)


def receipt(k,c,y,f,**labels):
    p,_=c.unmap(f["x"])
    out,_=k.metrics(p,y)
    return dict(labels,loss=float(out[0]),donor_margin=float(out[1]),stationarity=f["gradient"],
                rho=float(p[-1]),sigma=float(np.sqrt(out[2]/k.n)),
                quadratic=float(out[2]),logdet=float(out[3]),calls=f["calls"],iterations=f["iterations"],
                success=f["success"],termination=f["termination"],rounds=f["rounds"],seconds=f["seconds"],
                **dict(zip((f"x_{name}" for name in c.names),f["x"])),
                **{name:float(p[j]) for name,j in zip(c.names[:-1],c.active[:-1])})


def choose(rows,n=221):
    allowed=[r for r in rows if float(r["donor_margin"])>=0 and float(r["loss"])<1e12]
    if not allowed:
        return min(rows,key=lambda r:float(r["loss"]))
    floor=min(allowed,key=lambda r:float(r["loss"]))
    stable=[r for r in allowed if float(r["stationarity"])<=STATIONARY and
            likelihood_ratio(float(r["loss"]),float(floor["loss"]),n)<=.001]
    return min(stable,key=lambda r:float(r["loss"])) if stable else floor


def audit():
    k=model()
    c=Config(k)
    base=read(COUNTS/"real_best.tsv")[0]
    x=c.load_old(base)
    rows=[]
    for rho in (0.,.23,.806,.965):
        x[-1]=rho/RHO_MAX
        p,j=c.unmap(x)
        f,g=k.metrics(p,k.y,True)
        err=k.predict_logs(p)-k.y
        R=k.noise_covariance(rho)
        sign,d=np.linalg.slogdet(R)
        q=err@np.linalg.solve(R,err)
        assert sign==1
        obj=.5*q*np.exp(d/k.n)
        assert np.isclose(obj,f[0],rtol=5e-13)
        assert np.isclose(q,f[2],rtol=5e-13) and np.isclose(d,f[3],atol=1e-11)
        # One-sided noise derivative at zero, and feasible interior central
        # movements in all normalized coordinates for the gradient audit.
        xx=np.clip(x,3e-5,1-3e-5)
        p2,j2=c.unmap(xx)
        analytic=c.gradient(k.metrics(p2,k.y,True)[1],j2)
        h=1e-6
        numeric=[(k.metrics(c.unmap(xx+h*q0)[0],k.y)[0][0]-
                  k.metrics(c.unmap(xx-h*q0)[0],k.y)[0][0])/(2*h) for q0 in np.eye(len(xx))]
        error=max(abs(analytic-np.array(numeric))/np.maximum(np.abs(analytic),1.))
        assert error<3e-5,(rho,error)
        rows.append(dict(rho=rho,log_count_rows=k.n,passages=len(np.unique(k.groups)),
                         transitions=len(k.transition),rho_upper=RHO_MAX,
                         dense_loss_discrepancy=abs(obj-f[0]),max_scaled_gradient_error=error))
    write(HERE/"likelihood_validation.tsv",rows)
    stage1=HERE.parent.parent/"structural_identifiability/20261118_culture_experimental_counts_local_01"
    inputs=(COUNTS/"count_inputs.tsv",COUNTS/"real_best.tsv",COUNTS/"source_provenance.tsv",
            COUNTS/"profile_intervals.tsv",stage1/"exact_points.tsv",
            COUNTS/"count_kernel.cpp",COUNTS/"count_model.py")
    write(HERE/"input_hashes.tsv",[dict(path=str(p),sha256=sha256(p.read_bytes()).hexdigest()) for p in inputs])
    print("AR likelihood audited",k.n,114,max(r["dense_loss_discrepancy"] for r in rows),
          max(r["max_scaled_gradient_error"] for r in rows),flush=True)


if __name__=="__main__":
    audit()
