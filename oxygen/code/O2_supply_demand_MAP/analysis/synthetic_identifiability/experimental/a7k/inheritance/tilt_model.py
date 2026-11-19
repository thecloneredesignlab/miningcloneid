"""A7K composition sensitivity with the frozen AR(1) and karyotype model."""
import ctypes as C
from functools import lru_cache
import os
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import nnls
from scipy.special import xlogy

HERE=Path(__file__).resolve().parent
PREV=HERE.with_name("20261119_experimental_death_A7K_comparison_01")
sys.path.insert(0,str(PREV))
os.environ["CULTURE_DEATH_A7K_LIB"]=os.environ.get(
    "CULTURE_A7K_TV_LIB","/private/tmp/culture_a7k_tv_20261119.dylib")
from modal import Model as FrozenModel
from modal import BASE, CountConfig, EPS, FD, FI, GRID, THRESHOLD, pool, read, row_x, solve, write

TARGETS=("p_misseg","buffer_beta")
GROUP=("group_010","group_025","group_050")
BRANCH=("branch_010","branch_025","branch_050")
RHO=54
N_ORIGINAL_CONSTRAINTS=110
STATIONARY=.22


def cap(mode):
    return 0. if mode=="strict" else int(mode.split("_")[1])/1000.


class TiltModel(FrozenModel):
    def __init__(self,mode):
        super().__init__("stock_a7k")
        self.sensitivity=mode
        self.cap=cap(mode)
        self.lib.extension_enable_a7k_tilts.argtypes=[C.c_void_p]
        self.lib.extension_a7k_composition.argtypes=[C.c_void_p,FD,FD,C.c_void_p,
                                                     C.c_void_p,C.c_void_p]
        self.lib.extension_tv_adjoint.argtypes=[C.c_void_p,FD,C.c_int,FI,FD]
        assert self.lib.extension_enable_a7k_tilts(self.handle)

    def composition(self,p,detail=False):
        tv=np.empty(4)
        if detail:
            before=np.empty((4,133)); after=np.empty_like(before); seed=np.empty(4)
            pointers=[arr.ctypes.data for arr in (before,after,seed)]
        else:
            pointers=[None]*3
        ok=self.lib.extension_a7k_composition(self.handle,p,tv,*pointers)
        if not ok:
            if detail:
                raise ArithmeticError("Nonfinite A7K composition")
            return np.full(4,np.inf)
        return (tv,before,after,seed) if detail else tv

    def supply(self,p,gaussian=False,pick=None):
        original=np.empty(N_ORIGINAL_CONSTRAINTS)
        if pick is None:
            ids=np.empty(0,dtype=np.int32)
            grad=None
        else:
            ids=np.ascontiguousarray(pick,dtype=np.int32)
            grad=np.zeros((len(ids),55))
        ordinary=ids<N_ORIGINAL_CONSTRAINTS
        take=np.ascontiguousarray(ids[ordinary])
        raw=np.empty((len(take),54)) if len(take) else None
        ok=self.lib.ols_supply(self.handle,p,False,original,len(take),
                               take.ctypes.data if len(take) else None,
                               None if raw is None else raw.ctypes.data)
        if not ok:
            original.fill(-10.)
        margin=original
        if self.cap:
            margin=np.r_[original,1.-self.composition(p)/self.cap]
        if grad is not None and ok:
            if len(take):
                grad[ordinary,:54]=raw
            tv_select=np.ascontiguousarray(ids[~ordinary]-N_ORIGINAL_CONSTRAINTS)
            if len(tv_select):
                tv_grad=np.empty((len(tv_select),54))
                if not self.lib.extension_tv_adjoint(self.handle,p,len(tv_select),tv_select,tv_grad):
                    raise ArithmeticError("Nonfinite TV gradient")
                grad[~ordinary,:54]=-tv_grad/self.cap
        return margin,grad

    def parts(self,p,target,gradient=False):
        ar,_=super().metrics(p,target[:self.n])
        if ar[0]>=1e12:
            return ar,1e13,np.zeros(len(p)) if gradient else None,np.nan,np.nan
        channels,states=self.readouts(p)
        total=states.sum(axis=1,keepdims=True)
        probabilities=(1.-EPS)*states/total+EPS/133.
        data=target[self.n:].reshape(self.kary.shape)
        nll=-float(xlogy(data,probabilities).sum())
        saturated=-float(xlogy(data,data/data.sum(axis=1,keepdims=True)).sum())
        loss=self.n*np.log(ar[0])+2.*(nll-saturated)
        base_data=data[:2]
        base_sat=-float(xlogy(base_data,base_data/base_data.sum(axis=1,keepdims=True)).sum())
        base_loss=self.n*np.log(ar[0])+2.*(-float(xlogy(base_data,probabilities[:2]).sum())-base_sat)
        g=None
        if gradient:
            w=-(1.-EPS)*data/probabilities
            w=np.ascontiguousarray((w-(w*states).sum(axis=1,keepdims=True)/total)/total)
            g=np.zeros(55)
            if not self.lib.extension_deviance_adjoint(
                    self.handle,p,len(self.kary_rows),self.kary_rows,np.zeros_like(channels),2.*w,
                    self.n,self.indices,self.previous,self.gaps,target[:self.n],ar[2],g):
                raise ArithmeticError("Could not calculate joint gradient")
        return ar,loss,g,nll,base_loss

    def constrained(self,p,y,c,gradient=True):
        ar,f,g,nll,base=self.parts(p,y,gradient)
        margin=ar[1]
        if self.cap:
            margin=min(margin,float(np.min(1.-self.composition(p)/self.cap)))
        return f,g,margin


class Config:
    def __init__(self,k):
        self.base=CountConfig(k)
        self.k=k
        self.names=self.base.names.copy()
        self.active=self.base.active.copy()
        self.active[-1]=RHO
        self.dim=len(self.names)
        self.domain=12.+400.*k.cap
        self.kary_labels=[a["physical_passage"].split("NLS_")[1].split("_A7K")[0] for a in k.a7]
        self.centers=np.array([-1. if a["cohort"]=="2N" else 0. for a in k.a7])
        if not k.cap:
            self.groups=np.zeros(0,dtype=int)
            return
        if k.sensitivity.startswith("group_"):
            self.groups=np.array([int(a["cohort"]=="4N") for a in k.a7])
            labels=["2N","4N"]
        else:
            self.groups=np.arange(4)
            labels=self.kary_labels
        self.names += [f"tilt_{effect}_{g}" for g in labels
                       for effect in ("mean","shape")]

    def anchor(self,r):
        return np.r_[row_x(self.base,r),[float(r.get(f"x_{name}",.5))
                                       for name in self.names[self.dim:]]]

    def unmap(self,x):
        a,j=self.base.unmap(x[:self.dim])
        p=np.r_[a[:46],np.zeros(8),a[-1]]
        if len(self.groups):
            shifts=(self.domain*(2.*x[self.dim:]-1.)).reshape(-1,2)[self.groups,:].copy()
            # Center the two basis terms near the nominal donor separately by cohort.
            shifts[:,0]-=2.*self.centers*shifts[:,1]
            p[46:54]=shifts.ravel()
        return p,j

    def gradient(self,g,j):
        v=self.base.gradient(np.r_[g[:46],g[-1]],j)
        if not len(self.groups):
            return v
        extra=g[46:54].reshape(4,2).copy()
        extra[:,1]-=2.*self.centers*extra[:,0]
        accum=np.zeros((max(self.groups)+1,2))
        np.add.at(accum,self.groups,extra)
        return np.r_[v,2.*self.domain*accum.ravel()]


@lru_cache
def case(mode):
    k=TiltModel(mode)
    return k,Config(k),k.target()


def receipt(k,c,y,f,**labels):
    x=f["x"]
    p,_=c.unmap(x)
    ar,loss,g,nll,base=k.parts(p,y)
    tv,*_=k.composition(p,True)
    margin=min(ar[1],float(np.min(1.-tv/k.cap))) if k.cap else ar[1]
    return dict(labels,loss=loss,count_loss=float(ar[0]),kary_nll=nll,
                common_deviance=base,a7k_deviance=loss-base,
                donor_margin=float(ar[1]),all_constraint_margin=float(margin),
                stationarity=f["gradient"],rho=float(p[-1]),sigma=np.sqrt(ar[2]/k.n),
                success=f["success"],seconds=f["seconds"],calls=f.get("calls",""),
                **{f"x_{name}":v for name,v in zip(c.names,x)},
                **{n:p[j] for n,j in zip(c.names[:c.dim],c.active) if n!="rho"},
                **{f"tv_{n}":t for n,t in zip(c.kary_labels,tv)})


def stable(r):
    return float(r["all_constraint_margin"])>=0 and float(r["stationarity"])<=STATIONARY


def take(rows):
    viable=[r for r in rows if float(r["all_constraint_margin"])>=0.]
    if not viable:
        return min(rows,key=lambda r:float(r["loss"]))
    best=min(viable,key=lambda r:float(r["loss"]))
    matched=[r for r in viable if stable(r) and float(r["loss"])<=float(best["loss"])+.01]
    return min(matched or viable,key=lambda r:float(r["loss"]))


def checked_fit(k,c,y,start,**options):
    f=solve(k,c,y,start,**options)
    if not k.cap or f["supply"]<0 or not STATIONARY<f["gradient"]<30.:
        return f
    # The total-variation constraint can bind where one class changes sign.
    # Add the two one-sided normals at such points to the KKT cone check.
    x=f["x"]
    fixed=options.get("fixed") or {}
    active=np.array([j for j in range(len(x)) if j not in fixed])
    p,jac=c.unmap(x)
    margins,_=k.supply(p)
    binding=np.flatnonzero(margins<=5e-7)
    tv=binding[binding>=N_ORIGINAL_CONSTRAINTS]
    if not len(tv):
        return f
    _,old=k.supply(p,pick=binding)
    normals=[c.gradient(z,jac)[active] for z in old]
    for i,j in enumerate(active):
        if x[j]<5e-7 or x[j]>1.-5e-7:
            a=np.zeros(len(active)); a[i]=1. if x[j]<5e-7 else -1.
            normals.append(a)
    for index in tv:
        branch=index-N_ORIGINAL_CONSTRAINTS
        start=c.dim+2*c.groups[branch]
        for direction in ((-1.,0.),(1.,0.),(0.,-1.),(0.,1.),(-1.,-1.),(1.,1.),(-1.,1.),(1.,-1.)):
            z=x.copy()
            z[start:start+2]=np.clip(z[start:start+2]+3e-7*np.array(direction),0.,1.)
            p1,j1=c.unmap(z)
            m,g=k.supply(p1,pick=[index])
            if abs(m[index]-margins[index])<1e-4:
                u=c.gradient(g[0],j1)[active]
                if min(np.linalg.norm(u-a) for a in normals)>1e-3:
                    normals.append(u)
    g=c.gradient(k.parts(p,y,True)[2],jac)[active]
    try:
        norm=nnls(np.asarray(normals).T,g,maxiter=2000)[1]
    except RuntimeError:
        return f
    f["gradient"]=min(f["gradient"],float(norm))
    return f


def run_fit(job):
    mode,origin,x,fixed,limit,divisor=job
    k,c,y=case(mode)
    fit=checked_fit(k,c,y,x,fixed=fixed,limit=limit,divisor=divisor)
    return receipt(k,c,y,fit,mode=mode,origin=origin)
