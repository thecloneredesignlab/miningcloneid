"""Fixed-calendar comparison of existing death proportions and A7K karyotypes."""
from functools import lru_cache
import ctypes as C
import os
from pathlib import Path
import sys

import numpy as np
from scipy.special import betaln, digamma, expit, xlogy

HERE=Path(__file__).resolve().parent
RESULT=HERE.parent
STAGE3=RESULT/"20261119_experimental_counts_karyotype_practical_01"
PLAN=RESULT/"20261119_experimental_A19_followup_plan_01"
ALIGN=RESULT/"20261119_experimental_karyotype_alignment_01"
BASE=RESULT/"20261118_experimental_counts_ar1_practical_01"
sys.path.insert(0,str(BASE))
os.environ["CULTURE_AR1_LIB"]=os.environ.get(
    "CULTURE_DEATH_A7K_LIB","/private/tmp/culture_death_A7K_20261119.dylib")
from ar1_model import (CalendarAR1, Config as CountConfig, FD, FI, THRESHOLD, STATIONARY,
                       likelihood_ratio, pool, read, row_x, solve, write)
from practical import base_table

TARGETS=("p_misseg","buffer_beta","k_o_mis")
GRID=np.arange(22,155)
EPS=1e-12
DEATH_NAMES=("class_log_odds","shared_count_loading","misseg_dead_share","log10_concentration")
LOW=np.array([-4.,-4.,0.,0.])
SPAN=np.array([8.,8.,1.,6.])
MODES=("baseline","a7k","death","death_total","death_legacy","a7k_death",
       "stock_baseline","stock_a7k","stock_death","stock_death_total","stock_death_legacy",
       "stock_death_wide","stock_death_logwide")


class Model(CalendarAR1):
    def __init__(self,mode):
        super().__init__()
        self.mode=mode
        stock=mode.startswith("stock_")
        mode=mode.removeprefix("stock_")
        _,_,_,_,slots,weights,_,_,ox,meta,amounts,hist,_=base_table()
        indices={r["passage_id"]:i for i,r in enumerate(self.graph)}
        slot_index={r["record_id"]:i for i,r in enumerate(self.slots)}
        self.old_slots=len(slots)
        self.a7=read(PLAN/"A7K_matched_readouts.tsv")
        for a in self.a7:
            donor=self.slots[slot_index[a["donor_record"]]]
            assert float(donor["oxygen_pct"])==0. and int(donor["calendar_day"])==5
            o=int(np.flatnonzero(ox==float(a["oxygen_pct"]))[0])
            meta=np.vstack([meta,[o,int(a["cohort"]=="4N"),indices[donor["passage_id"]],
                                  int(donor["euler_step"]),80,-1]])
            amounts=np.r_[amounts,float(a["fixed_recorded_seed"])]
            slots=np.vstack([slots,[len(meta)-1,80,0]])
            weights=np.r_[weights,0.]
            slot_index[a["assay_origin"]]=len(slots)-1
        cells=read(ALIGN/"cell_alignment.tsv")
        sample=read(ALIGN/"sample_alignment.tsv")

        def histogram(origin):
            n=np.array([int(r["observed_N_22_autosomes"]) for r in cells
                        if r["assay_origin_id"]==origin])
            assert len(n)>0 and min(n)>=22 and max(n)<=154
            return np.bincount(n-22,minlength=133).astype(float)

        for r in sample:
            if r["role"]=="proxy_only":
                hist[int(r["cohort"]=="4N")]=histogram(r["assay_origin_id"])/20.
        ids=[r["assay_origin_id"] for r in sample if r["role"]=="direct_harvest_anchor"]
        self.a7_indices=np.ascontiguousarray([slot_index[r["assay_origin"]] for r in self.a7],
                                              dtype=np.int32)
        self.all_kary=ids+[r["assay_origin"] for r in self.a7]
        if mode=="a7k":
            ids=self.all_kary
        self.kary=np.ascontiguousarray([histogram(v) for v in ids])
        self.kary_rows=np.ascontiguousarray([slot_index[v] for v in ids],dtype=np.int32)
        self.all_kary_rows=np.ascontiguousarray([slot_index[v] for v in self.all_kary],dtype=np.int32)
        self.width=0.
        self.lib.ols_destroy.argtypes=[C.c_void_p]
        self.lib.ols_destroy(self.handle)
        self.handle=self.lib.ols_design(len(meta),len(ox),self.ols_n,ox,
                     np.ascontiguousarray(meta,dtype=np.int32),amounts,len(slots),
                     np.ascontiguousarray(slots,dtype=np.int32),weights,np.ascontiguousarray(hist))
        assert self.handle and self.lib.ols_nconstraints(self.handle)==110
        if stock:
            self.lib.extension_keep_original_supply.argtypes=[C.c_void_p,C.c_int]
            assert self.lib.extension_keep_original_supply(self.handle,4)
        self.nslots=len(slots)
        self.lib.extension_readouts.argtypes=[C.c_void_p,FD,C.c_int,FI,FD,FD]
        self.lib.extension_deviance_adjoint.argtypes=[C.c_void_p,FD,C.c_int,FI,FD,FD,
            C.c_int,FI,FI,FI,FD,C.c_double,FD]
        records=read(PLAN/"paired_death_readout_audit.tsv")
        if mode=="death_legacy":
            records=[r for r in records if r["in_legacy_time_window"]=="TRUE"]
        records=[r for r in records if r["classification_status"]=="classified_full_promoted"]
        assert len(records)==(149 if mode=="death_legacy" else 191)
        if mode=="a7k_death":
            ids=[r["assay_origin"] for r in self.a7]
            d=[int(r["dead_objects"]) for r in self.a7]
            n=[int(r["eligible_objects"]) for r in self.a7]
        else:
            ids=[r["record_id"] for r in records] if mode.startswith("death") else []
            d=[int(r["dead_count"]) for r in records] if ids else []
            n=[int(r["eligible_objects"]) for r in records] if ids else []
        self.death_rows=np.array([slot_index[i] for i in ids],dtype=int)
        index_of_capture={r:i for i,r in enumerate(self.indices)}
        self.death_shared=np.array([i in index_of_capture for i in self.death_rows],dtype=bool)
        self.death_log_index=np.array([index_of_capture.get(i,0) for i in self.death_rows],dtype=int)
        self.eligible=np.array(n,dtype=float)
        self.dead=np.array(d,dtype=float)
        self.use_death=bool(ids)
        self.kary_size=self.kary.size

    def supply(self,p,gaussian=False,pick=None):
        g,j=super().supply(p[:47],False,pick)
        if len(p)>47:
            j=np.pad(j,[(0,0)]*(j.ndim-1)+[(0,len(p)-47)])
        return g,j

    def readouts(self,p,active=None):
        active=self.kary_rows if active is None else active
        channels=np.empty((3,self.nslots))
        states=np.empty((len(active),133))
        if not self.lib.extension_readouts(self.handle,p,len(active),active,channels,states):
            raise ArithmeticError("Nonfinite trajectory")
        return channels,states

    def predict_logs(self,p):
        ch,_=self.readouts(p)
        return np.log(np.maximum(ch[:,self.indices].sum(axis=0)*1e12,1.))

    def target(self,log_counts=None,kary=None,dead=None):
        return np.ascontiguousarray(np.r_[self.y if log_counts is None else log_counts,
                  (self.kary if kary is None else kary).ravel(),
                  self.dead if dead is None else dead])

    def death_prediction(self,p,channels,target):
        d=self.death_rows
        live,hyp,mis=channels[:,d]
        intercept,loading,q,logconc=p[47:51]
        dead=hyp+q*mis+EPS
        other=live+(1.-q)*mis+EPS
        count=np.maximum((live+hyp+mis)*1e12,1.)
        residual=(target[self.death_log_index]-np.log(count))*self.death_shared
        odds=intercept+loading*residual+np.log(dead/other)
        probability=(1.-2.*EPS)*expit(odds)+EPS
        return probability,10.**logconc,dead,other,residual

    def parts(self,p,target,gradient=False):
        ar,_=super().metrics(p,target[:self.n])
        if ar[0]>=1e12:
            return ar, 1e13, np.zeros(len(p)) if gradient else None, np.nan, np.nan
        channels,states=self.readouts(p)
        total=states.sum(axis=1,keepdims=True)
        probabilities=(1.-EPS)*states/total+EPS/133.
        data=target[self.n:self.n+self.kary_size].reshape(self.kary.shape)
        nll_k=-float(xlogy(data,probabilities).sum())
        constants=-float(xlogy(data,data/data.sum(axis=1,keepdims=True)).sum())
        nll_d=0.
        forcing=np.zeros_like(channels)
        grad=np.zeros(len(p)) if gradient else None
        if self.use_death:
            obs=target[self.n+self.kary_size:]
            prob,kappa,dead,other,residual=self.death_prediction(p,channels,target)
            a=kappa*prob
            b=kappa*(1.-prob)
            n=self.eligible
            nll_d=float(np.sum(betaln(a,b)-betaln(obs+a,n-obs+b)))
            constants -= float(np.sum(xlogy(obs,obs/n)+xlogy(n-obs,1.-obs/n)))
            if gradient:
                ga=digamma(a)-digamma(a+obs)+digamma(kappa+n)-digamma(kappa)
                gb=digamma(b)-digamma(b+n-obs)+digamma(kappa+n)-digamma(kappa)
                dlog=(ga-gb)*kappa*(prob-EPS)*(1.-EPS-prob)/(1.-2.*EPS)
                slot=self.death_rows
                live,hyp,mis=channels[:,slot]
                loading,q=p[48:50]
                combined=live+hyp+mis
                common=-(loading*self.death_shared*(combined>EPS))/np.maximum(combined,EPS)
                forcing[:,slot]=dlog*np.vstack([common-1./other,common+1./dead,
                                               common+q/dead-(1.-q)/other])
                grad[47:51]=[dlog.sum(),dlog@residual,dlog@(mis/dead+mis/other),
                               np.log(10.)*kappa*np.sum(ga*prob+gb*(1.-prob))]
        loss=self.n*np.log(ar[0])+2.*(nll_k+nll_d-constants)
        if gradient:
            w=-(1.-EPS)*data/probabilities
            w=np.ascontiguousarray((w-(w*states).sum(axis=1,keepdims=True)/total)/total)
            if not self.lib.extension_deviance_adjoint(
                    self.handle,p,len(self.kary_rows),self.kary_rows,2.*forcing,2.*w,
                    self.n,self.indices,self.previous,self.gaps,target[:self.n],ar[2],grad[:47]):
                raise ArithmeticError("Could not evaluate the joint adjoint")
            grad[47:]*=2.
        return ar, loss, grad, nll_k, nll_d

    def constrained(self,p,y,c,gradient=True):
        ar,f,g,_,_=self.parts(p,y,gradient)
        return f,g,ar[1]


class Config:
    def __init__(self,k):
        self.base=CountConfig(k)
        self.k=k
        self.death=k.use_death
        self.death_low=LOW.copy()
        self.death_span=SPAN.copy()
        self.share_scale=0.
        self.log_share=k.mode=="stock_death_logwide"
        if k.mode in ("stock_death_wide","stock_death_logwide"):
            self.death_low[0]=-8.
            self.death_span[0]=16.
        if self.log_share:
            self.death_low[2]=-8.
            self.death_span[2]=8.
        self.names=self.base.names.copy()
        self.active=self.base.active.copy()
        if self.death:
            self.names+=list(DEATH_NAMES)
            self.active=np.r_[self.active,np.arange(47,51)]

    def anchor(self,r):
        out=row_x(self.base,r)
        if self.death:
            dx=np.array([float(r.get(f"x_{k}",v)) for k,v in zip(DEATH_NAMES,(.5,.5,.55,.33))])
            if all(k in r for k in DEATH_NAMES):
                d=np.array([float(r[k]) for k in DEATH_NAMES])
                if self.log_share:
                    d[2]=np.log10(max(d[2],1e-8))
                dx=np.clip((d-self.death_low)/self.death_span,0.,1.)
                if self.share_scale:
                    dx[2]=self.encode_share(d[2])
            out=np.r_[out,dx]
        return out

    def encode_share(self,q):
        return np.log1p(q/self.share_scale)/np.log1p(1./self.share_scale)

    def unmap(self,x):
        n=len(self.base.names)
        p,j=self.base.unmap(x[:n])
        if self.death:
            d=self.death_low+self.death_span*x[n:n+4]
            derivative=self.death_span.copy()
            if self.log_share:
                d[2]=10.**d[2]
                derivative[2]*=np.log(10.)*d[2]
            if self.share_scale:
                rate=np.log1p(1./self.share_scale)
                d[2]=self.share_scale*np.expm1(rate*x[n+2])
                derivative[2]=rate*(d[2]+self.share_scale)
            p=np.r_[p,d]
            j=np.r_[j,derivative]
        return p,j

    def gradient(self,g,j):
        return (np.r_[self.base.gradient(g[:47],j[:-4]),g[47:]*j[-4:]] if self.death
                else self.base.gradient(g,j))


@lru_cache
def case(mode):
    k=Model(mode)
    c=Config(k)
    return k,c,k.target()


def receipt(k,c,y,fit,**labels):
    x=fit["x"]
    p,j=c.unmap(x)
    ar,f,g,ky,de=k.parts(p,y)
    return dict(labels,loss=f,count_loss=ar[0],kary_nll=ky,death_nll=de,
        donor_margin=ar[1],stationarity=fit["gradient"],rho=p[46],sigma=np.sqrt(ar[2]/k.n),
        success=fit["success"],seconds=fit["seconds"],
        **{f"x_{n}":v for n,v in zip(c.names,x)},
        **{n:p[j] for n,j in zip(c.names,c.active) if n!="rho"})


def stable(row):
    return float(row["donor_margin"])>=0. and float(row["stationarity"])<=.22


def fixed_for(mode,c):
    # Fully classify the missegregation-associated object as dead only in the
    # specified alternate observation model.
    if mode.removeprefix("stock_")=="death_total":
        return {c.names.index("misseg_dead_share"):1.}
    if mode=="a7k_death":
        return {c.names.index("shared_count_loading"):.5,c.names.index("class_log_odds"):.5,
                c.names.index("misseg_dead_share"):1.,c.names.index("log10_concentration"):1./3}
    return {}
