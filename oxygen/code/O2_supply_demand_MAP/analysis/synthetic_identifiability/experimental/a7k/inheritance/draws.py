"""Archived strict draws and paired bounded A7K tilt challenges."""
import numpy as np
from scipy.optimize import brentq

from tilt_model import BASE,PREV,GRID,case,read

Z=(GRID-88.)/44.


def coefficients(p_before,group,delta,direction):
    """Choose a shared linear tilt reaching the cap on at least one sibling."""
    pair=p_before[group]
    def max_change(value):
        logweight=direction*value*Z
        weights=np.exp(logweight-logweight.max())
        z=pair*weights
        z=z/z.sum(axis=1,keepdims=True)
        return np.max(np.abs(z-pair).sum(axis=1)*.5)
    right=1.
    while max_change(right)<delta:
        right*=2.
        if right>128.:
            raise ArithmeticError("The requested tilt cannot reach its cap safely.")
    return direction*brentq(lambda value:max_change(value)-delta*(1-1e-8),0.,right)


def draw(group_mode,j,i,scenario="strict"):
    k,c,_=case("strict")
    truth=read(BASE/"truths.tsv")[j-1]
    x=c.anchor(truth)
    p,_=c.unmap(x)
    rho=float(truth["rho"]); sigma=float(truth["sigma"])
    rng=np.random.default_rng(np.random.SeedSequence([20261118,221,j,i]))
    z=rng.normal(size=k.n)
    noise=sigma*z
    for t in k.transition:
        a=rho**k.gaps[t]
        noise[t]=a*noise[t-1]+np.sqrt(1.-a*a)*sigma*z[t]
    logc=k.predict_logs(p)+noise
    tv=np.zeros(4)
    if scenario!="strict":
        q,_=case(group_mode)[:2]
        direction={"up":1.,"down":-1.}[scenario]
        _,before,_,_=k.composition(p,True)
        for g in (np.array([0,1]),np.array([2,3])):
            shift=coefficients(before,g,q.cap,direction)
            p[46+2*g]=shift
        tv=q.composition(p)
        assert 0.<np.max(tv)<=q.cap+1e-12
    _,states=k.readouts(p)
    pr=states/states.sum(axis=1,keepdims=True)
    rng=np.random.default_rng(np.random.SeedSequence([20261119,40,j,i]))
    sampled=np.stack([rng.multinomial(int(s),v) for s,v in zip(k.kary.sum(axis=1),pr)])
    return k.target(logc,sampled),x,tv
