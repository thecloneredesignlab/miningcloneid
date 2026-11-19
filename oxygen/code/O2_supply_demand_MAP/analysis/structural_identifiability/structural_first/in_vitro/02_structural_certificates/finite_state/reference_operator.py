"""Independent entrywise derivative reconstruction, not the signed-event loop."""
from math import comb
from flint import arb
from algebra import rates
from protocol import MECHANISMS

def entry(p,O,destination,mother,mask):
    (lam,mu,prob,logs),grad=rates(p,arb(mother)/44,O)
    dl,dm,dp,ds=grad;w=p['p_wgd'];m=abs(destination-mother)
    T=arb(0);dT=[arb(0) for name in MECHANISMS]
    if m<=mother and mask[m]=='1':
        pi=arb(comb(mother,m))*prob**m*(1-prob)**(mother-m)
        assert pi>arb('1e-8')
        T=(2 if m==0 else 1)*pi*((m*logs).exp() if m else 1)
        factor=m/prob-(mother-m)/(1-prob)
        dT=[T*(factor*dp[name]+m*ds[name]) for name in MECHANISMS]
    whole=int(destination==2*mother);diagonal=int(destination==mother)
    value=lam*((1-w)*T+w*whole-diagonal)
    gradient=[dl[name]*((1-w)*T+w*whole-diagonal)+lam*((1-w)*dT[j]+int(name=='p_wgd')*(whole-T)) for j,name in enumerate(MECHANISMS)]
    return value,gradient
