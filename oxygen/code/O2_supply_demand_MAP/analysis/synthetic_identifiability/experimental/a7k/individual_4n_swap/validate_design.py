"""Reproduce endpoint observations and verify the split model nests the shared model."""
import argparse
from pathlib import Path

import numpy as np

from branch_swap import BACKS, MODES, OUT, bounds, data, model, write


def audit(mode):
    k,shared,_=model(mode,"shared")
    _,split,_=model(mode,"split4N")
    archived=bounds(mode)
    output=[]
    for b in BACKS:
        for j in (1,2):
            for q in ("HH","LL"):
                err=[]
                for i in range(1,101):
                    row=archived[b,j,i,q]
                    y,_=data(mode,b,j,i,q)
                    x=shared.anchor(row)
                    z=split.anchor(row)
                    a=shared.unmap(x)[0]
                    c=split.unmap(z)[0]
                    assert np.array_equal(a,c)
                    dev=k.parts(a,y)[1]
                    discrepancy=abs(dev-float(row["loss"]))
                    assert discrepancy<1e-7
                    err.append(discrepancy)
                output.append(dict(mode=mode,background=b,case_id=j,cell=q,n=100,
                    max_difference=max(err),split_embedding_equal=True))
    y,_=data(mode,"up",1,1,"HH")
    x=split.anchor(archived["up",1,1,"HH"])
    p,jac=split.unmap(x)
    analytic=split.gradient(k.parts(p,y,True)[2],jac)
    gaps=[]
    for q in range(split.dim,len(x)):
        v=np.zeros(len(x))
        v[q]=1e-6
        left,right=(k.parts(split.unmap(x+sign*v)[0],y)[1] for sign in (-1,1))
        estimate=(right-left)/2e-6
        gaps.append(abs(estimate-analytic[q])/(1.+abs(estimate)))
    assert max(gaps)<5e-6
    write(OUT/f"design_replay_{mode}.tsv",output)
    print(mode,"800 equivalent control fits, max reproducibility error",
          max(r["max_difference"] for r in output),
          "maximum split tilt derivative check",max(gaps),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=MODES)
    audit(p.parse_args().mode)
