"""Audit all simulated estimates and apply the registered paired recovery gate."""
import argparse
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

from modal import HERE,BASE,TARGETS,case,read,stable,write
from recover_modal import draw


def run(mode):
    k,c,_=case(mode)
    p=HERE/f"recovery_{mode}.tsv"
    assigned=read(p)
    cases={(int(r["case_id"]),int(r["replicate"])):r for r in assigned}
    assert set(cases)=={(j,i) for i in range(1,101) for j in (1,2)}
    assert len(assigned)==len(cases)
    eligible=[r["parameter"] for r in read(HERE/f"intervals_{mode}.tsv") if r["precision_gate"]=="profile_pass"]
    verified=Counter()
    recovered=Counter()
    audit=[]
    for (j,i),row in sorted(cases.items()):
        y,truth=draw(mode,j,i)
        x=c.anchor(row)
        natural,_=c.unmap(x)
        ar,d,*_=k.parts(natural,y)
        discrepancy=abs(d-float(row["loss"]))
        assert discrepancy<1e-7,(j,i,d,discrepancy)
        accepted=row["validated"]=="True"
        if accepted:
            assert ar[1]>=-3e-12 and float(row["stationarity"])<=.22
            assert float(row["lr_agreement"])<=.10
            verified[j]+=1
            verified[j,"pilot"]+=int(i<=20)
        for name in eligible:
            index=c.names.index(name)
            error=abs(x[index]-truth[index])
            ok=accepted and error<=.10
            if ok:
                recovered[j,name]+=1
                recovered[j,name,"pilot"]+=int(i<=20)
            audit.append(dict(mode=mode,case_id=j,replicate=i,parameter=name,
                deviance_error=discrepancy,minimum_supply=ar[1],validated=accepted,
                normalized_error=error,within_tolerance=ok))
    assert min(verified[j,"pilot"] for j in (1,2))>=18
    rows=[]
    for j in (1,2):
        for name in eligible:
            val=recovered[j,name]
            ci=binomtest(val,100).proportion_ci()
            rows.append(dict(mode=mode,parameter=name,case_id=j,allocated=100,
                         numerically_valid=verified[j],recovered=val,rate=val/100,
                         binomial_ci_95_lower=ci.low,binomial_ci_95_upper=ci.high,
                         pilot_numerically_valid=verified[j,"pilot"],pilot_recovered=recovered[j,name,"pilot"]))
    write(HERE/f"recovery_metrics_{mode}.tsv",rows)
    write(HERE/f"recovery_audit_{mode}.tsv",audit)
    print(rows,flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("modes",nargs="+")
    for mode in p.parse_args().modes:
        run(mode)
