"""Release an improved common reference only after multiple exact OLS refits."""
import concurrent.futures as cf

import numpy as np

from practical import HERE, Configuration, read, write, row_x, solve_constrained, result_row
from inference import pool
from fit_stage import model


def trial(job):
    origin,x=job
    k=model()
    c=Configuration(k)
    r=solve_constrained(k,c,k.y,x,limit=80)
    return result_row(c,r,case="main",origin=origin)


if __name__=="__main__":
    k=model()
    c=Configuration(k)
    old=read(HERE/"main_constrained_best.tsv")[0]
    rows=read(HERE/"main_profiles_envelope.tsv")
    better=sorted([s for s in rows if s["feasible"]=="True" and float(s["loss"])<float(old["loss"])-1e-6],
                  key=lambda s:float(s["loss"]))
    assert better
    rnd=np.random.default_rng(1804)
    x=row_x(c,better[0])
    work=[("profile_1",x),("profile_2",np.clip(x+rnd.normal(0,.001,len(x)),0,1)),
          ("profile_3",row_x(c,better[1]))]
    result=[]
    with pool(3) as p:
        for task in cf.as_completed([p.submit(trial,r) for r in work]):
            res=task.result()
            result.append(res);write(HERE/"profile_reference_refits.tsv",result)
            print(res["origin"],res["loss"],res["max_projected_gradient"],res["min_aggregate_donor_margin_fraction"],flush=True)
    feasible=[s for s in result if s["min_aggregate_donor_margin_fraction"]>=0 and s["max_projected_gradient"]<=.025]
    best=min(feasible,key=lambda a:a["loss"],default=None)
    if best is not None and float(old["loss"])>best["loss"]:
        write(HERE/"main_constrained_best.tsv",[best])
        values,_,_,_,_=k.predict(c.unmap(row_x(c,best))[0])
        write(HERE/"real_residuals_constrained.tsv",
              [dict(passage_id=i,observed=y,predicted=b,residual=y-b,variance_multiplier=v)
               for i,y,b,v in zip(k.ids,k.y,values,k.var)])
    elif min(s["loss"] for s in result)<float(old["loss"])-1e-4:
        raise RuntimeError("improved reference not yet revalidated")
