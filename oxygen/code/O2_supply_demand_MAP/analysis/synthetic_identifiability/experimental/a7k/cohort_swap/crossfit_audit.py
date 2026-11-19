"""Check whether an available optimum is a better start on the reciprocal target."""
from concurrent.futures import as_completed
import argparse

import numpy as np

from cohort_swap import LAYOUT, MODES, OUT, archived, case, directory, load, pool, read, row_key, score, sample, task, work, write


def scan(mode):
    archived_rows=archived(mode)
    k,c,_=case(mode)
    hyb={(r["branch"],int(r["case_id"]),int(r["replicate"])):r
         for phase in ("pilot","expanded") for r in load(mode,phase)}
    assert len(hyb)==400
    all_targets={}
    for phase in ("pilot","expanded"):
        path=directory(mode,phase)
        manifest=read(path/"manifest.tsv")
        for r,y in zip(manifest,np.load(path/"frozen.npz")["target"]):
            all_targets[r["branch"],int(r["case_id"]),int(r["replicate"])]=y
    patches=[]
    checks=[]
    jobs=[]
    labels={}
    for j in (1,2):
        for i in range(1,101):
            rows=(archived_rows[j,i,"up","relaxed"],hyb["2N",j,i],
                  hyb["4N",j,i],archived_rows[j,i,"down","relaxed"])
            ys=(sample(mode,j,i,"up"),all_targets["2N",j,i],all_targets["4N",j,i],
                sample(mode,j,i,"down"))
            xs=[c.anchor(r) for r in rows]
            ps=[c.unmap(x)[0] for x in xs]
            for q,(row,y) in enumerate(zip(rows,ys)):
                vals=[k.parts(p,y)[1] for p in ps]
                assert abs(vals[q]-float(row["loss"]))<1e-7
                source=int(np.argmin(vals))
                gain=vals[q]-vals[source]
                cell=list(LAYOUT)[q]
                checks.append(dict(mode=mode,case_id=j,replicate=i,configuration=cell,
                    gain_from_cross_evaluation=gain,best_available_source=list(LAYOUT)[source]))
                if gain<=.025:
                    continue
                m=task(mode,j,i,"down" if cell=="LL" else "up","crosscheck",cell)
                med=np.clip(xs[source]+np.random.default_rng(i+519).normal(0.,.00014,len(xs[source])),0.,1.)
                jobs.append((m,y,xs[source],med,np.full(len(med),np.nan)))
                labels[row_key(m)]=(gain,row,vals[source],y)
    write(OUT/f"crossfit_checks_{mode}.tsv",checks)
    if not jobs:
        print(mode,"800 cross-configuration checks, no lower-deviance candidate",flush=True)
        return
    with pool(min(8,len(jobs))) as ex:
        for f in as_completed(ex.submit(work,q) for q in jobs):
            r=f.result()
            gain,previous,check,y=labels[row_key(r)]
            assert r["validated"] and float(r["loss"])<=check+.005
            p=c.unmap(c.anchor(r))[0]
            assert abs(k.parts(p,y)[1]-float(r["loss"]))<1e-7
            j=int(r["case_id"])
            a,b=(score(mode,j,v) for v in (previous,r))
            patches.append(dict(r,configuration=r["branch"],prior_error=a,new_error=b,
                                prior_loss=previous["loss"],changed=(abs(a)<=.1)!=(abs(b)<=.1)))
            print(mode,r["case_id"],r["replicate"],r["branch"],
                  float(previous["loss"])-float(r["loss"]),a,b,patches[-1]["changed"],flush=True)
    write(OUT/f"crossfit_rechecks_{mode}.tsv",patches)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=MODES)
    scan(p.parse_args().mode)
