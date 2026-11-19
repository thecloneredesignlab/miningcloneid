"""Promote less expensive feasible profile branches to unrestricted starts."""
from concurrent.futures import as_completed
import sys
import numpy as np

from modal import HERE,case,pool,read,write,stable
from fit import run_one,take


def run(mode):
    file=HERE/f"profile_points_{mode}.tsv"
    k,c,_=case(mode)
    base=read(HERE/f"best_{mode}.tsv")[0]
    rows=read(HERE/f"fits_{mode}.tsv")
    seed=[]
    seq=sorted((r for r in read(file) if float(r["donor_margin"])>=0
                and float(r["loss"])<float(base["loss"])+.2),key=lambda r:float(r["loss"]))
    for r in seq:
        x=c.anchor(r)
        if not seed or all(np.linalg.norm(x-v)>.04 for v in seed):
            seed.append(x)
        if len(seed)==5:
            break
    jobs=[(mode,f"profile_restart_{i}",z,{},125,34.) for i,z in enumerate(seed)]
    with pool(3) as executor:
        for f in as_completed(executor.submit(run_one,t) for t in jobs):
            rows.append(f.result())
    best=take(rows)
    assert stable(best)
    assert float(best["loss"])<=float(seq[0]["loss"])+.015
    write(HERE/f"fits_{mode}.tsv",rows)
    write(HERE/f"best_{mode}.tsv",[best])
    print(mode,best["loss"],best["donor_margin"],best["stationarity"],flush=True)


if __name__=="__main__":
    run(sys.argv[1])
