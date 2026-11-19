"""Resolve tolerance decisions around the strongest real-data AR(1) reference."""
import argparse
from concurrent.futures import as_completed

import numpy as np

from ar1_model import (HERE, STATIONARY, THRESHOLD, Config, choose,
                       likelihood_ratio, model, pool, read, receipt, row_x, solve, write)
from profile_ar1 import append_labels, envelope, load_profiles, old_bank, reference, update_reference


def task(job):
    name,q,start_points,scales=job
    k=model()
    c=Config(k)
    j=c.names.index(name)
    records=[]
    for start,d in zip(start_points,scales):
        f=solve(k,c,k.y,start,limit=91,divisor=d,fixed={j:q})
        r=receipt(k,c,k.y,f)
        records.append(r)
        qual=sorted((r for r in records if float(r["donor_margin"])>=0 and
                      float(r["stationarity"])<=STATIONARY),key=lambda r:float(r["loss"]))
        if len(qual)>=2 and likelihood_ratio(float(qual[1]["loss"]),float(qual[0]["loss"]),221)<=.06:
            break
    return append_labels(choose(records,k.n),name,q,"profile_edge_refit")


def run(workers,version):
    k=model()
    c=Config(k)
    f=reference()
    center=row_x(c,f)
    profiles=read(HERE/"profile_intervals.tsv")
    direct=load_profiles()
    bank=old_bank(c)
    global_bank=[row_x(c,r) for r in read(HERE/"observed_refits.tsv") if float(r["donor_margin"])>=0]
    jobs=[]
    for r in profiles:
        if r["profile_status"]=="outside_supported":
            continue
        name=r["parameter"]
        j=c.names.index(name)
        all_name=[v for v in direct if v["parameter"]==name and float(v["donor_margin"])>=0]
        chosen=[]
        for side,start_type in ((-1,"lower"),(1,"upper")):
            bound=center[j]+side*.10001
            if 0.<bound<1.:
                tests=[t for t in all_name if abs(float(t["coordinate"])-bound)<.003 and
                       float(t["stationarity"])<=STATIONARY and
                       likelihood_ratio(float(t["loss"]),float(f["loss"]),221)>THRESHOLD]
                if not tests:
                    chosen.extend((bound,center[j]+side*.087))
            endpoint=r[f"{start_type}_normalized"]
            if endpoint and 1e-6<float(endpoint)<1-1e-6 and (
                    r["profile_status"]=="unresolved" or r[f"{start_type}_type"]=="approximate_crossing"):
                chosen.append(float(endpoint))
        for value in sorted(set(chosen)):
            near_native=sorted(all_name,key=lambda r:(abs(float(r["coordinate"])-value),
                                                      float(r["loss"])))[:5]
            near_bank=bank[np.argsort(abs(bank[:,j]-value))[:25]]
            options=[center,*[row_x(c,r) for r in near_native],*global_bank,*near_bank]
            evaluation=[]
            for z in options:
                x=z.copy()
                x[j]=value
                out,_=k.metrics(c.unmap(x)[0],k.y)
                evaluation.append((out[1]<0,out[0] if out[1]>=0 else -out[1],x))
            evaluation.sort(key=lambda s:s[:2])
            initials=[]
            for _,_,q in evaluation:
                if all(np.linalg.norm(q-v)>.009 for v in initials):
                    initials.append(q)
                if len(initials)>=4:
                    break
            if not initials:
                initials=[center]
            if len(initials)==1:
                initials.append(center)
            jobs.append((name,value,initials,(18.,5.,35.,12.)))
    if jobs:
        path=HERE/f"profile_refined_{version}.tsv"
        rows=[]
        with pool(workers) as p:
            for f in as_completed([p.submit(task,j) for j in jobs]):
                v=f.result()
                rows.append(v)
                print("edge",v["parameter"],v["coordinate"],v["stationarity"],v["loss"],flush=True)
                write(path,rows)
        update_reference(workers)
    for r in envelope():
        print(r["parameter"],r["profile_status"],r["lower_normalized"],r["upper_normalized"],flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--workers",type=int,default=8)
    p.add_argument("--version",default="1")
    a=p.parse_args()
    run(a.workers,a.version)
