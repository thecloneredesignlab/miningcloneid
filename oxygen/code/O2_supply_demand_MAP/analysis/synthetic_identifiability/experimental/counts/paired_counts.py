"""Replay the frozen OLS simulations and fit the paired direct-count samples."""
import argparse
from collections import defaultdict
import concurrent.futures as cf
import zlib

import numpy as np

from count_model import HERE, OLS, model, optimize
from inference import likelihood_ratio, pool
from practical import Configuration, read, result_row, row_x, write
from recovery import CONVERGE_LR, KKT_TOL, PILOT_REQUIRED, SUPPLY_TOL


def replay():
    k=model()
    c=Configuration(k)
    truth={r["scenario"]:np.array([float(r[q]) for q in c.names]) for r in read(OLS/"synthetic_truth.tsv")}
    sigma=float(read(OLS/"truth_design.tsv")[0]["sigma_log_count"])
    old=defaultdict(dict)
    for r in read(OLS/"noisy_slope_inputs.tsv"):
        old[r["scenario"],int(r["dataset"])][r["passage_id"]]=float(r["ols_slope"])
    bases={}
    for name, x in truth.items():
        logs,diag=k.log_predict(c.unmap(x)[0])
        assert diag[1]>=-SUPPLY_TOL
        bases[name]=logs
    rng_inputs=[]
    checks=[]
    for i in range(1,101):
        for name,logs in bases.items():
            rng=np.random.default_rng(np.random.SeedSequence(
                [20261118,1 if name.startswith("A") else 2,i]))
            simulated=logs.copy()
            simulated[k.capture] += rng.normal(0,sigma,k.n)
            simulated = np.maximum(0, simulated)
            reconstructed=k.slopes(simulated)
            reference=np.array([old[name,i][p] for p in k.ids])
            error=max(abs(reference-reconstructed))
            assert error <= 1e-12, (name,i,error)
            checks.append(dict(scenario=name,dataset=i,seed=f"20261118:{1 if name.startswith('A') else 2}:{i}",
                               assumed_sigma=sigma,ols_slope_max_abs_discrepancy=error))
            rng_inputs.extend(dict(scenario=name,dataset=i,record_id=k.slots[j]["record_id"],log_count=simulated[j])
                              for j in k.indices)
    write(HERE/"paired_replay.tsv",checks)
    write(HERE/"paired_log_counts.tsv",rng_inputs)
    print("Paired replay",len(checks),"datasets; maximum OLS discrepancy",max(r["ols_slope_max_abs_discrepancy"] for r in checks),flush=True)


def loads():
    k=model()
    records=defaultdict(list)
    observed=defaultdict(list)
    for r in read(HERE/"paired_log_counts.tsv"):
        key=r["scenario"],int(r["dataset"])
        records[key].append(r["record_id"])
        observed[key].append(float(r["log_count"]))
    expected=[k.slots[j]["record_id"] for j in k.indices]
    assert len(records)==200 and all(names==expected for names in records.values())
    return {key:np.array(y) for key,y in observed.items()}


def judge(k,fits,incumbent=None):
    feasible=[f for f in fits if f["supply"]>=-SUPPLY_TOL and f["loss"]<1e9]
    if not feasible:
        return min(fits,key=lambda f:f["loss"]),False,np.inf
    floor=min(f["loss"] for f in feasible)
    if incumbent and incumbent["supply"]>=-SUPPLY_TOL:
        floor=min(floor,incumbent["loss"])
    stable=[f for f in feasible if f["gradient"]<=KKT_TOL]
    if incumbent and incumbent["supply"]>=-SUPPLY_TOL and incumbent["gradient"]<=KKT_TOL:
        stable.append(incumbent)
    stable.sort(key=lambda f:f["loss"])
    best=min(feasible,key=lambda f:f["loss"])
    if incumbent and incumbent["supply"]>=-SUPPLY_TOL and incumbent["loss"]<best["loss"]:
        best=incumbent
    if not stable:
        return best,False,np.inf
    gap=likelihood_ratio(stable[1]["loss"],stable[0]["loss"],k.n) if len(stable)>1 else np.inf
    if gap<=CONVERGE_LR and likelihood_ratio(stable[0]["loss"],floor,k.n)<=.001:
        return stable[0],True,gap
    return best,False,gap


def old_fit(c,r):
    if r is None:
        return None
    return dict(x=row_x(c,r), loss=float(r["loss"]), supply=float(r["min_aggregate_donor_margin_fraction"]),
                gradient=float(r["max_projected_gradient"]),calls=int(r["calls"]),
                iteration=int(r["iterations"]),scipy_success=r["scipy_success"]=="True",
                termination=r["termination"],seconds=float(r["seconds"]))


def worker(task):
    case,idx,y,paired,reference,last,alternate,phase=task
    k=model()
    c=Configuration(k)
    seed=zlib.crc32(f"{case}:{idx}:{phase}".encode())+67019
    rng=np.random.default_rng(seed)
    jolt=lambda x,s:np.clip(x+rng.normal(0,s,len(x)),0.,1.)
    incumbent=old_fit(c,last)
    if incumbent is None:
        starts=[(paired,10.,51),(jolt(reference,.004),8.,52)]
    else:
        base=incumbent["x"]
        starts=[(jolt(base,.00003),25.,70),(jolt(base,.0004),8.,71),
                (jolt(base,.004),2.,75),(jolt(alternate,.001),10.,76),(jolt(reference,.015),10.,80)]
    fits=[]
    got=False
    for start,scale,n in starts:
        fits.append(optimize(k,c,y,start,n,divisor=scale))
        best,got,gap=judge(k,fits,incumbent)
        if got:
            break
    if not got and incumbent is None:
        a=best["x"]
        for s,d in ((.00025,20.),(.0012,3.)):
            fits.append(optimize(k,c,y,jolt(a,s),57,divisor=d))
            best,got,gap=judge(k,fits)
            if got:
                break
    return result_row(c,best,scenario=case,dataset=idx,qualified=got,independent_LR_gap=gap,
                      attempts=len(fits),optimization_phase=phase)


def combine():
    groups=defaultdict(list)
    for f in sorted(HERE.glob("count_recovery_*.tsv")):
        if f.stem not in ("count_recovery_pilot","count_recovery_expanded") and not f.stem.startswith("count_recovery_repair_"):
            continue
        for r in read(f):
            groups[r["scenario"],int(r["dataset"])].append(r)
    result=[]
    for _,rows in sorted(groups.items()):
        candidates=[r for r in rows if float(r["min_aggregate_donor_margin_fraction"])>=-SUPPLY_TOL]
        floor=min(candidates or rows,key=lambda r:float(r["loss"]))
        confirmed=[r for r in candidates if r["qualified"]=="True"
                   and likelihood_ratio(float(r["loss"]),float(floor["loss"]),221)<=CONVERGE_LR]
        result.append(min(confirmed,key=lambda r:float(r["loss"])) if confirmed else floor)
    if result:
        write(HERE/"qualified_count_recovery.tsv",result)
    return result


def pilot_gate(rows, release=False):
    gate=[]
    for name in ("A_real_fit","B_weak_direction"):
        sample=[r for r in rows if r["scenario"]==name and int(r["dataset"])<=20]
        n=sum(r["qualified"]=="True" for r in sample)
        gate.append(dict(scenario=name,assigned=len(sample),qualified=n,threshold=PILOT_REQUIRED,
                         released=len(sample)==20 and n>=PILOT_REQUIRED))
    write(HERE/"pilot_gate.tsv",gate)
    ready=all(g["released"] for g in gate)
    path=HERE/"pilot_release.tsv"
    if ready and release and not path.exists():
        write(path,gate)
    print("Pilot gate",[(g["scenario"],g["assigned"],g["qualified"]) for g in gate],ready,flush=True)
    return ready


def batch(phase,workers,count,version):
    original=loads()
    k=model()
    c=Configuration(k)
    old={(r["scenario"],int(r["dataset"])):r for r in read(OLS/"qualified_recovery.tsv")}
    current=combine()
    existing={(r["scenario"],int(r["dataset"])):r for r in current}
    reference=row_x(c,read(HERE/"real_best.tsv")[0])
    if phase=="expanded":
        if not pilot_gate(current,release=True):
            raise RuntimeError("Count pilot must pass before expansion")
    filename=f"count_recovery_repair_{version}.tsv" if phase=="repair" else f"count_recovery_{phase}.tsv"
    output=HERE/filename
    results=read(output) if output.exists() else []
    finished={(r["scenario"],int(r["dataset"])) for r in results}
    jobs=[]
    high_quality=[r for r in current if r["qualified"]=="True"]
    for key,data in original.items():
        case,i=key
        if i>count or key in finished:
            continue
        last=existing.get(key)
        if phase=="repair":
            if last is None or last["qualified"]=="True":
                continue
        elif (phase=="pilot" and i>20) or (phase=="expanded" and i<=20):
            continue
        source=row_x(c,old[key])
        alternate=reference
        if phase=="repair" and high_quality:
            closest=min(high_quality,key=lambda q:np.linalg.norm(
                original[(q["scenario"],int(q["dataset"]))]-data))
            alternate=row_x(c,closest)
        jobs.append((case,i,data,source,reference,last,alternate,f"{phase}:{version}"))
    with pool(workers) as work:
        for f in cf.as_completed([work.submit(worker,t) for t in jobs]):
            r=f.result()
            results.append(r)
            write(output,results)
            print(r["scenario"],r["dataset"],r["qualified"],r["loss"],
                  r["max_projected_gradient"],r["independent_LR_gap"],flush=True)
    pilot_gate(combine())


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("phase",choices=("replay","pilot","expanded","repair","combine"))
    parser.add_argument("--workers",type=int,default=9)
    parser.add_argument("--max",type=int,choices=(20,100),default=100)
    parser.add_argument("--version",default="1")
    args=parser.parse_args()
    if args.phase=="replay":
        replay()
    elif args.phase=="combine":
        pilot_gate(combine())
    else:
        batch(args.phase,args.workers,args.max,args.version)
