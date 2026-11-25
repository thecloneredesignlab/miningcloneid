"""Run disjoint paired refits with independent first-wave gates and bounded CPU reservations."""

import argparse
from collections import Counter
from contextlib import redirect_stdout
from math import ceil, floor
from pathlib import Path
from time import monotonic
import fcntl
import os
import resource
import subprocess
import sys
import calibrate_bootstrap as batch
from endpoint_bootstrap_release import gate, PER_ANCHOR
from endpoint_zero_review import assess
from second_endpoint_prepare import FIRST, CALIBRATION_SECONDS, SETUP_SECONDS, seed_manifest
from recover_bootstrap_pairs import BudgetEnd
from scan import ARCHIVES, cpu, digest, read, write
from seal_endpoint_release import check_bounds

SCIENCE = ("second_endpoint_prepare.py","second_endpoint_workers.py")
MARGIN = 45


def state(root,slot):
    folder=root/f"anchor_{slot}"
    obs,=read(folder/"observed_receipt.tsv")
    name=obs["parameter"]
    rows=batch.read_pairs(folder,name)
    partial=list((folder/"pairs"/name).glob("replicate_*/partial.tsv"))
    return folder,name,rows,partial


def claim(folder,name,end):
    """Claim one unused frozen draw; pending work can never breach the failure stop."""
    parent=folder/"pairs"/name
    parent.mkdir(parents=True,exist_ok=True)
    with (folder/".controller.lock").open("a") as stream:
        fcntl.flock(stream,fcntl.LOCK_EX)
        pending=0
        failures=0
        occupied=set()
        for directory in parent.glob("replicate_*"):
            j=int(directory.name[10:])
            occupied.add(j)
            pair=read(directory/"pair.tsv")
            if pair:
                failures+=pair[0]["status"]!="verified"
            elif (directory/"partial.tsv").is_file():
                failures+=1
            else:
                pending+=1
        if failures+pending>=batch.FAIL_STOP:
            return None
        for j in range(1,end+1):
            if j not in occupied:
                target=parent/f"replicate_{j:03d}"
                target.mkdir(exist_ok=False)
                return j,target
        return None


def worker(out,slot,end,seconds):
    started=cpu()
    deadline=started+seconds-MARGIN
    folder,name,_,_=state(out,slot)
    batch.check_seal(out/"RUN.sha256")
    batch.check_seal(folder/"FREEZE.sha256")
    for r in read(out/"runtime_source.tsv"):
        assert digest(Path(__file__).with_name(r["module"]))==r["sha256"]
    frozen=read(folder/"seed_manifest.tsv")
    assert frozen==[{k:str(v) for k,v in s.items()} for s in seed_manifest(name,slot)]
    if end>25:
        g,=read(folder/"wave_gate.tsv")
        assert g["release_remaining"]=="True"
    fit=batch.PairedRefits(folder,deadline,name)
    while cpu()<deadline-10:
        item=claim(folder,name,end)
        if item is None:
            break
        j,target=item
        before=cpu()
        try:
            fit.pair(frozen[j-1],target)
            r,=read(target/"pair.tsv")
            assert r["seed"]==frozen[j-1]["seed"]
            batch.pair_receipt(folder,r)
        except BudgetEnd:
            write(target/"partial.tsv",[dict(replicate=j,seed=frozen[j-1]["seed"],reason="CPU reservation reached",
                                        cpu_seconds=cpu()-before)])
            break


def charged(out):
    return sum(float(r["cpu_seconds"]) for r in read(out/"worker_cpu.tsv"))


def spent_at(out,s):
    return sum(float(r["cpu_seconds"]) for r in read(out/"worker_cpu.tsv") if int(r["slot"])==s)


def run_group(out,stage,plans):
    # One process handles each work item. The coordinator handles the last item itself.
    assert 1<=len(plans)<=4
    ledger=read(out/"worker_cpu.tsv")
    assert charged(out)+sum(t for _,_,t in plans)<=CALIBRATION_SECONDS-10
    for s in {s for s,_,_ in plans}:
        assert spent_at(out,s)+sum(t for k,_,t in plans if k==s)<=PER_ANCHOR
    rfile=out/"dispatch"/(stage+".tsv")
    assert not rfile.exists()
    write(rfile,[dict(slot=s,replicate_limit=j,cpu_reservation=t) for s,j,t in plans])
    children=[]
    for i,(s,j,seconds) in enumerate(plans[:-1]):
        with (out/"dispatch"/f"{stage}_{i}.log").open("x") as log:
            p=subprocess.Popen([sys.executable,"-Bu",__file__,"worker",str(out),str(s),str(j),str(seconds)],
                               stdout=log,stderr=subprocess.STDOUT)
        children.append((p,i,s,seconds))
    s,j,seconds=plans[-1]
    t0=cpu()
    soft,hard=resource.getrlimit(resource.RLIMIT_CPU)
    status=0
    try:
        resource.setrlimit(resource.RLIMIT_CPU,(ceil(t0)+seconds,hard))
        with (out/"dispatch"/f"{stage}_{len(plans)-1}.log").open("x") as log, redirect_stdout(log):
            worker(out,s,j,seconds)
    except Exception:
        status=1
        raise
    finally:
        resource.setrlimit(resource.RLIMIT_CPU,(soft,hard))
        ledger.append(dict(stage=stage,worker=len(plans)-1,slot=s,reserved=seconds,cpu_seconds=cpu()-t0,
                           exit_code=status,role="Coordinator worker"))
        for proc,i,s,seconds in children:
            _,exit_status,usage=os.wait4(proc.pid,0)
            proc.returncode=os.waitstatus_to_exitcode(exit_status)
            ledger.append(dict(stage=stage,worker=i,slot=s,reserved=seconds,cpu_seconds=usage.ru_utime+usage.ru_stime,
                               exit_code=proc.returncode,role="Child worker"))
        write(out/"worker_cpu.tsv",ledger)
    assert all(int(r["exit_code"])==0 and float(r["cpu_seconds"])<=int(r["reserved"])+1 for r in ledger)


def freeze_runtime(out):
    assert len(read(out/"recovery_pairs.tsv"))==4 and (out/"prepare_cpu_1_2_3_4.tsv").is_file()
    batch.reproducibility_control()
    assert gate(23,25,1000,[20]*25)["release_remaining"]
    assert not gate(22,25,1000,[20]*25)["release_remaining"]
    assert not gate(25,25,1000,[200]*25)["release_remaining"]
    # Verify actual native dependencies and the archived observation bindings.
    batch.check_sources(out)
    batch.anonymize_provenance(out)
    code={r["code"]:r["sha256"] for r in read(FIRST/"code_receipt.tsv")}
    for k,h in code.items():
        assert digest(Path(__file__).with_name(k))==h
    for k in SCIENCE:
        code[k]=digest(Path(__file__).with_name(k))
    write(out/"runtime_source.tsv",[dict(module=k,sha256=h) for k,h in code.items()])
    write(out/"design.tsv",[dict(first_release_hash=digest(FIRST/"SHA256SUMS"),primary_slots=4,reserve_slots=2,
          calibration_cpu_seconds=CALIBRATION_SECONDS,overall_core_hours=25.46,per_anchor_cpu_seconds=PER_ANCHOR,
          max_processes_including_coordinator=4,threads=1,required_first_wave=25,minimum_verified=23,
          failure_stop=10,score="Strict donor inheritance, three-readout J",candidate_resolution=0.02,
          missing_midpoint_rule="One supported-side probe with the original bracket preserved",
          holdout="Flow and A19 specimen gates retained")])
    batch.write_seal(out/"RUN.sha256",[out/k for k in ("runtime_source.tsv","design.tsv","source_receipt.tsv","environment.tsv","recovery_calibration.tsv")])


def allocate(active,needed,maximum,available):
    """Reserve every forecast before dividing extra headroom among active anchors."""
    chosen=[]
    remain=available
    limits={}
    for s in active:
        if needed[s]<=min(maximum[s],remain):
            chosen.append(s)
            limits[s]=needed[s]
            remain-=needed[s]
    while remain>=len(chosen)>0:
        short=[s for s in chosen if limits[s]<maximum[s]]
        if not short:
            break
        share=max(1,remain//len(short))
        for s in short:
            inc=min(maximum[s]-limits[s],share,remain)
            limits[s]+=inc
            remain-=inc
    return chosen,limits


def plans(limits,ends):
    # Additional processes claim disjoint replicates on a shared anchor under a file lock.
    count=Counter({s:1 for s in limits})
    while count and sum(count.values())<4:
        s=max(count,key=lambda k:limits[k]/count[k])
        if limits[s]//(count[s]+1)<90: break
        count[s]+=1
    return [(s,ends[s],limits[s]//n) for s,n in count.items() for _ in range(n)]


def decide(out,slot,reason):
    folder,name,rr,partial=state(out,slot)
    anchor,=read(folder/"anchors.tsv")
    result=assess(name,anchor["fixed_value"],float(anchor["observed_statistic"]),rr)
    check_bounds({k:str(v) for k,v in result.items()},[float(r["statistic"]) for r in rr if r["status"]=="verified"])
    write(folder/"outcome.tsv",[dict(**result,partial_draws=len(partial),finish_reason=reason)])
    print(f'Slot {slot} at {anchor["fixed_value"]}: {result["decision"]}; {len(rr)}/199 processed. {reason}',flush=True)


def calibrate(out,slots,label):
    valid=[s for s in slots if (out/f"anchor_{s}"/"anchors.tsv").is_file()]
    if not valid:
        return
    available=floor(CALIBRATION_SECONDS-charged(out))-20
    limits={s:min(PER_ANCHOR-10,available//len(valid)) for s in valid}
    if min(limits.values())<90:
        for s in valid: decide(out,s,"Shared budget exhausted before first wave")
        return
    run_group(out,label+"_wave",plans(limits,{s:25 for s in valid}))
    active=[]
    for s in valid:
        folder,_,rr,partial=state(out,s)
        g=gate(sum(r["status"]=="verified" for r in rr),len(rr),spent_at(out,s),
               [float(r["pair_cpu_seconds"]) for r in rr])
        write(folder/"wave_gate.tsv",[dict(**g,paired_receipts_verified=True)])
        if g["release_remaining"] and not partial:
            active.append(s)
        else:
            decide(out,s,"First-wave quality, completeness or cost gate held")
    round_number=0
    while active:
        round_number+=1
        needed={}
        maximum={}
        ends={}
        available=floor(CALIBRATION_SECONDS-charged(out))-20
        for s in list(active):
            folder,name,rr,partial=state(out,s)
            done=len(rr)+len(partial)
            failures=sum(r["status"]!="verified" for r in rr)+len(partial)
            cap=floor(PER_ANCHOR-spent_at(out,s))-10
            if done==199 or failures>=batch.FAIL_STOP or cap<=90:
                reason=("All planned replicates reached a terminal state" if done==199 else
                        "Ten unresolved pairs" if failures>=10 else "Individual CPU cap")
                decide(out,s,reason)
                active.remove(s)
                continue
            g,=read(folder/"wave_gate.tsv")
            future=ceil(float(g["planning_cpu_seconds_per_pair"])*(199-done))+30
            if done>25:
                # Preserve the first-wave release; use actual earned budget on completed pairs.
                future=min(future,ceil(1.25*spent_at(out,s)/done*(199-done))+30)
            if future>cap and done==25:
                decide(out,s,"Individual CPU forecast exceeded the registered limit")
                active.remove(s)
                continue
            needed[s]=min(future,cap)
            maximum[s]=cap
            ends[s]=min(199,done+50)
        active,queued=list(active),[]
        if not active:break
        picked,limits=allocate(active,needed,maximum,available)
        if not picked:
            for s in active:decide(out,s,"Shared budget cannot cover remaining forecast")
            break
        run_group(out,f"{label}_main_{round_number}",plans(limits,ends))


def main_stage(out,phase):
    assert out.resolve().parent==ARCHIVES.resolve()
    start=cpu()
    base_local=sum(float(r["cpu_seconds"]) for r in read(out/"worker_cpu.tsv") if r["role"]=="Coordinator worker")
    with (out/".controller.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if phase=="primary":
            assert not (out/"RUN.sha256").exists()
            freeze_runtime(out)
            calibrate(out,[1,2,3,4],"primary")
        else:
            batch.check_seal(out/"RUN.sha256")
            slots=[]
            for s,parent in ((5,2),(6,1)):
                folder=out/f"anchor_{parent}"
                rows=read(folder/"outcome.tsv")
                if not rows or int(rows[0]["processed"])!=199:
                    continue
                slots.append(s)
            reserve=[]
            # Gate entry on the remaining shared account using the completed primary costs.
            available=CALIBRATION_SECONDS-charged(out)
            for s in slots:
                if s==5: donor=2
                elif read(out/"anchor_1"/"outcome.tsv")[0]["decision"]=="Inconclusive": donor=1
                else: donor=None
                forecast=(1.25*spent_at(out,donor) if donor else
                     1.25*float(next(r["worker_cpu_seconds"] for r in read(FIRST/"worker_resources.tsv") if r["parameter"]=="p_mis_base")))
                if min(PER_ANCHOR,ceil(forecast)+60)<=available:
                    reserve.append(s)
                    available-=min(PER_ANCHOR,ceil(forecast)+60)
            write(out/"reserve_release.tsv",[dict(slot=s,enabled=s in reserve,
                condition="Parent completed and CPU forecast fits" if s in reserve else "Parent or shared budget held")
                for s in (5,6)])
            if reserve:
                prep_used=float(read(out/"recovery_cpu.tsv")[0]["cpu_seconds"])+sum(float(read(p)[0]["cpu_seconds"]) for p in out.glob("prepare_cpu_*.tsv"))
                prep_limit=floor(SETUP_SECONDS-prep_used)-40
                assert prep_limit>30
                subprocess.run([sys.executable,"-Bu",str(Path(__file__).with_name("second_endpoint_prepare.py")),
                   "prepare",str(out),*map(str,reserve),"--cpu-seconds",str(prep_limit)],check=True)
                calibrate(out,reserve,"reserve")
        end_local=sum(float(r["cpu_seconds"]) for r in read(out/"worker_cpu.tsv") if r["role"]=="Coordinator worker")
        write(out/f"{phase}_controller_cpu.tsv",[dict(cpu_seconds=max(0,cpu()-(end_local-base_local)))])
        print(f"{phase.capitalize()} dispatch finished.",flush=True)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("phase",choices=("primary","reserve","worker"))
    cli.add_argument("out",type=Path)
    cli.add_argument("numbers",type=int,nargs="*")
    args=cli.parse_args()
    if args.phase=="worker":
        slot,end,seconds=args.numbers
        batch.install_cpu_cap(seconds)
        worker(args.out,slot,end,seconds)
    else:
        main_stage(args.out,args.phase)
