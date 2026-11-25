"""Advance sealed three-readout endpoint batches locally and publish their V6 panel."""

import argparse
import fcntl
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
import subprocess
import sys
import time
from zoneinfo import ZoneInfo

from scan import ARCHIVES, Run, cpu, digest, read, write
import calibrate_bootstrap as bootstrap
from endpoint_zero_review import PREVIOUS, SOURCE
from seal_endpoint_release import check_bounds
from v6_anchor import calibrate, prepare, seeds
from v6_panel import DEST, render
from v6_sides import FIRST, SECOND, PRIORITY, RESOLUTION, inventory, queue, source_points, split_candidate

HERE = Path(__file__).resolve().parent
OUT = ARCHIVES/"20261124_joint_support_v6_01"


def now():
    return datetime.now(ZoneInfo("America/New_York")).isoformat(timespec="seconds")


def status(stage,detail):
    row=dict(time_new_york=now(),stage=stage,detail=detail)
    write(OUT/"progress.tsv",[row])
    print(f'{row["time_new_york"]}  {stage}  {detail}',flush=True)


def dependencies():
    modules = {r["module"]:r["sha256"] for r in read(SECOND/"runtime_source.tsv")}
    for label,expected in modules.items():
        assert digest(HERE/label)==expected,label
    for label in ("v6_anchor.py","v6_endpoints.py","v6_sides.py","seal_endpoint_release.py"):
        modules[label]=digest(HERE/label)
    return [dict(module=k,sha256=v) for k,v in sorted(modules.items())]


def wait_for_parent():
    status("Awaiting parent seal","No new observed-data or bootstrap fits before the second release passes hash verification")
    while not (SECOND/"SHA256SUMS").is_file():
        done=SECOND/"reserve_controller_cpu.tsv"
        if done.is_file() and time.time()-done.stat().st_mtime>180:
            # Allow the predecessor's incumbent collector to finish before assuming
            # responsibility for the saved receipts. Do not dispatch prior fits.
            audit=SECOND/"final_audit_source.tsv"
            if not audit.is_file() or time.time()-audit.stat().st_mtime>180:
                cmd=[sys.executable,"-Bu",str(HERE/"seal_second_endpoints.py"),str(SECOND)]
                result=subprocess.run(cmd,check=False)
                if result.returncode and not (SECOND/"SHA256SUMS").is_file():
                    raise RuntimeError("Parent sealing failed. Inspect the saved parent receipts.")
        if not (SECOND/"SHA256SUMS").is_file():
            time.sleep(55)
    for path in (FIRST,SECOND,SOURCE,PREVIOUS):
        bootstrap.check_seal(path/"SHA256SUMS")
    slots=read(SECOND/"anchor_audit.tsv")
    assert set(range(1,5))<={int(r["slot"]) for r in slots}
    assert len(read(SECOND/"endpoint_brackets.tsv"))==28


def verify_setup():
    receipts=dependencies()
    if (OUT/"science_freeze.tsv").exists():
        assert receipts==read(OUT/"science_freeze.tsv"),"Runtime drift since V6 was frozen"
        bootstrap.check_seal(OUT/"INPUTS.sha256")
    else:
        bootstrap.reproducibility_control()
        bootstrap.check_sources(OUT)
        bootstrap.anonymize_provenance(OUT)
        write(OUT/"science_freeze.tsv",receipts)
        write(OUT/"parents.tsv",[dict(packet=p.name,sha256=digest(p/"SHA256SUMS")) for p in (SOURCE,PREVIOUS,FIRST,SECOND)])
        write(OUT/"design.tsv",[dict(contract="Strict donor inheritance, three-readout J",parameters=14,sides=28,
             resolution=RESOLUTION,pairs=199,first_wave=25,minimum_verified=23,unresolved_pair_stop=10,
             maximum_local_scientific_processes=4,numerical_threads_each=1,calibration_core_hours_per_anchor=7,
             maximum_calibration_core_hours_per_batch=28,scope="Sampled component containing the fitted reference",
             precision_stop="Undecided midpoint and both flanks or registered search resolution exhausted")])
        bootstrap.write_seal(OUT/"INPUTS.sha256",[OUT/f for f in ("science_freeze.tsv","parents.tsv","source_receipt.tsv","environment.tsv","design.tsv") if (OUT/f).is_file()])
    for r in read(OUT/"parents.tsv"):
        assert digest(ARCHIVES/r["packet"]/"SHA256SUMS")==r["sha256"]


def finish_slot(folder,spec,started,prior):
    checks=read(folder/"calibration.tsv")
    if checks:
        c,=checks
        pairs=bootstrap.read_pairs(folder,spec["parameter"])
        for r in pairs:
            bootstrap.pair_receipt(folder,r)
        check_bounds(c,[float(r["statistic"]) for r in pairs if r["status"]=="verified"])
        decision,reason,paired=c["disposition"],c["finish_reason"],float(c["paired_cpu_seconds"])
    else:
        decision,reason,paired="Inconclusive","Observed restricted fit failed agreement",0.
    write(folder/"result.tsv",[dict(parameter=spec["parameter"],value=spec["value"],coordinate=spec["coordinate"],
          decision=decision,reason=reason,source="Paired bootstrap" if checks else "Observed fit review",
          packet=str(folder.relative_to(ARCHIVES)),generation="V6",setup_cpu_seconds=prior,
          paired_cpu_seconds=paired,verification_cpu_seconds=max(0.,cpu()-started-prior-paired))])
    bootstrap.seal_results(folder)


def work(folder):
    with (folder/".controller.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (folder/"SHA256SUMS").is_file():
            bootstrap.check_seal(folder/"SHA256SUMS")
            return
        for r in read(OUT/"science_freeze.tsv"):
            assert digest(HERE/r["module"])==r["sha256"],r["module"]
        bootstrap.check_seal(folder.parent/"BATCH.sha256")
        spec,=read(folder/"design.tsv")
        started=cpu()
        if not (folder/"FREEZE.sha256").is_file():
            prepare(folder,spec,read(folder.parent/"all_seeds.tsv")[(int(spec["slot"])-1)*199:int(spec["slot"])*199],
                    read(folder/"neighbor_vectors.tsv"))
        bootstrap.check_seal(folder/"FREEZE.sha256")
        prep=cpu()-started
        if (folder/"anchors.tsv").exists() and not (folder/"calibration.tsv").exists():
            calibrate(folder)
        finish_slot(folder,spec,started,prep)


def saved_points():
    previous=[]
    for path in sorted(OUT.glob("batch_*/slot_*/result.tsv")):
        if not (path.parent/"SHA256SUMS").exists():
            continue
        record,=read(path)
        previous.append(record)
    return previous


def neighbor_vectors(name,previous,conf):
    rows=[]
    for folder in [*sorted(SECOND.glob("anchor_*")),*[ARCHIVES/r["packet"] for r in previous if r["parameter"]==name]]:
        vectors=folder/"generator_vectors.tsv"
        anchor=read(folder/"anchors.tsv")
        if vectors.is_file() and anchor and anchor[0]["parameter"]==name:
            rows.append(dict(coordinate=conf.target_coordinate(name,float(anchor[0]["fixed_value"])),
                             vector=str(vectors.relative_to(ARCHIVES.parents[2]))))
    return rows


def freeze_batch(folder,selected,previous,run):
    folder.mkdir(exist_ok=True)
    number=int(folder.name.split("_")[1])
    frozen=[]
    files=[]
    for slot,row in enumerate(selected,1):
        name=row["parameter"]
        place=folder/f"slot_{slot}"
        place.mkdir(exist_ok=True)
        write(place/"design.tsv",[dict(batch=number,slot=slot,parameter=name,side=row["side"],value=row["candidate"],
              coordinate=run.conf.target_coordinate(name,float(row["candidate"])),supported_at_dispatch=row["supported"],
              excluded_at_dispatch=row["excluded"],reason="Registered endpoint search")])
        write(place/"neighbor_vectors.tsv",neighbor_vectors(name,previous,run.conf))
        files.extend(place.glob("*.tsv"))
        frozen.extend(seeds(name,number,slot))
    write(folder/"all_seeds.tsv",frozen)
    all_ids=set()
    all_states=set()
    other_files=[FIRST/"seed_manifest.tsv",PREVIOUS/"seed_manifest.tsv",SOURCE/"bootstrap_pilot_pairs.tsv",
                 *SECOND.glob("anchor_*/seed_manifest.tsv"),*OUT.glob("batch_*/all_seeds.tsv")]
    for f in other_files:
        for r in read(f):
            seed=int(r["seed"])
            assert seed not in all_ids,"Duplicate seed across releases"
            all_ids.add(seed)
            h=r.get("state_sha256")
            if h:
                assert h not in all_states,"Duplicate generator entropy across releases"
                all_states.add(h)
    bootstrap.write_seal(folder/"BATCH.sha256",files+[folder/"all_seeds.tsv"])


def dispatch(folder):
    slots=[p for p in sorted(folder.glob("slot_*")) if not (p/"SHA256SUMS").is_file()]
    assert len(slots)<=4
    children=[]
    try:
        for p in slots[:-1]:
            with (folder/(p.name+".log")).open("a") as log:
                children.append(subprocess.Popen([sys.executable,"-Bu",str(Path(__file__).resolve()),"worker",str(p)],
                                                 stdout=log,stderr=subprocess.STDOUT))
        if slots:
            with (folder/(slots[-1].name+".log")).open("a") as log, redirect_stdout(log):
                work(slots[-1])
    finally:
        codes=[p.wait() for p in children]
    assert not any(codes),f"Child worker exit status {codes}"
    checks=[]
    for p in sorted(folder.glob("slot_*")):
        bootstrap.check_seal(p/"SHA256SUMS")
        checks.append(dict(slot=p.name,sha256=digest(p/"SHA256SUMS"),**read(p/"result.tsv")[0]))
    assert sum(float(r["paired_cpu_seconds"]) for r in checks)<=28*3600
    write(folder/"anchors_completed.tsv",checks)
    bootstrap.write_seal(folder/"CLOSED.sha256",[folder/"BATCH.sha256",folder/"anchors_completed.tsv"]+
                         [p/"SHA256SUMS" for p in sorted(folder.glob("slot_*"))])


def run_all():
    OUT.mkdir(exist_ok=True)
    with (OUT/".controller.lock").open("a") as lock:
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print("The local V6 controller already owns this release.",flush=True)
            return
        if (OUT/"FINAL.sha256").is_file():
            bootstrap.check_seal(OUT/"FINAL.sha256")
            return
        wait_for_parent()
        verify_setup()
        run=Run(SOURCE)
        while True:
            folders=sorted(OUT.glob("batch_*"))
            if folders and not (folders[-1]/"CLOSED.sha256").is_file() and (folders[-1]/"BATCH.sha256").is_file():
                status("Calibrating",folders[-1].name)
                dispatch(folders[-1])
                continue
            previous=saved_points()
            points=source_points(run.conf,previous)
            sides=inventory(run,points)
            selected=queue(sides)
            write(OUT/"endpoint_status.tsv",sides)
            fields=("parameter","value","coordinate","decision","source","packet","reason")
            write(OUT/"tested_values.tsv",[{k:r[k] for k in fields} for r in points])
            render(OUT,sides,points,run.conf,final=not selected)
            if not selected:
                costs=[r for folder in folders for r in read(folder/"anchors_completed.tsv")]
                write(OUT/"resources.tsv",[dict(batches=len(folders),anchors=len(costs),
                      paired_core_hours=sum(float(r["paired_cpu_seconds"]) for r in costs)/3600,
                      preparation_core_hours=sum(float(r["setup_cpu_seconds"]) for r in costs)/3600,
                      verification_core_hours=sum(float(r["verification_cpu_seconds"]) for r in costs)/3600,
                      html_sha256=digest(DEST))])
                bootstrap.write_seal(OUT/"FINAL.sha256",[OUT/f for f in ("INPUTS.sha256","endpoint_status.tsv",
                      "tested_values.tsv","panel_a_data.tsv","resources.tsv")]+[f/"CLOSED.sha256" for f in folders])
                status("Complete","All 28 sides reached a registered endpoint, domain contact, or recorded stopping rule")
                return
            if folders and not (folders[-1]/"BATCH.sha256").exists():
                folder=folders[-1]
            else:
                folder=OUT/f"batch_{len(folders)+1:03d}"
            freeze_batch(folder,selected,previous,run)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("action",choices=("run","worker","test"))
    cli.add_argument("slot",nargs="?",type=Path)
    args=cli.parse_args()
    if args.action=="test":
        assert len(set(PRIORITY))==14 and RESOLUTION==.02
        assert split_candidate(.3,.5,[])==.4
        assert split_candidate(.3,.309,[]) is None
        assert abs(split_candidate(.3,.5,[.37])-.435)<1e-12
        for i in (1,2,3,4):
            test=seeds("test",1,i)
            assert len(test)==len({r["state_sha256"] for r in test})==199
        bootstrap.reproducibility_control()
        print("Passed endpoint selection, seed and conservative-tail controls.")
    elif args.action=="worker":
        work(args.slot)
    else:
        try:
            run_all()
        except Exception as exc:
            status("Review required",str(exc))
            raise
