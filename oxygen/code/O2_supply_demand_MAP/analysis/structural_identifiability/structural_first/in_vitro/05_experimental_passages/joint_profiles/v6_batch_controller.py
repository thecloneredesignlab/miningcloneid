"""Continue the frozen V6 run with a manifest-based slot inventory at collection."""

from contextlib import redirect_stdout
import fcntl
from pathlib import Path
import subprocess
import sys

import v6_endpoints as run


def slot_directories(folder):
    # Each slot's frozen seed range, rather than its log filename, identifies a job.
    run.bootstrap.check_seal(folder/"BATCH.sha256")
    ids=sorted({int(r["anchor_code"]) for r in run.read(folder/"all_seeds.tsv")})
    assert 1<=len(ids)<=4 and ids==list(range(1,len(ids)+1))
    slots=[folder/f"slot_{i}" for i in ids]
    for i,p in zip(ids,slots):
        design,=run.read(p/"design.tsv")
        assert int(design["slot"])==i
    return slots


def dispatch(folder):
    slots=slot_directories(folder)
    pending=[p for p in slots if not (p/"SHA256SUMS").is_file()]
    children=[]
    for p in pending:
        if (p/"pairs").exists() and not (p/"calibration.tsv").is_file():
            raise RuntimeError(f"Reconcile incomplete CPU and paired-fit receipts in {p.name} before resuming.")
    try:
        for p in pending[:-1]:
            with (folder/(p.name+".log")).open("a") as log:
                children.append(subprocess.Popen([sys.executable,"-Bu",str(run.HERE/"v6_endpoints.py"),
                                "worker",str(p)],stdout=log,stderr=subprocess.STDOUT))
        if pending:
            with (folder/(pending[-1].name+".log")).open("a") as log,redirect_stdout(log):
                run.work(pending[-1])
    finally:
        codes=[child.wait() for child in children]
    assert all(code==0 for code in codes),codes
    checks=[]
    for p in slots:
        run.bootstrap.check_seal(p/"SHA256SUMS")
        result,=run.read(p/"result.tsv")
        checks.append(dict(slot=p.name,sha256=run.digest(p/"SHA256SUMS"),**result))
    assert sum(float(r["paired_cpu_seconds"]) for r in checks)<=28*3600
    run.write(folder/"anchors_completed.tsv",checks)
    run.bootstrap.write_seal(folder/"CLOSED.sha256",[folder/"BATCH.sha256",folder/"anchors_completed.tsv",
                             *[p/"SHA256SUMS" for p in slots]])
    run.bootstrap.check_seal(folder/"CLOSED.sha256")


def resume():
    # Wait for the active owner and leave all of its workers and log handles alone.
    with (run.OUT/".controller.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (run.OUT/"FINAL.sha256").is_file():
            run.bootstrap.check_seal(run.OUT/"FINAL.sha256")
            return
        first=run.OUT/"batch_001"
        run.bootstrap.check_seal(run.OUT/"INPUTS.sha256")
        if not (first/"CLOSED.sha256").is_file():
            for p in slot_directories(first):
                with (p/".controller.lock").open("a") as slot_lock:
                    fcntl.flock(slot_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    run.bootstrap.check_seal(p/"SHA256SUMS")
        receipt=run.OUT/"collection_controller.tsv"
        record=dict(original_controller_sha256=run.digest(run.HERE/"v6_endpoints.py"),
                    collection_controller_sha256=run.digest(Path(__file__)),
                    change="Collect only slots registered in the frozen batch seed manifest",
                    worker="Frozen v6_endpoints.py worker",first_batch="batch_001")
        if receipt.is_file():
            assert run.read(receipt)==[record]
        else:
            run.write(receipt,[record])
    # run_all acquires the same exclusive owner lock and independently rechecks the
    # frozen scientific runtime. All new numerical calls still use its original work().
    run.dispatch=dispatch
    try:
        run.run_all()
    except Exception as exc:
        run.status("Review required",str(exc))
        raise


if __name__=="__main__":
    resume()
