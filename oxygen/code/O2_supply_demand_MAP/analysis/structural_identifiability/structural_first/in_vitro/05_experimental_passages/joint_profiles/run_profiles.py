"""Validate the local engineering pilot and dispatch the observed J scans."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys

from scan import DESIGN, PILOT, Run, fingerprint, independent_score, read, write
from validate import safe_interior
import numpy as np


def verify_interface(run):
    model, conf = run.model, run.conf
    def event_name(slot):
        return model.slot_ids[slot].replace("SUM159_NLS","SUM-159_NLS")
    mask = {(r["modality"],r["event"]):r for r in read(DESIGN/"event_mask.tsv")}
    checks=[]
    for assay,slots in (("corrected_count",model.indices),("death_fraction",model.classifier_slots),
                        ("karyotype",model.kary_rows)):
        ids=[event_name(s) for s in slots]
        expected={name for kind,name in mask if kind==assay}
        assert expected==set(ids) and len(ids)==len(expected)
        if assay=="death_fraction":
            assert all(float(mask[assay,e]["denominator"])==n for e,n in zip(ids,model.eligible))
            assert all(event_name(i)==event_name(j) for i,j in
                       zip(model.indices[model.death_log_index],model.classifier_slots))
        if assay=="karyotype":
            assert all(int(mask[assay,e]["sample_size"])==k for e,k in zip(ids,model.kary.sum(axis=1)))
        checks.append(dict(assay=assay,events=len(ids),exact_id_order_bound=True,passed=True))
    write(run.out/"observer_binding.tsv",checks)

    x = safe_interior(model,conf)[0]
    p,jac = conf.unmap(x)
    independent_score(model,p)
    source = fingerprint(model,p)
    analytic = conf.gradient(model.score(p,True)["gradient"],jac)
    rows=[]
    for i,name in enumerate(conf.names):
        options=[]
        for h in (1e-6,4e-6,2.5e-7,1.5e-5):
            l,r=x.copy(),x.copy()
            l[i]-=h
            r[i]+=h
            a,b=(conf.unmap(u)[0] for u in (l,r))
            same=fingerprint(model,a)==source==fingerprint(model,b)
            fd=(model.score(b)["loss"]-model.score(a)["loss"])/(2*h)
            error=abs(fd-analytic[i])/max(1.,abs(fd),abs(analytic[i]))
            options.append((not same,error,h,fd))
            if same and error<=5e-4:
                break
        changed,error,h,fd=min(options)
        rows.append(dict(parameter=name,analytic=analytic[i],finite_difference=fd,
                         step=h,scaled_error=error,transition_support_unchanged=not changed,
                         passed=not changed and error<=5e-4))
    write(run.out/"gradient_checks.tsv",rows)
    assert all(r["passed"] for r in rows)


def verify_pilot(run):
    baseline, = read(run.out/"baseline.tsv")
    assert baseline["status"]=="verified" and int(baseline["agreeing_starts"])>=2
    checks=[]
    for p in PILOT:
        rows=read(run.out/"profiles"/p/"pilot_grid.tsv") or read(run.out/"profiles"/p/"grid.tsv")
        assert len(rows)>=2 and all(r["status"]=="verified" for r in rows),p
        pairs=read(run.out/"bootstrap_pilot"/p/"pairs.tsv")
        assert pairs and max(float(r["likelihood_scale_agreement"]) for r in pairs)<2e-7,p
        valid = sum(r["pair_status"]=="verified" for r in pairs)
        assert valid >= 1, f"No fully verified paired refit for {p}"
        checks.append(dict(parameter=p,verified_measured_grid_points=len(rows),
                           paired_refits_recorded=len(pairs),verified_pairs=valid,
                           first_verified_pair_cpu_seconds=next(r["total_cpu_seconds"] for r in pairs
                                                               if r["pair_status"]=="verified"),passed=True))
    write(run.out/"pilot_release.tsv",checks)


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument("out",type=Path)
    cli.add_argument("--workers",type=int,default=4)
    cli.add_argument("--validate-only",action="store_true")
    cli.add_argument("--refine",action="store_true")
    args=cli.parse_args()
    run=Run(args.out)
    verify_interface(run)
    if args.validate_only:
        return
    verify_pilot(run)
    phase="refine" if args.refine else "profile"
    print(f"Four-parameter engineering gate passed. Dispatching {phase}.",flush=True)
    def launch(name):
        with (args.out/f"{phase}_{name}.log").open("a") as log:
            subprocess.run([sys.executable,"-B",str(Path(__file__).with_name("scan.py")),
                            phase,str(args.out),name], stdout=log,stderr=subprocess.STDOUT,check=True)
        print(f"Collected {name}",flush=True)

    with ThreadPoolExecutor(max_workers=args.workers) as workers:
        list(workers.map(launch,run.conf.targets))
    run.reconcile()


if __name__=="__main__":
    main()
