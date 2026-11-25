"""Freeze and replay the exact strict-inheritance J point used for Panel A."""
import argparse
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from math import exp, isclose, log
import os
from pathlib import Path
import sys
from time import process_time

from paths import FIT, FROZEN_DESIGN, PROFILES, ROOT, SYNTHETIC, output_directory

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
os.environ["JOINT_KERNEL"] = "1"
sys.path.insert(0, str(PROFILES))
from profile_model import ProfileModel, ProfileConfig
from runtime import read, write
import numpy as np

ALL = SYNTHETIC / "20261123_experimental_joint_all_mechanism_profiles_01"
HTML = ROOT / "oxygen/results/invitro_passage_structural_identifiability_v5.html"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def canonical(event):
    return event.replace("SUM159_NLS", "SUM-159_NLS")


def freeze(out):
    start = process_time()
    out.mkdir(parents=True, exist_ok=True)
    if (out/"freeze_resources.tsv").exists():
        raise ValueError("The fitted point is already frozen.")
    inputs = [PROFILES / name for name in ("final_baselines.tsv", "parameter_trajectories.tsv",
              "unrestricted_anchor_receipts.tsv", "profile_model.py", "runtime.py", "step03.dylib")]
    inputs += [FIT / name for name in ("selected_primary.tsv", "constraint_margins.tsv", "joint.dylib",
               "joint_kernel.cpp", "a7k_event_contract.tsv", "karyotype_distributions.tsv", "readout_predictions.tsv")]
    inputs += [ALL / name for name in ("baseline_vector.tsv", "domain.tsv")]
    inputs += [FROZEN_DESIGN / name for name in ("passage_graph.tsv", "event_mask.tsv", "saved_endpoint_signal_replay.tsv")]
    receipts = []
    for path in inputs:
        index = path.parent / "SHA256SUMS"
        hashes = {line.split(maxsplit=1)[1].lstrip("*"): line.split()[0] for line in index.read_text().splitlines()}
        actual = digest(path)
        assert actual == hashes[path.name], path
        receipts.append(dict(path=str(path.relative_to(ROOT)), sha256=actual, basis="Sealed source file"))
    receipts.append(dict(path=str(HTML.relative_to(ROOT)), sha256=digest(HTML), basis="Presentation at freeze"))
    write(out/"source_receipts.tsv", receipts)
    anchor, = [r for r in read(PROFILES/"final_baselines.tsv") if r["contract"] == "J"]
    assert anchor["baseline_status"] == "verified"
    rows = [r for r in read(PROFILES/"parameter_trajectories.tsv") if r["contract"] == "J"
            and r["point"] == "unrestricted_anchor" and r["attempt_id"] == anchor["attempt_id"]]
    point = {r["parameter"]:r for r in rows}
    assert len(point) == len(rows) == 21
    old, = read(FIT/"selected_primary.tsv")
    later, = read(ALL/"baseline_vector.tsv")
    domain = {r["parameter"]: r for r in read(ALL/"domain.tsv")}
    model = ProfileModel("J")
    conf = ProfileConfig(model)
    assert set(conf.names) == set(point) == set(domain)
    x = np.array([float(point[n]["normalized"]) for n in conf.names])
    p, _ = conf.unmap(x)
    active = []
    for i,n in enumerate(conf.names):
        r = point[n]
        v = float(r["natural"])
        assert isclose(p[conf.active[i]], v, abs_tol=1e-15, rel_tol=2e-14), n
        p[conf.active[i]] = v
        limit = domain[n]
        lo, hi = float(limit["natural_lower"]), float(limit["natural_upper"])
        assert lo <= v <= hi and 0 <= x[i] <= 1, n
        assert isclose(v, float(later[n]), rel_tol=1e-10, abs_tol=1e-14)
        face = "lower" if v==lo else "upper" if v==hi else "interior"
        direction = {"lower":"Increase or stay", "upper":"Decrease or stay", "interior":"Both signs locally"}[face]
        role = "Mechanism" if i < 14 else ("Signal observer" if n in ("eta_H_2N","eta_M_2N","class_log_odds","misseg_dead_share") else "Observation law")
        active.append(dict(parameter=n, role=role, native_slot=int(conf.active[i]), natural_source=r["natural"],
                           exact_binary_value=str(Fraction.from_float(v)), normalized=r["normalized"],
                           lower=lo, upper=hi, reported_face=r["bound"] or "interior",
                           exact_face=face, permitted_direction=direction, earlier_fit=old[n], later_reference=later[n]))
    p[16:18] = p[14:16]
    score = model.score(p)
    assert abs(score["loss"]-float(anchor["loss"])) < 1e-6
    write(out/"fitted_point.tsv", active)
    write(out/"boundary_directions.tsv", [r for r in active if r["reported_face"]!="interior"])
    write(out/"native_vector.tsv", [dict(slot=i,value=repr(float(v)),exact_binary_value=str(Fraction.from_float(float(v)))) for i,v in enumerate(p)])
    write(out/"fit_replay.tsv", [dict(attempt=anchor["attempt_id"], registered_J=anchor["loss"], replay_J=score["loss"],
          sigma_profiled=score["sigma"], earlier_attempt=old["attempt_id"], later_reference_J=later["loss"])])
    graph = read(FROZEN_DESIGN/"passage_graph.tsv")
    groups = defaultdict(list)
    for g in graph:
        if g["donor_passage"] and not g["historical_supply"].startswith("stock transfer"):
            groups[g["donor_passage"], round(20*float(g["donor_model_day"]))].append(g)
    order = {r["passage_id"]:i for i,r in enumerate(graph)}
    keys = sorted(groups, key=lambda k:(order[k[0]],k[1]))
    gaps,_=model.supply(p)
    oldgaps={int(r["constraint_index"]):float(r["value"]) for r in read(FIT/"constraint_margins.tsv") if r["case"]=="selected_primary"}
    assert len(keys)==len(gaps)==len(oldgaps)==110 and np.argmin(gaps)==49
    ledger=[]
    for i,((donor,tick),gap) in enumerate(zip(keys,gaps)):
        children=groups[donor,tick]
        assert gap>0 and abs(gap-oldgaps[i])<1e-10*max(1.,abs(gap))
        ledger.append(dict(index=i,donor=donor,tick=tick,day=tick/20,donor_record=children[0]["donor_record"],
            timestamp=children[0]["donor_timestamp"],children="|".join(c["passage_id"] for c in children),
            demand_cells=sum(int(float(c["fixed_seed_live_objects"])) for c in children),
            native_relative=gap,earlier_fit_relative=oldgaps[i]))
    write(out/"native_supply.tsv",ledger)
    def predictions(vector):
        channels,states=model.readouts(vector)
        signals={}
        eps=1e-12
        for slot in model.indices:
            event=canonical(model.slot_ids[slot])
            signals["C|"+event]=sum(channels[:,slot])*1e12
        for slot in model.death_rows:
            live,h,m=channels[:,slot]
            g=vector[55]+log((h+vector[57]*m+eps)/(live+(1-vector[57])*m+eps))
            signals["D|"+canonical(model.slot_ids[slot])]=eps+(1-2*eps)/(1+exp(-g))
        for event,state in zip(model.all_kary,states):
            probabilities=(1-eps)*state/state.sum()+eps/133
            assert abs(sum(probabilities)-1)<1e-12
            for n,v in enumerate(probabilities,22):
                signals[f"K|{canonical(event)}|{n}"]=v
        return signals
    signals=predictions(p)
    previous=p.copy()
    for i,n in enumerate(conf.names):
        previous[conf.active[i]]=float(old[n])
    previous[16:18]=previous[14:16]
    prior_replay=predictions(previous)
    saved={}
    for r in read(FROZEN_DESIGN/"saved_endpoint_signal_replay.tsv"):
        kind={"visible_count":"C","death_probability_zero_residual":"D"}[r["readout"]]
        saved[kind+"|"+r["event"]]=float(r["predicted"])
    for r in read(FIT/"karyotype_distributions.tsv"):
        if r["mode"]=="primary":
            saved[f"K|{canonical(r['event_id'])}|{r['autosomes']}"]=float(r["fitted_probability"])
    assert len(signals)==1218 and signals.keys()==saved.keys()
    evidence=[]
    for name,v in signals.items():
        reference=prior_replay[name]
        discrepancy=abs(reference/saved[name]-1) if name.startswith("C|") else abs(reference-saved[name])
        assert discrepancy<1e-10,(name,discrepancy)
        evidence.append(dict(output=name,native_diamond=v,earlier_archived=saved[name],
                             earlier_point_replay_error=discrepancy,between_fit_change=v-reference))
    write(out/"signal_replay.tsv",evidence)
    write(out/"freeze_resources.tsv",[dict(cpu_seconds=process_time()-start)])
    print("Frozen",anchor["attempt_id"],"with the original 21 coordinates.",ledger[49],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=output_directory)
    freeze(parser.parse_args().out)
