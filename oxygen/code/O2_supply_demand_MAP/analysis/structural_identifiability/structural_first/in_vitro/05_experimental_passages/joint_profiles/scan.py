"""Fit fixed-coordinate profiles on the frozen three-readout J contract."""

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
from math import inf, isclose
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import monotonic

ROOT = next(p for p in Path(__file__).resolve().parents if (p / ".git").exists())
ARCHIVES = ROOT / "oxygen/results/synthetic_identifiability"
STAGE3 = ARCHIVES / "20261122_experimental_joint_identifiability_step03_01"
REVIEW = ARCHIVES / "20261122_experimental_joint_identifiability_limited_review_01"
DESIGN = ARCHIVES / "20261122_in_vitro_joint_structural_steps01_03_01"

sys.dont_write_bytecode = True
for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[key] = "1"
sys.path[:0] = [str(REVIEW), str(STAGE3)]
from runtime import AGREEMENT, BASE_REPAIR, SCORE_TOLERANCES, STATIONARY, STEP2, cpu
from profile_model import GRIDS, SUSPENDED, ProfileConfig, ProfileModel, independently
from profile_run import MATERIAL_SEPARATION, apply_fixed, choose_start, known
from review_support import polished_solve
from solve_local import difference, exact_margins, normalized_constraint, stationarity
from tv_review import fingerprint
import numpy as np
import scipy

PILOT = ("p_misseg", "buffer_beta", "k_o_mis", "n_O")
ARCHIVED = "archived Step 3"
SCORE_KEYS = ("loss", "count_component", "death_component", "kary_component", "sigma")
STEPS = (1e-5, 1e-6, 1e-7)


def read(path):
    if not path.is_file():
        return []
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, rows, append=False):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists() and path.stat().st_size > 0
    pending = path if append else path.with_suffix(".tmp")
    with pending.open("a" if append else "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n")
        if not append or not exists:
            writer.writeheader()
        writer.writerows(rows)
        stream.flush()
    if not append:
        pending.replace(path)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


class Config(ProfileConfig):
    @property
    def targets(self):
        return self.names[:14]

    def natural(self, x, name):
        i = self.idx[name]
        return float(self.unmap(x)[0][self.active[i]])

    def target_coordinate(self, name, value):
        if name in self.transforms:
            return float(self.to_x(name, value))
        archived = self.base.original.base.base
        i = archived.names.index(name)
        mode, scale = archived.transform[i], archived.scale[i]
        encoded = np.log10(value) if mode == "log10" else (
            np.log1p(value/scale) if mode == "log1p" else value)
        return float(np.clip((encoded-archived.shift[i])/archived.span[i], 0., 1.))

    def constraints(self, name, value):
        if name is None:
            return {}
        assert name in self.targets
        fixed = super().fixed(name, value) if name in SUSPENDED else {
            self.idx[name]: self.target_coordinate(name, value)}
        x = apply_fixed(np.zeros(len(self.names)), fixed)
        assert isclose(self.natural(x, name), value, abs_tol=1e-13, rel_tol=1e-12), (name, value)
        if value == 0:
            assert self.natural(x, name) == 0
        return fixed

    def vector(self, row):
        return np.array([float(row["x_" + name]) for name in self.names])


def check_sources(out):
    inputs = []
    for packet in (STAGE3, STEP2, REVIEW, DESIGN):
        for line in (packet / "SHA256SUMS").read_text().splitlines():
            expected, name = line.split(maxsplit=1)
            path = packet / name
            actual = digest(path)
            assert actual == expected, path
            inputs.append(dict(source=str(path.relative_to(ROOT)), sha256=actual))
    listed = {row["source"] for row in inputs}
    for row in read(STAGE3 / "source_hashes.tsv"):
        label = row["path"]
        if " at " in label:
            repo, name = label.split(" at ", 1)
            actual = sha256(subprocess.check_output(["git", "-C", repo, "show", name])).hexdigest()
            key = name
        else:
            path = Path(label)
            actual = digest(path)
            key = str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)
        assert actual == row["actual"], key
        if key not in listed:
            inputs.append(dict(source=key, sha256=actual))
            listed.add(key)
    write(out / "source_receipt.tsv", inputs)
    write(out / "environment.tsv", [dict(item=k,value=v) for k,v in (
        ("git", subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()),
        ("platform", platform.platform()), ("python", platform.python_version()),
        ("numpy", np.__version__), ("scipy", scipy.__version__), ("numerical_threads_per_process", 1))])


def independent_score(model, p, saved=None):
    checked = independently(model, p)
    current = model.score(p)
    for key, tolerance in SCORE_TOLERANCES.items():
        assert abs(current[key] - checked[key]) <= tolerance, (key, current[key], checked[key])
    for key in ("count_component", "death_component", "kary_component", "sigma"):
        assert isclose(current[key], checked[key], rel_tol=1e-11, abs_tol=1e-7), key
    if saved is not None:
        for key in ("loss","count_component","death_component","kary_component"):
            assert isclose(checked[key], float(saved[key]), abs_tol=1e-6, rel_tol=0), (key, saved)
    return checked


def local_audit(model, conf, x, fixed, diag):
    p, jac = conf.unmap(x)
    raw = model.score(p, True)
    g = conf.gradient(raw["gradient"], jac)
    old = raw["loss"]
    margins = exact_margins(model, p)
    active = diag["active"]
    face = [i for i in active if x[i] <= 5e-7 or x[i] >= 1-5e-7]
    dirs = [(conf.names[i], np.eye(len(x))[i]*(1 if x[i] < .5 else -1)) for i in face]
    residual = np.zeros(len(x))
    residual[active] = -diag["residual"]
    dirs.append(("stationarity_residual", residual))
    interior = np.array([i for i in active if i not in face], dtype=int)
    if len(interior):
        donors = np.flatnonzero(margins <= 5e-7)
        _, matrix = normalized_constraint(model, conf, p, jac, interior, donors, [])
        projected = g[interior]
        if len(matrix):
            _, s, vt = np.linalg.svd(matrix, full_matrices=True)
            basis = vt[int(np.sum(s > 1e-9*max(1.,s[0]))):].T
            projected = basis @ (basis.T @ projected)
        tangent = np.zeros(len(x))
        tangent[interior] = -projected
        dirs.append(("active_tangent", tangent))
    ledger = []
    descent = False
    for name, d in dirs:
        norm = np.linalg.norm(d)
        if norm == 0.:
            continue
        drops = []
        for step in STEPS:
            y = apply_fixed(x + step*d/norm, fixed)
            q = conf.unmap(y)[0]
            margin = float(min(exact_margins(model, q)))
            if margin >= 0 and not np.array_equal(x,y):
                delta = model.score(q)["loss"]-old
                drops.append(delta)
                ledger.append(dict(direction=name,step=step,minimum_margin=margin,objective_change=delta))
        descent |= sum(v < -1e-6 for v in drops) >= 2
    return not descent, ledger


def check_zero(model, conf, x, target, fixed):
    carrier = SUSPENDED[target]
    assert conf.natural(x, target) == 0.0 and conf.idx[carrier] in fixed
    p = conf.unmap(x)[0]
    old = independent_score(model, p)
    margins = exact_margins(model, p)
    max_delta = 0.
    for v in (0., .5, 1.):
        y = x.copy()
        y[conf.idx[carrier]] = v
        q = conf.unmap(y)[0]
        assert conf.natural(y, target) == 0.0
        alt = independently(model, q)
        errors = [abs(old["loss"]-alt["loss"]), *(float(np.max(abs(old[k]-alt[k])))
                  for k in ("channels","fraction","probability")),
                  float(np.max(abs(margins-exact_margins(model,q))))]
        max_delta = max(max_delta, *errors)
    assert max_delta <= 1e-9, (target, max_delta)
    return max_delta


def zero_face_audit(model, conf, x, fixed, loss):
    """Probe a source mechanism releasing from zero at different carrier values."""
    unresolved = []
    for target, carrier in SUSPENDED.items():
        i, j = conf.idx[target], conf.idx[carrier]
        if i in fixed or x[i] != 0.:
            continue
        for value in (0., .5, 1.) if j not in fixed else (x[j],):
            for step in (1e-6,1e-5):
                y = x.copy()
                y[j], y[i] = value, step
                p = conf.unmap(y)[0]
                if min(exact_margins(model,p)) >= 0 and model.score(p)["loss"] < loss-1e-6:
                    unresolved.append(f"{target}:{value:g}")
                    break
    return " ".join(unresolved)


class Run:
    def __init__(self, out):
        self.out = out
        self.model = ProfileModel("J")
        self.conf = Config(self.model)
        self.history = read(STAGE3 / "profile_endpoints.tsv")
        self.archived_points = {r["attempt_id"]: r for r in self.history if r["contract"]=="J"}
        self.archived_starts = defaultdict(dict)
        for row in read(STAGE3 / "starts.tsv"):
            if row["attempt_id"] in self.archived_points:
                self.archived_starts[row["attempt_id"]][row["name"]] = float(row["normalized"])

    def start(self, row):
        return np.array([float(row["start_"+n]) for n in self.conf.names])

    def log(self, folder, row):
        write(folder/"attempts.tsv", [row], append=True)

    def accept(self, x, start, fixed, target, value, label, source, origin, metrics, saved=None):
        model, conf = self.model, self.conf
        p = conf.unmap(x)[0]
        measured = independent_score(model, p, saved)
        exact = float(min(exact_margins(model,p)))
        diag = stationarity(model, conf, x, fixed)
        points_match = target is None or isclose(conf.natural(x,target), value, rel_tol=1e-12, abs_tol=1e-13)
        assert points_match
        zero_error = check_zero(model,conf,x,target,fixed) if target in SUSPENDED and value==0. else ""
        audit, audits = local_audit(model,conf,x,fixed,diag) if exact>=0 and diag["norm"]<=STATIONARY else (False,[])
        unusual = zero_face_audit(model,conf,x,fixed,measured["loss"]) if audit else ""
        faces = " ".join(f"{conf.names[i]}:{'lower' if v <= 5e-7 else 'upper'}"
                         for i,v in enumerate(x) if i not in fixed and min(v,1.-v)<=5e-7)
        inactive = " ".join(carrier for t,carrier in SUSPENDED.items()
                            if conf.natural(x,t)==0.)
        qualified = exact>=0 and diag["norm"]<=STATIONARY and audit and not unusual
        row = dict(attempt=label, parameter=target or "", value="" if value is None else value,
                   source=source, origin=origin, outcome="qualified" if qualified else "provisional",
                   details="", minimum_margin=exact, stationarity=diag["norm"], local_check=audit,
                   zero_invariance_error=zero_error, zero_release_descent=unusual,
                   inactive_at_endpoint=inactive, active_box_faces=faces, active_constraints=" ".join(diag["labels"]),
                   native_support=" ".join(f"{m:x}" for m in fingerprint(model,p)),
                   **metrics, **{key:measured[key] for key in SCORE_KEYS},
                   **{"x_"+name:float(v) for name,v in zip(conf.names,x)},
                   **{name:float(p[i]) for name,i in zip(conf.names,conf.active)},
                   **{"start_"+name:float(v) for name,v in zip(conf.names,start)})
        for a in audits:
            a["attempt"] = label
        return row, audits

    def execute(self, folder, label, target, value, preferred, used, variant):
        state = {r["attempt"]:r for r in read(folder/"attempts.tsv")}
        if label in state:
            return state[label]
        conf, model = self.conf,self.model
        fixed = conf.constraints(target,value)
        seed = int.from_bytes(sha256(label.encode()).digest()[:4],"little")
        begin, wall = cpu(),monotonic()
        chosen = choose_start(model,conf,fixed,preferred,used,inf,seed)
        if chosen is None:
            row = self.failed(label,target,value,"start_search","No distinct feasible start found",cpu()-begin,monotonic()-wall)
            self.log(folder,row)
            return row
        source, start, _ = chosen
        try:
            y, fit, _, polish = polished_solve(model,conf,start,fixed,inf,variant)
            row, audits = self.accept(y,start,fixed,target,value,label,source,"fit",
                   dict(cpu_seconds=cpu()-begin, wall_seconds=monotonic()-wall, calls=fit.calls, iterations=fit.iterations,
                        solver_rounds=fit.rounds, solver_message=fit.stop, observation_polish_steps=polish))
            row.update(cpu_seconds=cpu()-begin,wall_seconds=monotonic()-wall)
            if audits:
                write(folder/"direction_checks.tsv",audits,append=True)
        except (ArithmeticError, np.linalg.LinAlgError) as error:
            row = self.failed(label,target,value,"fit",str(error),cpu()-begin,monotonic()-wall,source,start)
        self.log(folder,row)
        print(f"{label} {row['outcome']} {row['loss']} KKT={row['stationarity']} CPU={float(row['cpu_seconds']):.1f}s",flush=True)
        return row

    def failed(self,label,target,value,origin,message,seconds,wall,source="",start=None):
        return dict(attempt=label,parameter=target or "",value="" if value is None else value,source=source,
                    origin=origin,outcome="failed",details=message,minimum_margin="",stationarity="",local_check="",
                    zero_invariance_error="",zero_release_descent="",inactive_at_endpoint="",active_box_faces="",
                    active_constraints="",native_support="",cpu_seconds=seconds,wall_seconds=wall,calls="",
                    iterations="",solver_rounds="",solver_message="",observation_polish_steps="",
                    **{key:"" for key in SCORE_KEYS},**{"x_"+n:"" for n in self.conf.names},
                    **{n:"" for n in self.conf.names},
                    **{"start_"+n:"" if start is None else float(start[i]) for i,n in enumerate(self.conf.names)})

    def import_archive(self,folder,saved):
        target = saved["parameter"] or None
        value = float(saved["fixed_value"]) if target else None
        label = "archive_"+saved["attempt_id"]
        if any(r["attempt"]==label for r in read(folder/"attempts.tsv")):
            return
        conf = self.conf
        begin, wall = cpu(),monotonic()
        start = np.array([self.archived_starts[saved["attempt_id"]][n] for n in conf.names])
        row,audits = self.accept(conf.vector(saved),start,conf.constraints(target,value),target,value,label,
                               saved["source"],ARCHIVED,dict(cpu_seconds=0.,wall_seconds=0.,calls=0,iterations=0,
                                    solver_rounds=0,solver_message="Independent replay",observation_polish_steps=0),
                               saved=saved)
        row.update(cpu_seconds=cpu()-begin,wall_seconds=monotonic()-wall)
        if audits:
            write(folder/"direction_checks.tsv",audits,append=True)
        self.log(folder,row)

    def verdict(self, rows, minimum=2):
        fits = [r for r in rows if r["loss"]!="" and float(r["minimum_margin"])>=0]
        floor = min(fits,key=lambda r:float(r["loss"])) if fits else None
        qualified = [r for r in sorted(fits,key=lambda r:float(r["loss"])) if
                     r["outcome"]=="qualified" and float(r["loss"])<=float(floor["loss"])+AGREEMENT]
        separate = []
        keys = ["start_"+n for n in self.conf.names]
        for r in qualified:
            x = [float(r[k]) for k in keys]
            if all(difference(x,[float(prior[k]) for k in keys])>MATERIAL_SEPARATION for prior in separate):
                separate.append(r)
        status = "verified" if len(separate)>=minimum else ("provisional" if floor else "feasibility_unresolved")
        return status,floor,separate

    def baseline(self):
        folder = self.out/"baseline"
        selected, = [r for r in read(STAGE3/"final_baselines.tsv") if r["contract"]=="J"]
        old = self.archived_points[selected["attempt_id"]]
        assert selected["baseline_status"]=="verified"
        assert len(self.model.y)==225 and len(self.model.dead)==195 and self.model.kary.shape==(6,133)
        assert self.model.kary.sum()==104
        assert len(exact_margins(self.model,self.conf.unmap(self.conf.vector(old))[0]))==110
        assert all(r["included"]=="False" for r in read(STAGE3/"admission.tsv"))
        for row in self.archived_points.values():
            if not row["parameter"] and row["qualified"]=="True" and float(row["loss"])<=float(old["loss"])+AGREEMENT:
                self.import_archive(folder,row)
        self.exact_zero_preflight(self.conf.vector(old))
        genuine = [(k,self.archived_start_vector(k)) for k in self.archived_starts
                   if k.startswith("J_base_") and k in self.archived_points]
        bank = genuine + known(self.conf)
        used = []
        for i in range(3):
            row = self.execute(folder,f"baseline_fresh_{i+1}",None,None,bank,used,i)
            if row["start_"+self.conf.names[0]]!="":
                used.append(self.start(row))
        self.finish_baseline()
        self.freeze_domain()

    def finish_baseline(self):
        folder = self.out/"baseline"
        rows = read(folder/"attempts.tsv")
        status,floor,agree = self.verdict(rows)
        assert status=="verified", (status, floor, agree)
        assert any(r["origin"]=="fit" for r in agree), "No fresh baseline fit in reproduced region"
        hit = agree[0]["active_box_faces"].split()
        write(self.out/"baseline.tsv",[dict(status=status,agreeing_starts=len(agree),
                    minimum_feasible=floor["loss"],selected_attempt=agree[0]["attempt"],
                    selected_loss=agree[0]["loss"],mechanistic_box_faces=" ".join(k for k in hit if k.split(":")[0] in self.conf.targets))])
        write(self.out/"baseline_vector.tsv",[agree[0]])
        print(f"Baseline {status}: {agree[0]['loss']}, {len(agree)} distinct qualified starts.",flush=True)

    def archived_start_vector(self,key):
        return np.array([self.archived_starts[key][n] for n in self.conf.names])

    def exact_zero_preflight(self, anchor):
        data = []
        for target in SUSPENDED:
            fixed = self.conf.constraints(target,0.)
            trial = apply_fixed(anchor, fixed)
            p = self.conf.unmap(trial)[0]
            max_delta = check_zero(self.model,self.conf,trial,target,fixed)
            data.append(dict(parameter=target,exact_native_zero=self.conf.natural(trial,target)==0.0,
                             inactive_carrier=SUSPENDED[target],maximum_difference=max_delta,passed=True))
        write(self.out/"exact_zero_preflight.tsv",data)

    def freeze_domain(self):
        data=[]
        for i,n in enumerate(self.conf.names):
            role = "mechanism" if n in self.conf.targets else "refitted_observation"
            lo,hi=(self.conf.natural(np.full(len(self.conf.names),x),n) for x in (0.,1.))
            data.append(dict(parameter=n,role=role,natural_lower=lo,natural_upper=hi,
                             exact_zero_carrier=SUSPENDED.get(n,"")))
        write(self.out/"domain.tsv",data)

    def baseline_x(self):
        row, = read(self.out/"baseline_vector.tsv")
        return self.conf.vector(row)

    def grid(self,name,pilot=False):
        c = self.conf
        anchor = self.baseline_x()
        x = anchor[c.idx[name]]
        def natural(q):
            a=anchor.copy()
            a[c.idx[name]]=q
            return c.natural(a,name)
        if pilot:
            values = {"p_misseg": (0.,.005), "buffer_beta": (0.,.25),
                      "k_o_mis": (natural(0.),natural(.045)),
                      "n_O": (natural(x),natural(max(.0,x-.07)))}[name]
            return sorted(set(values))
        coordinates = list(np.linspace(0.,1.,10))
        extra = list(GRIDS.get(name,()))
        if not any(abs(natural(x)-v)<=5e-12*max(1.,abs(v)) for v in extra):
            extra.append(natural(x))
        for v in extra:
            q = c.target_coordinate(name,v)
            coordinates = [z for z in coordinates if abs(z-q)>1e-9]
        return sorted(set(extra+[natural(q) for q in coordinates]))

    def series(self,name,pilot,refine=False):
        assert name in self.conf.targets
        folder = self.out/"profiles"/name
        grid = self.grid(name,pilot)
        if refine:
            planned = read(folder/"refinement_plan.tsv")
            if planned:
                grid = [float(r["value"]) for r in planned]
            else:
                existing = read(folder/"grid.tsv")
                assert len(existing)>=10
                z = self.baseline_x()
                center = z[self.conf.idx[name]]
                grid = [float(r["value"]) for r in existing+read(folder/"pilot_grid.tsv")]
                for sign in (-1,1):
                    options=[float(r["normalized"]) for r in existing
                             if (float(r["normalized"])-center)*sign>1e-8]
                    if options:
                        adjacent=min(options,key=lambda value:abs(value-center))
                        if abs(adjacent-center)>.055:
                            z[self.conf.idx[name]]=(center+adjacent)/2.
                            grid.append(self.conf.natural(z,name))
                grid=sorted(set(grid))
                write(folder/"refinement_plan.tsv",
                      [dict(value=v,baseline_coordinate=center) for v in grid])
        for saved in self.archived_points.values():
            if saved["parameter"]==name and float(saved["fixed_value"]) in grid:
                self.import_archive(folder,saved)
        base = self.baseline_x()
        bank = [("accepted J baseline",base)]+known(self.conf)
        current = read(folder/"attempts.tsv")
        def states(value):
            return [r for r in current if float(r["value"])==value]
        def best(value):
            return self.verdict(states(value))[1]
        def neighbors(value, side):
            other = sorted((v for v in grid if (v-value)*side>0),key=lambda v:abs(v-value))
            return [(f"neighbor {name} {v:.16g}",self.conf.vector(row))
                    for v in other if (row:=best(v)) is not None][:2]
        focus = self.conf.natural(base,name)
        ascending = sorted(grid)
        center_out = sorted(grid,key=lambda v:abs(self.conf.target_coordinate(name,v)-self.conf.target_coordinate(name,focus)))
        schedule = (("outward",center_out,0),("upward",ascending,-1),("downward",ascending[::-1],1))
        for label,order,side in schedule:
            for v in order:
                found = states(v)
                status, floor, agree = self.verdict(found)
                if status=="verified":
                    continue
                used = [self.start(r) for r in found if r["start_"+self.conf.names[0]]!=""]
                both = neighbors(v,-1)+neighbors(v,1)
                chosen = neighbors(v,side) if side else both
                preferred = (chosen+bank+both) if label!="downward" else (chosen+bank[1:]+bank[:1]+both)
                point = f"{name}_{format(v,'.17g')}_{label}"
                self.execute(folder,point,name,v,preferred,used,len(used))
                current = read(folder/"attempts.tsv")
        for v in grid:
            found = states(v)
            status,floor,agree = self.verdict(found)
            if status=="verified" or floor is None:
                continue
            used = [self.start(r) for r in found if r["start_"+self.conf.names[0]]!=""]
            preferred = [(f"feasible candidate {floor['attempt']}",self.conf.vector(floor))]+neighbors(v,-1)+neighbors(v,1)+bank
            label = f"{name}_{format(v,'.17g')}_reserve"
            self.execute(folder,label,name,v,preferred,used,len(used))
            current = read(folder/"attempts.tsv")
        self.collect_parameter(name,grid)
        if pilot:
            self.collect_parameter(name,grid,"pilot_grid.tsv")

    def collect_parameter(self,name,grid,filename="grid.tsv"):
        folder = self.out/"profiles"/name
        current = read(folder/"attempts.tsv")
        baseline, = read(self.out/"baseline.tsv")
        loss = float(baseline["selected_loss"])
        data=[]
        for v in grid:
            subset = [r for r in current if float(r["value"])==v]
            status,floor,agree = self.verdict(subset)
            chosen = agree[0] if status=="verified" else floor
            data.append(dict(parameter=name,value=v,normalized=self.conf.target_coordinate(name,v),status=status,
                exact_zero=v==0.,baseline_loss=loss,loss=chosen["loss"] if chosen else "",
                rise=float(chosen["loss"])-loss if chosen else "",
                feasible_floor=floor["loss"] if floor else "", agreeing_starts=len(agree),
                attempt=chosen["attempt"] if chosen else "",stationarity=chosen["stationarity"] if chosen else "",
                active_box_faces=chosen["active_box_faces"] if chosen else "",
                minimum_margin=chosen["minimum_margin"] if chosen else "",
                count_component=chosen["count_component"] if chosen else "",
                death_component=chosen["death_component"] if chosen else "",
                kary_component=chosen["kary_component"] if chosen else "",
                inactive_at_endpoint=" ".join(carrier for t,carrier in SUSPENDED.items()
                                               if chosen and float(chosen[t])==0.),
                fits=sum(r["origin"]=="fit" for r in subset),
                provenance="Replayed and newly fit" if any(r["origin"]==ARCHIVED for r in subset) else "New fit"))
        write(folder/filename,data)
        counts={s:sum(r["status"]==s for r in data) for s in ("verified","provisional","feasibility_unresolved")}
        print(f"{name}: {len(data)} grid points, {counts}.",flush=True)

    def reconcile(self):
        rows = [r for n in self.conf.targets for r in read(self.out/"profiles"/n/"grid.tsv")]
        current = [r for n in self.conf.targets for r in read(self.out/"profiles"/n/"attempts.tsv")]
        baseline, = read(self.out/"baseline.tsv")
        base = float(baseline["selected_loss"])
        lower = sorted((r for r in current if r["loss"] and
                       float(r["loss"]) < base-BASE_REPAIR), key=lambda r:float(r["loss"]))
        if lower:
            bank = [(r["attempt"],self.conf.vector(r)) for r in lower[:4]]
            bank += known(self.conf)
            f = self.out/"baseline"
            count = sum(r["attempt"].startswith("baseline_repair_") for r in read(f/"attempts.tsv"))
            used=[]
            for i in range(3):
                r = self.execute(f,f"baseline_repair_{count+i+1}",None,None,bank,used,i)
                if r["start_"+self.conf.names[0]]!="":
                    used.append(self.start(r))
            self.finish_baseline()
            baseline, = read(self.out/"baseline.tsv")
            assert float(baseline["selected_loss"]) <= float(lower[0]["loss"])+BASE_REPAIR, "Lower feasible score remains unresolved"
            for n in self.conf.targets:
                prior = [float(r["value"]) for r in rows if r["parameter"]==n]
                if prior:
                    self.collect_parameter(n,prior)
        rows = [r for n in self.conf.targets for r in read(self.out/"profiles"/n/"grid.tsv")]
        write(self.out/"profile_grid.tsv",rows)


def main():
    command = argparse.ArgumentParser()
    command.add_argument("stage", choices=("baseline","pilot","profile","refine","reconcile"))
    command.add_argument("out",type=Path)
    command.add_argument("parameter",nargs="?")
    args = command.parse_args()
    if args.stage=="baseline":
        check_sources(args.out)
    run = Run(args.out)
    if args.stage=="baseline":
        run.baseline()
    elif args.stage=="reconcile":
        run.reconcile()
    else:
        assert args.stage!="pilot" or args.parameter in PILOT
        run.series(args.parameter,args.stage=="pilot",args.stage=="refine")


if __name__=="__main__":
    main()
