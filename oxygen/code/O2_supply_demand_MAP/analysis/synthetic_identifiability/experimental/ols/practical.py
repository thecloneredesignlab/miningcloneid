"""Conditional experimental OLS likelihood with a source-checked native adjoint."""
import argparse
from collections import Counter, defaultdict
import csv
import ctypes as C
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import nnls
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
STRUCT = RESULTS / "structural_identifiability"
DESIGN = STRUCT / "20261118_culture_experimental_inheritance_01"
STEP2 = STRUCT / "20261118_culture_experimental_OLS_step02"
CONTRACT = HERE.parent / "20261113_culture_contract_steps01_03"
OLD = STRUCT / "20261117_culture_actual_schedule_step02"
MECHANISMS = ("lam_max", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
              "buffer_smax", "buffer_beta", "buffer_n_exp", "alpha_o2",
              "gamma_growth", "mu_hp", "gamma_mu", "O2_crit", "n_O")
INITIALS = ("init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N")
ETA = ("eta_H_2N", "eta_M_2N", "eta_H_4N", "eta_M_4N")
FD = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
FI = np.ctypeslib.ndpointer(dtype=np.int32, flags="C_CONTIGUOUS")
FU = np.ctypeslib.ndpointer(dtype=np.uint64, flags="C_CONTIGUOUS")


def read(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows):
    rows = list(rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    temp.replace(path)


def base_table():
    ordered = read(DESIGN / "hist/passage_design.tsv")
    passages = {p["passage_id"]: i for i, p in enumerate(ordered)}
    obs = read(STEP2 / "experimental_ols_reproduction.tsv")
    ids = [r["passage_id"] for r in obs]
    slope_index = {p: i for i, p in enumerate(ids)}
    target = np.array([float(x["recorded_ols_per_day"]) for x in obs])
    slots = sorted(read(STEP2 / "ols_calendar_map.tsv"),
                   key=lambda p: (passages[p["passage_id"]], int(p["euler_step"])))
    outcomes = []
    weights = []
    groups = []
    last = defaultdict(int)
    variances = np.zeros(len(obs))
    for slot in slots:
        p = slot["passage_id"]
        t = int(slot["euler_step"])
        k = slope_index[p]
        w = int(slot["ols_weight_numerator"]) / int(slot["ols_weight_denominator"])
        last[p] = max(last[p], t)
        outcomes.append((passages[p], t, k))
        weights.append(w)
        groups.append(k)
        variances[k] += w*w if t else 0.
    for p in ordered:
        if p["parent_id"]:
            parent = p["parent_id"]
            last[parent] = max(last[parent], int(p["donor_step"]))
            assert passages[parent] < passages[p["passage_id"]]
    roots = {x["child_passage_id"]: x for x in read(DESIGN / "hist/restart_policy.tsv")
             if x["new_policy"].startswith("inherit")}
    tilts = tuple(roots)
    assert len(tilts) == 12
    oxygens = sorted({float(x["oxygen_pct"]) for x in ordered})
    meta = []
    counts = []
    for p in ordered:
        key = p["passage_id"]
        i = passages.get(p["parent_id"], -1)
        donor_time = int(p["donor_step"]) if i >= 0 else -1
        meta.append((oxygens.index(float(p["oxygen_pct"])), p["cohort"] == "4N",
                     i, donor_time, last[key], tilts.index(key) if key in tilts else -1))
        counts.append(float(p["fixed_live_cells"]))
    hist = np.zeros((2, 133))
    for p in read(DESIGN / "hist/dp_proxy_histograms.tsv"):
        hist[int(p["cohort"] == "4N"), int(p["chromosome_count"]) - 22] = int(p["cells"]) / int(p["total_cells"])
    assert len(ordered) == 116 and len(obs) == 114 and len(slots) == 335 and np.all(variances > 0)
    assert np.allclose(hist.sum(axis=1), 1.)
    return (ordered, ids, target, slots, np.array(outcomes, dtype=np.int32), np.array(weights),
            np.array(groups), variances, np.array(oxygens), np.array(meta, dtype=np.int32),
            np.array(counts), hist, tilts)


def source_param(tilts):
    # Keep the bounds of the OLS rank study. Use the existing fitting transforms.
    roles = {r["parameter"]: r for r in read(CONTRACT / "parameter_roles.tsv")}
    raw = {r["parameter"]: r for r in read(OLD / "interior_384/exact_parameters.tsv")}
    names = (*MECHANISMS, *ETA, *INITIALS,
             *(f"tilt_{t}:{p}" for p in tilts for t in ("mean", "dispersion")))
    from fractions import Fraction
    rows = []
    for key in names:
        if key in raw:
            r = raw[key]
            q = roles[key]
            bounds = tuple(float(Fraction(r[b])) for b in ("lower", "upper"))
            v = float(Fraction(r["exact_value"]))
            method, scale = q["transform"], float(q["scale"])
        elif key in ETA:
            bounds = (0., 1.)
            v = 1/3 if key.startswith("eta_H") else 2/3
            method, scale = "identity", 1.
        else:
            bounds, v = (-0.01, 0.01), 0.
            method, scale = "identity", 1.
        rows.append(dict(parameter=key, lower=bounds[0], upper=bounds[1],
                         transform=method, scale=scale, witness=v))
    return rows


class Kernel:
    def __init__(self):
        (self.graph, self.ids, self.y, self.slots, outcomes, self.weights, self.group, self.var,
         self.oxygens, meta, counts, hist, self.tilts) = base_table()
        self.rows = source_param(self.tilts)
        self.anchor = np.array([r["witness"] for r in self.rows])
        self.n = len(self.ids)
        self.scale = np.sqrt(self.var)
        self.capture = np.array([int(r["euler_step"]) > 0 for r in self.slots])
        self.lib = C.CDLL(os.environ.get("CULTURE_OLS_LIB", "/private/tmp/culture_ols_fast_constraints.dylib"))
        lib = self.lib
        lib.ols_design.argtypes = [C.c_int, C.c_int, C.c_int, FD, FI, FD, C.c_int, FI, FD, FD]
        lib.ols_design.restype = C.c_void_p
        lib.ols_eval.argtypes = [C.c_void_p, FD, C.c_int, C.c_void_p, FD, FD, C.c_void_p, FD, FU]
        lib.ols_eval.restype = C.c_int
        lib.ols_loss.argtypes = [C.c_void_p, FD, C.c_int, FD, FD, C.c_int, FD, FD]
        lib.ols_loss.restype = C.c_int
        lib.ols_destroy.argtypes = [C.c_void_p]
        self.handle = lib.ols_design(len(meta), len(self.oxygens), self.n, self.oxygens, meta, counts,
                                     len(outcomes), outcomes, self.weights, np.ascontiguousarray(hist))
        if not self.handle:
            raise RuntimeError("Could not initialize experimental OLS design")
        lib.ols_constrained_loss.argtypes = [C.c_void_p, FD, C.c_int, FD, FD, C.c_int, FD, C.c_void_p]
        lib.ols_constrained_loss.restype = C.c_int
        lib.ols_supply.argtypes = [C.c_void_p, FD, C.c_int, FD, C.c_int, C.c_void_p, C.c_void_p]
        lib.ols_supply.restype = C.c_int
        lib.ols_nconstraints.argtypes = [C.c_void_p]
        self.n_constraints = lib.ols_nconstraints(self.handle)
        assert self.n_constraints == 110

    def predict(self, p, gaussian=False, combination=None):
        logc = np.empty(len(self.slots))
        raw = np.empty(len(self.slots))
        diag = np.empty(2)
        mask = np.empty(len(self.oxygens), dtype=np.uint64)
        if combination is None:
            u, g = None, None
        else:
            u = np.ascontiguousarray(self.weights * combination[self.group])
            g = np.empty(46)
        ok = self.lib.ols_eval(self.handle, np.ascontiguousarray(p), gaussian,
                               None if u is None else u.ctypes.data, logc, raw,
                               None if g is None else g.ctypes.data, diag, mask)
        if not ok:
            raise ArithmeticError("Non-finite trajectory or zero donor")
        return np.bincount(self.group, self.weights*logc, minlength=self.n), raw, diag, mask, g

    def evaluate(self, p, target, conf):
        diag = np.empty(3)
        d = np.empty(46)
        ok = self.lib.ols_loss(self.handle, np.ascontiguousarray(p), conf.gaussian, target,
                                conf.precision, conf.precision.ndim == 2, diag, d)
        if not ok:
            return 1e6, np.zeros(46), np.nan
        return diag[0], d, diag[2]

    def constrained(self, p, target, conf, gradient=True):
        f = np.empty(2)
        g = np.empty(46) if gradient else None
        okay = self.lib.ols_constrained_loss(
            self.handle, np.ascontiguousarray(p), conf.gaussian, target,
            conf.precision, conf.precision.ndim == 2, f, None if g is None else g.ctypes.data
        )
        if not okay:
            return 1e6, g if g is None else np.zeros(46), -10.
        return f[0], g, f[1]

    def supply(self, p, gaussian, pick=None):
        margin = np.empty(self.n_constraints)
        pick = np.ascontiguousarray([] if pick is None else pick,dtype=np.int32)
        grad = np.empty((len(pick),46))
        okay = self.lib.ols_supply(self.handle,np.ascontiguousarray(p),gaussian,margin,len(pick),
                                    pick.ctypes.data,grad.ctypes.data if len(pick) else None)
        if not okay:
            return np.full(self.n_constraints,-10.),np.zeros_like(grad)
        return margin,grad


class Configuration:
    def __init__(self, kernel, gaussian=False, split=False, tilt=False, precision=None):
        self.gaussian, self.split, self.tilt = gaussian, split, tilt
        self.active = np.array(list(range(14)) + (list(range(14,18)) if split else [14,15])
                               + (list(range(18,22)) if gaussian else [])
                               + (list(range(22,46)) if tilt else []))
        self.root = kernel.anchor.copy()
        self.precision = np.ascontiguousarray(1/kernel.var if precision is None else precision)
        self.parameters = [kernel.rows[k] for k in self.active]
        self.names = [r["parameter"] for r in self.parameters]
        self.transform = [r["transform"] for r in self.parameters]
        self.scale = np.array([r["scale"] for r in self.parameters])
        lo, hi = (self.enc(np.array([r[b] for r in self.parameters])) for b in ("lower", "upper"))
        self.shift, self.span = lo, hi-lo
        self.log = np.array([t == "log10" for t in self.transform])
        self.log1p = np.array([t == "log1p" for t in self.transform])
        self.x_anchor = self.normal(self.root[self.active])

    def enc(self, natural):
        y = natural.copy()
        for j, t in enumerate(self.transform):
            if t == "log10":
                y[j] = np.log10(y[j])
            elif t == "log1p":
                y[j] = np.log1p(y[j]/self.scale[j])
        return y

    def normal(self, natural):
        return (self.enc(natural) - self.shift) / self.span

    def unmap(self, x):
        z = self.shift + self.span*x
        v = z.copy()
        v[self.log] = 10.**z[self.log]
        v[self.log1p] = self.scale[self.log1p]*np.expm1(z[self.log1p])
        deriv = self.span.copy()
        deriv[self.log] *= np.log(10.)*v[self.log]
        deriv[self.log1p] *= v[self.log1p] + self.scale[self.log1p]
        p = self.root.copy()
        p[self.active] = v
        if not self.split:
            p[16:18] = p[14:16]
        return p, deriv

    def gradient(self, d, transform_jac):
        ar = d[self.active]
        if not self.split:
            ar[14:16] += d[16:18]
        return ar*transform_jac


def projected(g, x):
    return x-np.clip(x-g, 0, 1)


def optimize(kernel, conf, target, x0, limit=400, fixed=None, yield_score=False):
    target = np.ascontiguousarray(target)
    start = np.clip(np.asarray(x0), 0, 1).copy()
    active = np.arange(len(start))
    if fixed is not None:
        start[fixed[0]] = fixed[1]
        active = active[active != fixed[0]]
    calls = 0
    best = [np.inf, None, None, np.nan]

    def evaluate(x):
        nonlocal calls
        calls += 1
        a = start.copy()
        a[active] = x
        p, j = conf.unmap(a)
        f, d, supply = kernel.evaluate(p, target, conf)
        g = conf.gradient(d, j)
        if f < best[0]:
            best[:] = [f, a.copy(), g.copy(), supply]
        return f, g[active]

    started = time.monotonic()
    res = minimize(evaluate, start[active], method="L-BFGS-B", jac=True,
                   bounds=[(0.,1.)] * len(active),
                   options=dict(maxiter=limit, maxfun=limit*3, gtol=1e-7, ftol=1e-11, maxls=24, maxcor=20))
    f, a, grad, supply = best
    if a is None:
        raise RuntimeError(str(res.message))
    return dict(x=a, loss=f, supply=supply, gradient=float(np.max(np.abs(projected(grad[active],a[active])))),
                calls=calls, iteration=res.nit, scipy_success=bool(res.success), termination=str(res.message),
                seconds=time.monotonic()-started)


def starting_points(kernel, config, seed, n=10, center=None):
    rng = np.random.default_rng(seed)
    center = config.x_anchor.copy() if center is None else np.asarray(center).copy()
    design = qmc.LatinHypercube(len(center), seed=rng).random(n)
    starts = [center, np.clip(center + rng.normal(0,.07,len(center)), 0,1)]
    for i, offset in enumerate(design):
        # Global low-discrepancy coverage plus the anchored basin.
        a = offset if i % 3 == 2 else 0.30*offset + 0.70*center
        starts.append(a)
    return starts


def solve_constrained(kernel, conf, target, x0, limit=100, fixed=None,
                      best_iterate_tolerance=1e-8):
    """Generate violated supply constraints, then validate all 110 margins."""
    start = np.clip(np.asarray(x0),0,1).copy()
    active = np.arange(len(start))
    if fixed is not None:
        start[fixed[0]] = fixed[1]
        active = active[active != fixed[0]]
    target = np.ascontiguousarray(target)
    t = time.monotonic()
    calls = 0
    stages = 0
    iterations = 0
    best = None
    slack = 1e-8

    def coords(x):
        a = start.copy()
        a[active] = x
        p,j = conf.unmap(a)
        return p,j,a

    def objective(x):
        nonlocal calls,best
        calls += 1
        p,j,a = coords(x)
        f,g,supply = kernel.constrained(p,target,conf)
        if supply >= -slack and (best is None or f<best[0]):
            best = (f,a.copy(),supply)
        return f,conf.gradient(g,j)[active]

    margin,_ = kernel.supply(conf.unmap(start)[0],conf.gaussian)
    smallest = np.argsort(margin)[:4]
    chosen = set(smallest[margin[smallest] < .08])
    search = start[active].copy()
    for _ in range(6):
        pick = np.array(sorted(chosen),dtype=np.int32)

        def gfun(x):
            p,_,_ = coords(x)
            g,_ = kernel.supply(p,conf.gaussian)
            return g[pick]-slack

        def gjac(x):
            p,j,_ = coords(x)
            _,g = kernel.supply(p,conf.gaussian,pick)
            return np.asarray([conf.gradient(z,j)[active] for z in g])

        res = minimize(objective,search,jac=True,bounds=[(0,1)]*len(active),
                       constraints=[dict(type="ineq",fun=gfun,jac=gjac)] if len(pick) else (),
                       method="SLSQP",options=dict(maxiter=limit,ftol=2e-9))
        stages+=1
        iterations+=res.nit
        search = np.clip(res.x,0,1)
        p,j,a = coords(search)
        margin,_=kernel.supply(p,conf.gaussian)
        new=set(np.flatnonzero(margin < -slack))-chosen
        if not new:
            break
        chosen.update(sorted(new,key=lambda i:margin[i])[:4])

    f,_,sup=kernel.constrained(p,target,conf,gradient=False)
    if best is not None and (sup < -slack or best[0] < f-best_iterate_tolerance):
        f,a,sup=best
        p,j=conf.unmap(a)
    f,g,sup=kernel.constrained(p,target,conf)
    margin,_ = kernel.supply(p,conf.gaussian)
    near=np.flatnonzero(margin<=5e-7)
    grad=conf.gradient(g,j)[active]
    _,grad_g=kernel.supply(p,conf.gaussian,near)
    directions=[conf.gradient(z,j)[active] for z in grad_g]
    for i,index in enumerate(active):
        if a[index] < 5e-7 or a[index] > 1-5e-7:
            v=np.zeros(len(active))
            v[i]=1 if a[index]<.5 else -1
            directions.append(v)
    if directions:
        _,stationarity=nnls(np.array(directions).T,grad,maxiter=1000)
    else:
        stationarity=np.linalg.norm(grad)
    return dict(x=a,loss=f,supply=sup,gradient=float(stationarity),calls=calls,
                iteration=iterations,scipy_success=bool(res.success),
                termination=f"{res.message}; {stages} constraint rounds; {len(near)} active supply",
                seconds=time.monotonic()-t)


def result_row(conf, fit, **labels):
    p, _ = conf.unmap(fit["x"])
    return dict(labels, loss=fit["loss"], min_aggregate_donor_margin_fraction=fit["supply"],
                max_projected_gradient=fit["gradient"], calls=fit["calls"],
                iterations=fit["iteration"], scipy_success=fit["scipy_success"],
                termination=fit["termination"], seconds=fit["seconds"],
                **{f"x_{name}": float(x) for name,x in zip(conf.names,fit["x"])},
                **{name: p[j] for name,j in zip(conf.names,conf.active)})


def row_x(conf, row):
    return np.array([float(row[f"x_{k}"]) for k in conf.names])


def validate(kernel):
    import verification
    return verification.check(kernel)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("validate","bench"))
    args = parser.parse_args()
    k = Kernel()
    if args.command == "validate":
        validate(k)
    else:
        cfg = Configuration(k)
        t = time.monotonic()
        p,_ = cfg.unmap(cfg.x_anchor)
        for _ in range(30):
            k.evaluate(p, k.y, cfg)
        print("seconds per complete objective with gradient", (time.monotonic()-t)/30, flush=True)
