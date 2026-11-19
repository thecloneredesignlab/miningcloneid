"""Direct likelihood for the audited post-inoculation counts."""
import ctypes as C
from datetime import date
from hashlib import sha256
import os
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
OLS = HERE.with_name("20261118_experimental_OLS_practical_01")
sys.path.insert(0, str(OLS))
os.environ["CULTURE_OLS_LIB"] = os.environ.get(
    "CULTURE_COUNT_LIB", "/private/tmp/culture_count_practical_01.dylib"
)
from practical import FD, FI, Kernel, STEP2, Configuration, read, solve_constrained, write

MAPPING = STEP2.with_name("20261118_culture_experimental_OLS_step01") / "count_event_mapping.tsv"
REPO = HERE.parents[3]
SNAPSHOT = REPO.parent / "soft_couping_org/data/InVitroData/cloneid_passaging_sum159_snapshot_20260731.tsv"


class CountKernel(Kernel):
    def __init__(self):
        super().__init__()
        self.ols_n, self.ols_y = self.n, self.y.copy()
        self.indices = np.flatnonzero(self.capture).astype(np.int32)
        self.n = len(self.indices)
        self.var = np.ones(self.n)
        self.scale = self.var.copy()
        self.lib.count_loss.argtypes = [C.c_void_p, FD, C.c_int, C.c_int, FI, FD, FD,
                                      C.c_int, FD, C.c_void_p]
        self.lib.count_loss.restype = C.c_int
        self.observations = read_observations(self)
        self.counts = np.array([float(r["corrected_count"]) for r in self.observations])
        self.observed_logs = np.log(self.counts)
        self.y = np.ascontiguousarray(self.observed_logs[self.capture])

    def slopes(self, log_counts):
        return np.bincount(self.group, self.weights * log_counts, minlength=self.ols_n)

    def log_predict(self, p, gaussian=False):
        _, counts, diagnostics, _, _ = self.predict(p, gaussian)
        return np.log(np.maximum(counts, 1.)), diagnostics

    def constrained(self, p, target, conf, gradient=True):
        f = np.empty(2)
        g = np.empty(46) if gradient else None
        okay = self.lib.count_loss(self.handle, np.ascontiguousarray(p), conf.gaussian, self.n,
                                  self.indices, target, conf.precision, conf.precision.ndim == 2,
                                  f, None if g is None else g.ctypes.data)
        if not okay:
            return 1e12, np.zeros(46) if gradient else None, -10.
        return f[0], g, f[1]


def read_observations(k):
    mapped = {r["record_id"]:r for r in read(MAPPING)}
    raw = read(SNAPSHOT)
    raw_by_id = {r["id"]:r for r in raw}
    assert len(raw_by_id) == len(raw) and len(mapped) == len(k.slots) == 335
    graph = {r["passage_id"]:r for r in k.graph}
    observations = []
    for slot in k.slots:
        record, passage = slot["record_id"], slot["passage_id"]
        row, source = mapped[record], raw_by_id[record]
        count = float(source["correctedCount"])
        assert np.isfinite(count) and count > 1
        assert float(row["db_correctedCount"]) == count
        assert row["count_timestamp_as_recorded"] == source["date"]
        age = (date.fromisoformat(source["date"][:10]) -
               date.fromisoformat(raw_by_id[passage]["date"][:10])).days
        assert age == int(slot["calendar_day"]) == int(row["ols_calendar_day"])
        assert 20 * age == int(slot["euler_step"]) and row["passage_id"] == passage
        assert slot["model_segment_id"] == row["model_segment_id"]
        if age == 0:
            assert record == passage and count == float(graph[passage]["fixed_live_cells"])
            assert source["event"] == "seeding"
        else:
            assert source["event"] == "harvest" and source["passaged_from_id1"] == passage
        observations.append(dict(record_id=record, passage_id=passage, cohort=slot["cohort"],
                                 oxygen_pct=slot["oxygen_pct"], calendar_day=age,
                                 corrected_count=count, fitted_count=bool(age)))
    return observations


MODEL = None


class ScaledObjective:
    def __init__(self, kernel, divisor):
        self.kernel = kernel
        self.supply = kernel.supply
        self.divisor = divisor

    def constrained(self, p, target, conf, gradient=True):
        f, g, margin = self.kernel.constrained(p, target, conf, gradient)
        return f/self.divisor, g/self.divisor if gradient else None, margin


def optimize(k, c, target, start, limit, fixed=None, divisor=10.):
    fit = solve_constrained(ScaledObjective(k, divisor), c, target, start, limit=limit, fixed=fixed,
                            best_iterate_tolerance=1e-6/divisor)
    fit["loss"] *= divisor
    fit["gradient"] *= divisor
    return fit


def model():
    global MODEL
    if MODEL is None:
        MODEL = CountKernel()
    return MODEL


def audit():
    k = model()
    c = Configuration(k)
    assert len(k.y) == 221 and k.ols_n == 114
    slopes = k.slopes(k.observed_logs)
    discrepancy = float(np.max(np.abs(slopes-k.ols_y)))
    assert discrepancy < 1e-12
    x = np.clip(c.x_anchor, .015, .985)
    p,jac = c.unmap(x)
    log_pred,diagnostics = k.log_predict(p)
    assert diagnostics[0] == 0
    assert np.max(abs(log_pred[~k.capture]-k.observed_logs[~k.capture])) < 1e-13
    gradient = c.gradient(k.constrained(p,k.y,c)[1],jac)
    numerical=[]
    distance = 1e-6
    for j in range(len(x)):
        direction = np.eye(len(x))[j] * distance
        values = [k.constrained(c.unmap(x + sign*direction)[0],k.y,c,False)[0] for sign in (-1,1)]
        numerical.append((values[1]-values[0])/(2*distance))
    scale = np.maximum(1.,abs(gradient))
    gradient_error = np.max(np.abs(gradient-numerical)/scale)
    assert gradient_error < 5e-5, (gradient_error,gradient,numerical)
    f,_,margin=k.constrained(p,k.y,c,False)
    assert np.isclose(f,0.5*np.sum((log_pred[k.capture]-k.y)**2),rtol=1e-13)
    write(HERE/"count_inputs.tsv", k.observations)
    write(HERE/"input_audit.tsv", [dict(records=len(k.slots),seed_records=k.ols_n,fitted_records=k.n,
                                       reproduced_OLS_max_abs_error=discrepancy,
                                       scaled_gradient_max_error=gradient_error,
                                       reference_donor_margin=margin)])
    sources = (SNAPSHOT, MAPPING, STEP2/"ols_calendar_map.tsv", STEP2/"experimental_ols_reproduction.tsv",
               OLS/"main_constrained_best.tsv", OLS/"parameter_practical_verdicts.tsv",
               OLS/"qualified_recovery.tsv", OLS/"synthetic_truth.tsv", OLS/"truth_design.tsv",
               OLS/"noisy_slope_inputs.tsv", OLS/"ols_kernel.cpp", OLS/"practical.py")
    write(HERE/"source_provenance.tsv", [dict(file=str(p),sha256=sha256(p.read_bytes()).hexdigest())
                                         for p in sources])
    print("Audit",len(k.slots),len(k.y),discrepancy,gradient_error,flush=True)


if __name__ == "__main__":
    audit()
