"""Matched AR(1) fits with an optional likelihood for two day-three karyotypes."""
import ctypes as C
import os
from pathlib import Path
import sys

import numpy as np
from scipy.special import xlogy

HERE = Path(__file__).resolve().parent
BASE = HERE.with_name("20261118_experimental_counts_ar1_practical_01")
ALIGNMENT = HERE.with_name("20261119_experimental_karyotype_alignment_01")
sys.path.insert(0, str(BASE))
os.environ["CULTURE_AR1_LIB"] = os.environ.get(
    "CULTURE_JOINT_LIB", "/private/tmp/culture_count_karyotype_20261119.dylib"
)
from ar1_model import (CalendarAR1, Config, FD, FI, RHO_MAX, STATIONARY, THRESHOLD,
                       likelihood_ratio, pool, read, row_x, solve, write)

TARGETS = ("p_misseg", "buffer_beta", "k_o_mis", "buffer_smax", "buffer_n_exp")
GRID = np.arange(22, 155)
OBSERVATION_EPSILON = 1e-12


class Experiment(CalendarAR1):
    def __init__(self, convention="autosomes", width=0.):
        super().__init__()
        self.convention = convention
        self.width = width
        lib = self.lib
        lib.joint_set_histograms.argtypes = [C.c_void_p, FD]
        lib.joint_live_states.argtypes = [C.c_void_p, FD, C.c_int, FI, FD]
        lib.joint_state_adjoint.argtypes = [C.c_void_p, FD, C.c_int, FI, FD, FD]
        field = {"autosomes": "observed_N_22_autosomes",
                 "archived": "observed_N_archived"}[convention]
        cells = read(ALIGNMENT / "cell_alignment.tsv")
        samples = read(ALIGNMENT / "sample_alignment.tsv")
        direct = [r for r in samples if r["role"] == "direct_harvest_anchor"]
        assert len(direct) == 2
        indices = {r["passage_id"]: j for j, r in enumerate(self.graph)}
        pairs = []
        for r in direct:
            day = int(r["harvest_anchored_assay_calendar_day"])
            slot, = [s for s in self.slots if s["record_id"] == r["assay_origin_id"]]
            assert slot["passage_id"] == r["physical_passage_id"]
            assert int(slot["calendar_day"]) == day == 3
            assert slot["cohort"] == "4N" and float(slot["oxygen_pct"]) == 0.
            pairs.append((indices[slot["passage_id"]], int(slot["euler_step"])))
        self.pairs = np.ascontiguousarray(pairs, dtype=np.int32)
        self.origin_ids = [r["assay_origin_id"] for r in direct]

        def hist(origin):
            n = [int(r[field]) for r in cells if r["assay_origin_id"] == origin]
            assert len(n) == 20 and min(n) >= GRID[0] and max(n) <= GRID[-1]
            return np.bincount(np.array(n)-GRID[0], minlength=len(GRID)).astype(float)

        proxies = sorted((r for r in samples if r["role"] == "proxy_only"),
                         key=lambda r: r["cohort"])
        assert [r["cohort"] for r in proxies] == ["2N", "4N"]
        self.proxies = np.ascontiguousarray([hist(r["assay_origin_id"])/20 for r in proxies])
        assert lib.joint_set_histograms(self.handle, self.proxies)
        self.kary = np.ascontiguousarray([hist(r["assay_origin_id"]) for r in direct])
        if width:
            H = np.exp(-.5*((GRID[:, None]-GRID)/width)**2)
            self.smoothing = H/H.sum(axis=0, keepdims=True)
        else:
            self.smoothing = None

    def prediction(self, p):
        state = np.empty_like(self.kary)
        if not self.lib.joint_live_states(self.handle, np.ascontiguousarray(p), 2, self.pairs, state):
            raise ArithmeticError("Could not extract a viable population at the aligned harvest")
        prob = state/state.sum(axis=1, keepdims=True)
        prob = prob if self.smoothing is None else prob @ self.smoothing.T
        prob = (1.-OBSERVATION_EPSILON)*prob+OBSERVATION_EPSILON/len(GRID)
        return prob, state

    def kary_loss(self, p, hist, gradient=False):
        prob, state = self.prediction(p)
        loss = -float(xlogy(hist, prob).sum())
        if not gradient:
            return loss, None
        partial = -(1.-OBSERVATION_EPSILON)*hist/prob
        if self.smoothing is not None:
            partial = partial @ self.smoothing
        mass = state.sum(axis=1, keepdims=True)
        g = np.ascontiguousarray((partial-(partial*state).sum(axis=1, keepdims=True)/mass)/mass)
        adjoint = np.zeros(47)
        if not self.lib.joint_state_adjoint(self.handle, p, 2, self.pairs, g, adjoint[:46]):
            raise ArithmeticError("Could not differentiate the aligned harvest probabilities")
        return loss, adjoint

    def target(self, log_counts=None, kary=None):
        count = self.y if log_counts is None else log_counts
        return np.ascontiguousarray(count if kary is None else np.r_[count, kary.ravel()])

    def evaluate_joint(self, p, target, gradient=False):
        count, g_count = super().metrics(p, target[:self.n], gradient)
        if len(target) == self.n or count[0] >= 1e12:
            return count, g_count, count[0], 0.
        hist = target[self.n:].reshape(self.kary.shape)
        neglog, g_kary = self.kary_loss(p, hist, gradient)
        # Subtract the saturated multinomial term for a stable objective scale.
        # Its constant cancels in all LR comparisons on the same dataset.
        saturated = -xlogy(hist, hist/hist.sum(axis=1, keepdims=True)).sum()
        multiplier = np.exp(2*(neglog-saturated)/self.n)
        loss = count[0]*multiplier
        g = None if not gradient else multiplier*(g_count+2.*count[0]/self.n*g_kary)
        return count, g, loss, neglog

    def constrained(self, p, target, conf, gradient=True):
        out, grad, loss, _ = self.evaluate_joint(p, target, gradient)
        return loss, grad, out[1]


def fit_receipt(k, c, y, fit, **labels):
    p, _ = c.unmap(fit["x"])
    out, _, loss, kary_nll = k.evaluate_joint(p, y)
    return dict(labels, loss=float(loss), count_loss=float(out[0]), kary_nll=kary_nll,
                donor_margin=float(out[1]), stationarity=fit["gradient"], rho=float(p[-1]),
                sigma=float(np.sqrt(out[2]/k.n)), success=fit["success"],
                termination=fit["termination"], seconds=fit["seconds"],
                **{f"x_{n}": float(v) for n,v in zip(c.names,fit["x"])},
                **{n: float(p[j]) for n,j in zip(c.names[:-1],c.active[:-1])})


def stable(r):
    return float(r["donor_margin"]) >= 0 and float(r["stationarity"]) <= STATIONARY


def audit():
    rows = []
    baseline = read(BASE/"observed_best.tsv")[0]
    for case, width in (("archived",0.), ("autosomes",0.), ("autosomes",1.)):
        k = Experiment(case, width)
        c = Config(k)
        x = row_x(c,baseline)
        p, _ = c.unmap(x)
        y = k.target(kary=k.kary)
        out, _, f, nk = k.evaluate_joint(p,y)
        prob, _ = k.prediction(p)
        assert np.allclose(prob.sum(axis=1), 1.)
        assert abs(nk+np.sum(xlogy(k.kary,prob))) < 1e-11
        if case == "archived":
            assert abs(out[0]-float(baseline["loss"])) < 1e-11
        # Check the full likelihood and numerical derivatives at an interior
        # point, and confirm that a state adjoint did not alter the cached fit.
        x = np.clip(x,.002,.998)
        p,j = c.unmap(x)
        count,g,f_joint,nk = k.evaluate_joint(p,y,True)
        f_check = k.evaluate_joint(p,y)[2]
        assert abs(f_check-f_joint) < 1e-12
        expected = float(k.n*np.log(count[0]/out[0]) +
                         2*(nk-k.evaluate_joint(c.unmap(row_x(c,baseline))[0],y)[3]))
        assert abs(likelihood_ratio(f_joint, f, k.n)-expected) < 1e-10
        h = 2e-6
        numeric = np.array([(k.evaluate_joint(c.unmap(x+h*v)[0],y)[2]
                            - k.evaluate_joint(c.unmap(x-h*v)[0],y)[2])/(2*h)
                           for v in np.eye(len(x))])
        analytic = c.gradient(g,j)
        error = max(abs(analytic-numeric)/np.maximum(1.,abs(numeric)))
        assert error < 3e-4, (error, numeric, analytic)
        rows.append(dict(convention=case,chromosome_sd=width,cells=int(k.kary.sum()),
                         count_at_stage2_vector=out[0],joint_at_stage2_vector=f,
                         max_scaled_gradient_error=error))
    write(HERE/"likelihood_validation.tsv",rows)
    print(rows,flush=True)


if __name__ == "__main__":
    audit()
