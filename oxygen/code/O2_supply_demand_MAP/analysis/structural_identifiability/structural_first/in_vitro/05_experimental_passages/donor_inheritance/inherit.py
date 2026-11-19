"""Conditioned inheritance for the experimental calendar-day culture design."""
import argparse
from collections import Counter, defaultdict
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

from flint import arb, arb_mat, ctx
import numpy as np

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "20261118_culture_experimental_OLS_rank_01"
sys.path.insert(0, str(PREVIOUS))
import certify_experimental as cert
from algebra import exact, encode, initial
from protocol import INITIALS, MECHANISMS, read, table

N = cert.N
DT = cert.DT
ROOT_IDS = {f"SUM-159_NLS_{cohort}_dp_seed": cohort for cohort in ("2N", "4N")}
HIST_IDS = {"2N": "SUM-159_NLS_2N_A5_seed", "4N": "SUM-159_NLS_4N_A2_seed"}


def build_design(out):
    with TemporaryDirectory() as scratch:
        original, _, _ = cert.load_design(Path(scratch))
    path = cert.RESULTS / "20261118_culture_experimental_fixed_inoculum_01/ols_shared/cut_edges.tsv"
    cut = {r["child_passage_id"]: r for r in read(path)}
    psgs = dict(original)
    links = []
    for child, row in cut.items():
        if child in ROOT_IDS:
            links.append(dict(row, new_policy="independent_dp_proxy", donor_step=""))
            continue
        parent = row["donor_passage_id_in_crosswalk"]
        donors = [t for t, record, _ in psgs[parent].samples if record == row["donor_record_id"]]
        cert.need(len(donors) == 1, f"unmatched donor {child}")
        psgs[child] = replace(psgs[child], parent=parent, donor_tick=donors[0])
        links.append(dict(row, new_policy="inherit_recorded_donor_composition_no_gap_integration", donor_step=donors[0]))
    remaining = dict(psgs)
    order = {}
    while remaining:
        ready = [p for p in remaining.values() if p.parent is None or p.parent in order]
        cert.need(bool(ready), "ancestry cycle")
        for p in ready:
            order[p.id] = remaining.pop(p.id)
    demand = defaultdict(Counter)
    for child in order.values():
        if child.parent:
            demand[child.parent][child.donor_tick] += child.seed_count
    cert.need(len(order) == 116 and sum(len(p.samples) for p in order.values()) == 335 and
              {p.id for p in order.values() if not p.parent} == set(ROOT_IDS), "wrong graph")
    table(out / "restart_policy.tsv", links)
    table(out / "passage_design.tsv", [
        dict(passage_id=p.id, cohort=p.cohort, parent_id=p.parent or "",
             donor_step=p.donor_tick, fixed_live_cells=p.seed_count, oxygen_pct=p.oxygen,
             observation_count=len(p.samples)) for p in order.values()
    ])
    return order, demand, set(cut) - set(ROOT_IDS)


def points(out):
    p = cert.parameter_point(out)
    rows = read(cert.PRODUCTION / "interior_384/exact_parameters.tsv")
    for r in rows:
        if r["parameter"] in INITIALS:
            p[r["parameter"]] = exact(r["exact_value"])
    return p


def root_histograms(out):
    raw = read(cert.ARCHIVE / "invitro_observed_kary.tsv")
    hists = {}
    rows = []
    for cohort, donor in HIST_IDS.items():
        counts = Counter(int(r["observed_kary_N"]) for r in raw if r["passage_id"] == donor)
        cert.need(sum(counts.values()) == 20 and set(counts) <= set(cert.model.GRID), f"unsupported proxy {donor}")
        hists[cohort] = arb_mat([[exact(Fraction(counts[k], 20))] for k in cert.model.GRID])
        for k in cert.model.GRID:
            rows.append(dict(cohort=cohort, archival_donor_passage_id=donor, chromosome_count=k,
                             cells=counts[k], total_cells=20))
    table(out / "dp_proxy_histograms.tsv", rows)
    return hists


def operator(p, oxygen):
    cert.model.rates = cert.fixed_rates
    ox = exact(oxygen)
    op = cert.model.operator(p, ox)
    m = len(MECHANISMS)
    dt = exact(DT)
    mu = op["mu"]
    bio, dbio = [], []
    for n, v, ledger in zip(cert.model.GRID, op["values"], op["ledger"]):
        (lam, _, prob, logs), (dl, _, dp, ds) = cert.fixed_rates(p, arb(n) / 44, ox)
        s = logs.exp()
        base = 1 - prob * (1 - s)
        f = 1 - base**n
        w = p["p_wgd"]
        rate = 2 * lam * (1 - w) * f
        cert.need(rate > 0 and (rate - ledger["biological"]).contains(0), "biological dead stock mismatch")
        bio.append(rate)
        dbio.append([
            2 * ((dl[k] * (1 - w) - lam * int(k == "p_wgd")) * f
                 - lam * (1 - w) * n * base**(n-1) * (-(1-s) * dp[k] + prob * s * ds[k]))
            for k in MECHANISMS
        ])
    diagonal = [1 - dt * mu[j] for j in range(N)]
    cert.need(all(d - dt * v["lam"] > 0 for d, v in zip(diagonal, op["values"])), "clamp branch")
    g = op["G"] * dt
    derivatives = arb_mat([[op["dG"][j][i, k] * dt for k in range(N)]
                           for i in range(N) for j in range(m)])
    return (g, derivatives, diagonal, arb_mat([mu]) * dt,
            arb_mat([[x * dt for x in v] for v in op["dmu"]]),
            arb_mat([bio]) * dt, arb_mat(dbio).transpose() * dt,
            op["mask_sha256"])


def columns_for(cohort, ordered, restored, root_kind):
    columns = list(MECHANISMS) + [f"eta_H_{cohort}", f"eta_M_{cohort}"]
    if root_kind == "gaussian":
        columns += [f"init_mean_{cohort}", f"init_sd_{cohort}"]
    for p in ordered.values():
        if p.cohort == cohort and p.id in restored:
            columns += [f"tilt_mean:{p.id}", f"tilt_dispersion:{p.id}"]
    return columns


def normalized(x, j, L):
    pi = x / L
    size_j = cert.sums(j)
    return pi, arb_mat([[(j[i, k] - pi[i, 0] * size_j[k]) / L for k in range(j.ncols())] for i in range(N)])


def add_tilt(pi, j, columns, passage_id):
    z = [exact(Fraction(n - 88, 44)) for n in cert.model.GRID]
    for prefix, design in (("tilt_mean", z), ("tilt_dispersion", [v*v for v in z])):
        idx = columns.index(f"{prefix}:{passage_id}")
        mean = sum(pi[i, 0] * design[i] for i in range(N))
        for i in range(N):
            j[i, idx] += pi[i, 0] * (design[i] - mean)


def propagate(x, j, wanted, operator_data, decay, eta_h, eta_m):
    G, derivative, diagonal, H, dH, B, dB, _ = operator_data
    m, dimension = len(MECHANISMS), j.ncols()
    # Passive weighted dead stock. Its derivative includes both detection efficiencies.
    dead, dj = arb(0), [arb(0) for _ in range(dimension)]
    min_margin = arb(1)
    result = {}
    final = max(wanted)
    for tick in range(final + 1):
        L = sum(x[i, 0] for i in range(N))
        crowd = 1 - L
        cert.need(L * cert.CAPACITY > exact("1e-12") and crowd > 0, f"Euler branch at {tick}")
        min_margin = min(min_margin, crowd.lower())
        if tick in wanted:
            c = L + dead
            cert.need(c * cert.CAPACITY > 1, f"OLS log floor at {tick}")
            result[tick] = x, j, L, c, [v + w for v, w in zip(cert.sums(j), dj)]
        if tick == final:
            break
        gx, gj = G * x, G * j
        gd, hd, bd = derivative * x, dH * x, dB * x
        live_j = cert.sums(j)
        h, b = (H * x)[0, 0], (B * x)[0, 0]
        hj, bj = H * j, B * j
        dead = decay * dead + eta_h * h + eta_m * crowd * b
        dj = [decay * dj[k] + eta_h * hj[0, k] + eta_m * (crowd * bj[0, k] - live_j[k] * b)
              for k in range(dimension)]
        for k in range(m):
            dj[k] += eta_h * hd[k, 0] + eta_m * crowd * bd[k, 0]
        dj[m] += h
        dj[m + 1] += crowd * b
        j = arb_mat([
            [diagonal[i] * j[i, k] + crowd * gj[i, k] - gx[i, 0] * live_j[k]
             + (crowd * gd[i*m + k, 0] - x[i, 0] * dH[k, i] if k < m else 0)
             for k in range(dimension)] for i in range(N)
        ])
        x = arb_mat([[diagonal[i] * x[i, 0] + crowd * gx[i, 0]] for i in range(N)])
    return result, min_margin


def local_rows(ordered, demand, restored, histograms, p, ops, root_kind, cohort, out):
    cols = columns_for(cohort, ordered, restored, root_kind)
    dim = len(cols)
    table(out / f"{cohort}_columns.tsv", [dict(column=k) for k in cols])
    states = {}
    proof = []
    rows = {"ols": {}, "counts": {}}
    samples = []
    for psg in ordered.values():
        if psg.cohort != cohort:
            continue
        if psg.parent:
            x, j, L, _, _ = states[psg.parent][psg.donor_tick]
            pi, pi_j = normalized(x, j, L)
            if psg.id in restored:
                add_tilt(pi, pi_j, cols, psg.id)
        elif root_kind == "hist":
            pi, pi_j = histograms[cohort], arb_mat(N, dim)
        else:
            pi, pi_j = initial(p, cohort, cols, N0=1)
        a = exact(Fraction(psg.seed_count, cert.CAPACITY))
        ticks = {t for t, _, _ in psg.samples} | set(demand.get(psg.id, ()))
        try:
            values, margin = propagate(pi * a, pi_j * a, ticks, ops[psg.oxygen],
                                       1 - exact(DT) * p["k_clear"], exact("1/3"), exact("2/3"))
        except ValueError as exc:
            raise ValueError(f"{psg.id}: {exc}") from exc
        for tick, requested in demand.get(psg.id, {}).items():
            supply = values[tick][2] * cert.CAPACITY - requested
            cert.need(supply > 0, f"unmet fixed combined inoculum at {psg.id}, day {tick * DT}: {supply}")
            proof.append(dict(passage_id=psg.id, day=str(tick*DT), fixed_child_demand=requested,
                              lower_live_surplus_cells=str(supply.lower()), min_logistic_gap=str(margin)))
        states[psg.id] = {t: values[t] for t in demand.get(psg.id, ())}
        if psg.samples:
            t0 = values[0]
            cert.need((t0[3]-a).contains(0) and all(v.contains(0) for v in t0[4]), "fixed seed invariant")
            rows["ols"][psg.id] = [sum(exact(w) * values[t][4][k] / values[t][3]
                                       for t, _, w in psg.samples[1:]) for k in range(dim)]
            for tick, key, _ in psg.samples:
                samples.append(dict(passage_id=psg.id, record_id=key, euler_step=tick,
                                    expected_visible_witness_cells=str(values[tick][3] * cert.CAPACITY)))
                if tick:
                    # Scaled counts have exactly the same rank as cell counts.
                    rows["counts"][key] = values[tick][4]
    table(out / f"{cohort}_supply.tsv", proof)
    table(out / f"{cohort}_counts_replay.tsv", samples)
    for phase, result in rows.items():
        # Save one derivative run for later phase gating and visibility configurations.
        table(out / f"{cohort}_{phase}_derivatives.tsv",
              (dict(output_id=name, column=k, **encode(col))
               for name, r in result.items() for k, col in zip(cols, r)))
    return rows


def run(root_kind, bits, out):
    out.mkdir(exist_ok=True, parents=True)
    ctx.prec = bits
    ordered, demand, restored = build_design(out)
    p = points(out)
    hist = root_histograms(out)
    oxygens = sorted({r.oxygen for r in ordered.values()}, key=Fraction)
    ops = {o: operator(p, o) for o in oxygens}
    saved_masks = {r["oxygen_pct"]: r["mask_sha256"] for r in read(cert.PRODUCTION / "interior_384/operator_branches.tsv")}
    cert.need(all(op[-1] == saved_masks[o] for o, op in ops.items()), "changed live model operator mask")
    for cohort in ("2N", "4N"):
        local_rows(ordered, demand, restored, hist, p, ops, root_kind, cohort, out)
        print("Derivative propagation complete", cohort, root_kind, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, choices=("hist", "gaussian"))
    parser.add_argument("--bits", type=int, default=384)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.root, args.bits, args.out)
