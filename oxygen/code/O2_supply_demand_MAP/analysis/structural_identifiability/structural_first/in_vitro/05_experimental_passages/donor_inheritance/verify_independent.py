"""Independent event recurrence, separate dead stocks, complex-step sensitivities."""
import argparse
from collections import defaultdict
from fractions import Fraction

import numpy as np
from flint import arb, ctx

import inherit as study
from algebra import decode
from protocol import INITIALS, MECHANISMS, read, table

HERE = study.HERE
GRID = np.array(study.cert.model.GRID)
Z = (GRID - 88) / 44
N = len(GRID)
H = 1e-24


def group(rows, key):
    d = defaultdict(list)
    for r in rows:
        d[r[key]].append(r)
    return d


def make_operator(oxygen, p):
    # Independent implementation of the symmetric ordinary and WGD event routing.
    weight = p["O2_crit"]**p["n_O"] / (p["O2_crit"]**p["n_O"] + oxygen**p["n_O"]) if oxygen else 1.
    u = GRID / 44
    mu = p["mu_hp"] * weight * u**p["gamma_mu"]
    lam = p["lam_max"] / (1 + p["alpha_o2"] * weight * u**p["gamma_growth"])
    prob = p["p_mis_base"] + p["p_misseg"] * mu/(p["k_o_mis"]+mu)
    survive = p["buffer_smax"] * np.exp(-p["buffer_beta"] * u**(-p["buffer_n_exp"]))
    wgd = p["p_wgd"]
    live = np.zeros((N, N), dtype=np.complex128)
    b = np.zeros(N, dtype=np.complex128)
    for i, n in enumerate(GRID):
        n = int(n)
        ordinary = lam[i] * (1-wgd)
        b[i] = 2 * ordinary * (1 - (1 - prob[i] + prob[i] * survive[i])**n)
        live[i,i] -= lam[i]
        if 2*n <= GRID[-1]:
            live[2*n-GRID[0],i] += lam[i]*wgd
        m = (1-prob[i])**n
        for j in range(n+1):
            if m.real >= 1e-8:
                add = ordinary * m * survive[i]**j
                if n-j >= GRID[0]:
                    live[i-j,i] += add
                if n+j <= GRID[-1]:
                    live[i+j,i] += add
            m *= (n-j)/(j+1) * prob[i]/(1-prob[i])
    return live, mu, b


def replay(root_kind, cohort, base_params, arcs, design, sampled, perturbation=None, reusable=None):
    p = dict(base_params)
    if perturbation in p:
        p[perturbation] += 1j*H
    eta_h = 1/3 + (1j*H if perturbation == f"eta_H_{cohort}" else 0)
    eta_b = 2/3 + (1j*H if perturbation == f"eta_M_{cohort}" else 0)
    out = HERE / root_kind
    if reusable is None or perturbation in MECHANISMS:
        oxygen_levels = {x["oxygen_pct"] for x in design.values()}
        ops = {o: make_operator(float(o), p) for o in oxygen_levels}
    else:
        ops = reusable
    hist_id = study.HIST_IDS[cohort]
    if root_kind == "hist":
        frequencies = {int(r["chromosome_count"]): int(r["cells"])/int(r["total_cells"])
                       for r in read(out/"dp_proxy_histograms.tsv") if r["archival_donor_passage_id"] == hist_id}
        root = np.array([frequencies[int(k)] for k in GRID])
    else:
        root = np.exp(-0.5*((GRID-p[f"init_mean_{cohort}"])/p[f"init_sd_{cohort}"])**2)
        root /= root.sum()
    counts, ols, donor_states = {}, {}, {}
    for passage_id, row in design.items():
        if row["cohort"] != cohort:
            continue
        if row["parent_id"]:
            x = donor_states[row["parent_id"], int(row["donor_step"])]
            pi = x/x.sum()
        else:
            pi = root
        if perturbation == f"tilt_mean:{passage_id}":
            pi = pi * np.exp(1j*H*Z)
            pi /= pi.sum()
        if perturbation == f"tilt_dispersion:{passage_id}":
            pi = pi * np.exp(1j*H*Z**2)
            pi /= pi.sum()
        x = pi.astype(np.complex128) * int(row["fixed_live_cells"])/study.cert.CAPACITY
        direct = mitotic = 0.
        local = sampled.get(passage_id, ())
        times = {int(row["euler_step"]): row for row in local}
        donor_days = arcs.get(passage_id, ())
        need = set(times) | set(donor_days)
        g, mu, b = ops[row["oxygen_pct"]]
        slope = 0.
        for tick in range(max(need)+1):
            L = x.sum()
            if tick in donor_days:
                donor_states[passage_id, tick] = x
            if tick in times:
                slot = times[tick]
                total = L + eta_h*direct + eta_b*mitotic
                counts[slot["record_id"]] = total
                slope += np.log(total * study.cert.CAPACITY) * (
                    int(slot["ols_weight_numerator"])/int(slot["ols_weight_denominator"]))
            crowd = 1-L
            direct += 0.05*(mu @ x - p["k_clear"]*direct)
            mitotic += 0.05*(crowd*(b @ x) - p["k_clear"]*mitotic)
            x = x + 0.05*(crowd*(g @ x) - mu*x)
        if local:
            ols[passage_id] = slope
    return counts, ols, ops


def run(root_kind=None, selected_cohort=None):
    ctx.prec = 256
    sampled = group(read(study.cert.STEP2 / "ols_calendar_map.tsv"), "passage_id")
    data = study.cert.PRODUCTION/"interior_384/exact_parameters.tsv"
    names = tuple(MECHANISMS) + INITIALS
    p = {r["parameter"]: float(Fraction(r["exact_value"])) for r in read(data) if r["parameter"] in names}
    p["k_clear"] = 1e-4
    verification = []
    for kind in ((root_kind,) if root_kind else ("hist", "gaussian")):
        out = HERE / kind
        design = {r["passage_id"]: r for r in read(out/"passage_design.tsv")}
        arcs = defaultdict(set)
        for row in design.values():
            if row["parent_id"]:
                arcs[row["parent_id"]].add(int(row["donor_step"]))
        for cohort in ((selected_cohort,) if selected_cohort else ("2N", "4N")):
            counts, _, ops = replay(kind, cohort, p, arcs, design, sampled)
            interval = {r["record_id"]: float(arb(r["expected_visible_witness_cells"]).mid())
                        for r in read(out/f"{cohort}_counts_replay.tsv")}
            assert counts.keys() == interval.keys()
            forward_error = max(abs(v.real*study.cert.CAPACITY/interval[k] - 1) for k,v in counts.items())
            assert forward_error < 2e-13, forward_error
            gradient = group(read(out/f"{cohort}_ols_derivatives.tsv"), "column")
            for k, values in gradient.items():
                _, slopes, _ = replay(kind, cohort, p, arcs, design, sampled, k, ops)
                saved = {r["output_id"]: float(decode(r).mid()) for r in values}
                assert slopes.keys() == saved.keys()
                discrepancy = max(abs(slopes[q].imag/H - v) for q,v in saved.items())
                scale = max(abs(v) for v in saved.values())
                scaled_error = discrepancy / scale if scale > 0 else discrepancy
                assert scaled_error < 1e-7, (kind, cohort, k, scaled_error)
                verification.append(dict(dp_root=kind, cohort=cohort, column=k,
                                         max_count_relative_error=forward_error,
                                         derivative_max_error_over_max_column=float(scaled_error)))
            print(kind, cohort, len(gradient), "derivative columns checked", flush=True)
    destination = HERE / "verification.tsv"
    if root_kind and selected_cohort:
        destination = HERE / root_kind / f"{selected_cohort}_independent_verification.tsv"
    table(destination, verification)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", choices=("hist", "gaussian"))
    parser.add_argument("--cohort", choices=("2N", "4N"))
    args = parser.parse_args()
    if bool(args.root) != bool(args.cohort):
        parser.error("--root and --cohort must be given together")
    run(args.root, args.cohort)
