"""Independent complex-step observer validation and archived minor recheck."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import sys
from time import process_time

import numpy as np
from flint import arb_mat, ctx

import joint_rank as rank
from algebra import decode, independent_minor_check
from verify_independent import make_operator
from paths import FROZEN_DESIGN, PROFILES, output_directory

H = 1e-24
EPS = 1e-12
GRID = np.array(rank.CHROMOSOMES)


def forward(case, raw, root_set=None, perturbation=None, cache=None, design_dir=FROZEN_DESIGN):
    p = {k: complex(float(Fraction(v))) for k,v in raw.items()}
    if perturbation and len(perturbation)==1:
        p[perturbation[0]] += 1j*H
    graph, obs, demand, wanted = rank.design(case, design_dir)
    proxy = rank.proxy_states()
    rates, state, values, margins = {}, {}, {}, []
    epoch_counts = defaultdict(set)
    for r in graph:
        if r["parent"]:
            epoch_counts[r["parent"]].add(r["donor_tick"])
    for row in graph:
        node, parent, cohort, root = (row[k] for k in ("passage_id","parent","cohort","root"))
        if root_set and root not in root_set:
            continue
        oxygen = float(row["oxygen_pct"])
        if oxygen not in rates:
            if cache is not None and (not perturbation or perturbation[-1] not in rank.MECHANISMS):
                if oxygen not in cache:
                    cache[oxygen] = make_operator(oxygen,p)
                rates[oxygen] = cache[oxygen]
            else:
                rates[oxygen] = make_operator(oxygen,p)
        G,mu,division_death = rates[oxygen]
        if parent:
            donor = state[parent,row["donor_tick"]]
            pi = donor / donor.sum()
        elif case == "gaussian":
            pi = np.exp(-.5*((GRID-p[f"init_mean_{cohort}"])/p[f"init_sd_{cohort}"])**2)
            pi /= pi.sum()
        else:
            root_profile = rank.free_profile(cohort) if case in ("free","elapsed_free","free_seeds") else proxy[cohort]
            pi = np.array([float(q) for q in root_profile],dtype=np.complex128)
            if perturbation and perturbation[0]==root and len(perturbation)==2 and perturbation[1].startswith("pi_"):
                pi[int(perturbation[1][3:])-22] += 1j*H
                pi[-1] -= 1j*H
        if case == "tv" and "_A7K_" in node:
            z = (GRID-(44 if cohort=="2N" else 88))/44
            pi = pi*np.exp(p[f"tilt_mean_{cohort}"]*z + p[f"tilt_shape_{cohort}"]*z*z)
            pi /= pi.sum()
        seed = row["seed"]/rank.CAPACITY
        bh = bm = 0j
        live_seed = seed
        if case == "free_seeds" and row["donor_passage"]:
            live_seed *= .99
            if perturbation == (root,"seed:"+node):
                live_seed += 1j * H * seed
            bh = bm = (seed-live_seed)/2
        live = pi*live_seed
        by_tick = defaultdict(list)
        for r in obs[node]:
            by_tick[r["tick"]].append(r)
        final = max(wanted[node])
        q = p["misseg_dead_share"]
        for tick in range(final+1):
            L = live.sum()
            assert 0 < L.real < 1
            if tick in epoch_counts[node]:
                state[node,tick] = live
            if tick in demand[node]:
                margins.append((node,tick,L.real*rank.CAPACITY-demand[node][tick]))
            for obs_row in by_tick[tick]:
                kind, event = rank.MODALITY[obs_row["modality"]],obs_row["event"]
                if kind=="C":
                    values[kind,event] = np.log(L+bh+bm)
                elif kind=="D":
                    odds = (bh+q*bm+EPS)/(L+(1-q)*bm+EPS)
                    values[kind,event] = p["class_log_odds"]+np.log(odds)
                else:
                    for n,v in zip(GRID,live/L):
                        values[kind,event,int(n)] = v
            if tick == final:
                break
            dt = float(row["dt_last"] if row["final"] is not None and tick+1 == final else rank.DT)
            crowd = 1-L
            bh = (1-dt*1e-4)*bh + dt*p["eta_H"]*(mu@live)
            bm = (1-dt*1e-4)*bm + dt*crowd*p["eta_M"]*(division_death@live)
            live = (1-dt*mu)*live + dt*crowd*(G@live)
    return values,margins


def physical(key,value):
    if key[0]=="C":
        return float(np.exp(value.real)*rank.CAPACITY)
    if key[0]=="D":
        return float(EPS+(1-2*EPS)/(1+np.exp(-value.real)))
    return float((1-EPS)*value.real+EPS/rank.N)


def compare(case, raw, roots, blocks, prediction, design_dir, all_local=False):
    start = process_time()
    reference,_ = forward(case,raw,roots,design_dir=design_dir)
    by_type = defaultdict(list)
    for key,value in prediction.items():
        expected = float(value)
        actual = physical(key,reference[key])
        difference = abs(actual / expected - 1) if key[0]=="C" else abs(actual-expected)
        by_type[key[0]].append(difference)
    checks = [dict(check="observer_value",coordinate=k,n=len(v),error=max(v),threshold=1e-10,passed=max(v)<1e-10)
              for k,v in sorted(by_type.items())]
    cache = {}
    for root,(cols,J) in blocks.items():
        selected = [j for j,k in enumerate(cols) if all_local or not k.startswith("pi_") or k in ("pi_22","pi_46","pi_88","pi_153")]
        for j in selected:
            name = cols[j]
            is_local = name.startswith(("pi_","seed:"))
            step = (root,name) if is_local else (name,)
            vals,_ = forward(case,raw,{root},step,cache,design_dir=design_dir)
            for kind in ("C","D","K"):
                keys = [k for k in J if k[0]==kind]
                if not keys:
                    continue
                analytic_values = np.array([float(J[k][j].mid()) for k in keys])
                actual = np.array([vals[k].imag/H for k in keys])
                scale = max(float(np.max(np.abs(analytic_values))),1e-10)
                error = float(np.max(np.abs(analytic_values-actual)))/scale
                checks.append(dict(check="derivative_"+kind,coordinate=root+"|"+name,n=len(keys),
                                   error=error,threshold=1e-8,passed=error<1e-8))
    print(f"Independent observer validation took {process_time()-start:.3f} CPU seconds",flush=True)
    return checks


def archived_controls(out):
    ctx.prec = 768
    fixed = out/"archived_free_control"
    result = defaultdict(list)
    for row in rank.read(fixed/"root_minors.tsv"):
        result[row["root_id"]].append(row)
    results = []
    for name,rows in result.items():
        row_ids = list(dict.fromkeys(r["output_id"] for r in rows))
        col_ids = list(dict.fromkeys(r["chromosome"] for r in rows))
        M = arb_mat([[decode(r) for r in rows[i*len(col_ids):(i+1)*len(col_ids)]] for i in range(len(row_ids))])
        passed,bound = independent_minor_check(M)
        results.append(dict(source="replayed_free_count",minor=name,order=M.nrows(),check="exact_rational_preconditioner",
                            bound=str(bound),passed=passed))
    assert sum(r["order"] for r in results)==221
    strict = rank.STRUCT/"20261118_culture_experimental_counts_local_01/certified_minors.tsv"
    rows = [r for r in rank.read(strict) if r["case"]=="original_interior"]
    M = arb_mat([[decode(r) for r in rows[i*16:(i+1)*16]] for i in range(16)])
    passed,bound = independent_minor_check(M)
    results.append(dict(source="archived_inherited_counts",minor="original_interior",order=16,
                        check="exact_rational_preconditioner",bound=str(bound),passed=passed))
    return results


def controls(design_dir):
    ctx.prec = 384
    out = []

    def record(name, passed, detail):
        out.append(dict(test=name,passed=bool(passed),detail=detail))

    _, bound = rank.precondition(arb_mat([[1,2],[3,4]]))
    record("nonsingular_matrix",bound<1,f"Outward inverse residual {float(bound):.3g}")
    _, bound = rank.precondition(arb_mat([[1,2],[1,2]]))
    record("duplicated_output_rows",not bound<1,"Duplicate row has no nonsingular two-row minor")
    _, bound = rank.precondition(arb_mat([[1,0],[2,0]]))
    record("disabled_measurement_column",not bound<1,"Zero observer column cannot add rank")
    raw = rank.witness("inherit","original")
    numeric = {k:float(Fraction(v)) for k,v in raw.items()}
    for zero, carrier, changed in (("p_misseg","k_o_mis",.007),("buffer_beta","buffer_n_exp",7.)):
        first = dict(numeric)
        first[zero] = 0.
        second = dict(first)
        second[carrier] = changed
        oxygen = sorted({float(r["oxygen_pct"]) for r in rank.read(design_dir/"passage_graph.tsv")})
        agrees = all(np.array_equal(a,b) for o in oxygen for a,b in zip(make_operator(o,first),make_operator(o,second)))
        record(f"exact_zero_{zero}",agrees,f"{carrier} leaves live, hypoxic-death and division-death rates unchanged")
    full, margins = forward("inherit",raw,design_dir=design_dir)
    histograms = defaultdict(list)
    for k,v in full.items():
        if k[0]=="K":
            histograms[k[1]].append(v)
    normalization = max(abs(sum(v)-1.) for v in histograms.values())
    record("karyotype_normalization",len(histograms)==6 and normalization<1e-12,f"Six normalized 133-state distributions, max error {normalization:.3g}")
    _, blocked = forward("inherit",dict(raw,lam_max="0.005"),design_dir=design_dir)
    record("infeasible_donor",all(r[2]>0 for r in margins) and any(r[2]<0 for r in blocked),
           "Historical supply slack changes sign in an intentional low-growth control")
    previous = PROFILES/"count_deletion_leakage_check.tsv"
    row = rank.read(previous)[0]
    record("count_deletion_leakage_archived",all(v=="True" for k,v in row.items() if k!="contract"),
           "Archived zero-loading score and gradient invariant to changes in excluded count values")
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls",action="store_true")
    parser.add_argument("--archive",action="store_true")
    parser.add_argument("--design",type=Path,default=FROZEN_DESIGN)
    parser.add_argument("--out",type=output_directory,required=True)
    args = parser.parse_args()
    if args.controls:
        receipts = archived_controls(args.out)
        rank.write(args.out/"control_replays.tsv",receipts)
        assert all(row["passed"] for row in receipts)
        fixtures = controls(args.design)
        rank.write(args.out/"structural_fixtures.tsv",fixtures)
        assert all(row["passed"] for row in fixtures)
    if args.archive:
        saved = {(rank.MODALITY[r["readout"]],r["event_id"].replace("SUM159_NLS","SUM-159_NLS")):float(r["predicted"])
                 for r in rank.read(rank.FIT/"readout_predictions.tsv") if r["mode"]=="primary"
                 and r["readout"]=="corrected_count"}
        values,margins = forward("inherit",rank.witness("inherit","archived_fit"),design_dir=args.design)
        errors = [abs(physical(k,values[k])/v-1) for k,v in saved.items()]
        assert len(errors)==225 and max(errors)<1e-10
        print(f"Archived count parity: {max(errors):.3g}",flush=True)
