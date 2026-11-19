"""Certify nuisance-only full row rank for the cut-lineage culture design.

The nominal values select a structural witness. They are not parameter estimates.
Rows/columns are reordered by physical ancestry for the determinant factorization.
"""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from fractions import Fraction
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = HERE.parent
STEP1 = RESULTS / "20261118_culture_experimental_OLS_step01"
STEP2 = RESULTS / "20261118_culture_experimental_OLS_step02"
PRODUCTION = RESULTS / "20261117_culture_actual_schedule_step02"
ARCHIVE = RESULTS / "20261109_seed144_wp0/source/seed144"
sys.path.insert(0, str(ROOT / "oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/culture_stageB_finite"))
sys.path.insert(0, str(PRODUCTION))
import model
from algebra import exact, encode, nonzero
from certify import fixed_rates
from protocol import MECHANISMS, read, table

N = len(model.GRID)
CAPACITY = 10**12
DT = Fraction(1, 20)


@dataclass(frozen=True)
class Passage:
    id: str
    cohort: str
    model_id: str
    oxygen: str
    parent: str | None
    donor_tick: int | None
    seed_count: int
    samples: tuple


def need(condition, message):
    if not condition:
        raise ValueError(message)


def timestep(date, start):
    days = (datetime.fromisoformat(date).date() - datetime.fromisoformat(start).date()).days
    need(days >= 0, f"negative calendar age {date}, seed {start}")
    return int(days / DT)


def load_design(out):
    crosswalk = {r["passage_id"]: r for r in read(RESULTS / "20261117_culture_actual_schedule_step01/passage_crosswalk.tsv")}
    followup = read(STEP2 / "passage_gap_followup.tsv")
    raw_counts = {r["record_id"]: Fraction(r["db_correctedCount"]) for r in read(STEP1 / "count_event_mapping.tsv")}
    calendar = defaultdict(list)
    for row in read(STEP2 / "ols_calendar_map.tsv"):
        calendar[row["passage_id"]].append(row)
    for rows in calendar.values():
        rows.sort(key=lambda r: int(r["euler_step"]))
    cuts = {}
    for edge in followup:
        long_gap = edge["recorded_gap_over_7_days_for_review"] == "TRUE"
        note = "maintenance" in (edge["donor_comment_as_recorded"] + " " + edge["seed_comment_as_recorded"]).lower()
        if long_gap or note:
            cuts[edge["child_passage_id"]] = (
                edge | {"cut_due_to_audit_gap": long_gap, "cut_due_to_maintenance_note": note}
            )
    need(len(cuts) == 14 and len(calendar) == 114, "registered design changed")
    table(out / "cut_edges.tsv", cuts.values())

    def retained_parent(passage_id):
        if passage_id in cuts:
            return None
        parent = crosswalk[passage_id]["physical_parent_passage_id"]
        return None if parent in ("", "NA") else parent

    order = []
    visited = set()

    def visit(passage_id):
        if passage_id in visited:
            return
        parent = retained_parent(passage_id)
        if parent is not None:
            visit(parent)
        order.append(passage_id)
        visited.add(passage_id)

    for passage_id in sorted(calendar):
        visit(passage_id)
    need(len(order) == 116 and {p for p in order if retained_parent(p) is None} == set(cuts),
         "expected 114 OLS passages, two retained upstream passages and 14 independent roots")
    passages = {}
    targets = defaultdict(list)
    for passage_id in order:
        row = crosswalk[passage_id]
        parent = retained_parent(passage_id)
        donor_tick = None
        if parent is not None:
            donor = row["direct_parent_record_id"]
            donor_tick = timestep(row["direct_parent_record_time"], crosswalk[parent]["seeding_time"])
            if parent in calendar:
                need(any(r["record_id"] == donor and int(r["euler_step"]) == donor_tick for r in calendar[parent]),
                     f"recorded donor does not match parent OLS calendar {donor}")
            targets[parent].append(donor_tick)
        samples = []
        for sample in calendar.get(passage_id, ()):
            tick = int(sample["euler_step"])
            need(Fraction(tick) * DT == int(sample["calendar_day"]), "noncalendar observation")
            need(sample["model_segment_id"] == row["production_segment_id"], "wrong model segment")
            weight = Fraction(int(sample["ols_weight_numerator"]), int(sample["ols_weight_denominator"]))
            samples.append((tick, sample["record_id"], weight))
        count = Fraction(row["seed_corrected_count"])
        need(count.denominator == 1 and count > 0, f"missing seed count {passage_id}")
        if samples:
            need(samples[0][0] == 0 and raw_counts[samples[0][1]] == count and len(samples) >= 2,
                 f"invalid OLS group {passage_id}")
            need(sum(w for _, _, w in samples) == 0 and sum(t * DT * w for t, _, w in samples) == 1,
                 f"wrong OLS weights {passage_id}")
        passages[passage_id] = Passage(passage_id, row["cohort"], row["production_segment_id"],
                                       row["assigned_oxygen_pct"], parent, donor_tick, int(count), tuple(samples))
    need(sum(len(p.samples) for p in passages.values()) == 335, "missing raw outputs")
    need(len({(p.cohort, p.model_id) for p in passages.values()}) == 116, "active model segment alias")
    roots = {}
    uses = {p: Counter(ticks) for p, ticks in targets.items()}
    for passage in passages.values():
        roots[passage.id] = roots[passage.parent] if passage.parent else passage.id
    table(out / "active_passages.tsv", (
        dict(passage_id=p.id, cohort=p.cohort, model_segment_id=p.model_id, oxygen_pct=p.oxygen,
             root_id=roots[p.id], retained_parent_id=p.parent or "", parent_donor_step=p.donor_tick,
             seed_corrected_count=p.seed_count, sample_steps=";".join(str(s[0]) for s in p.samples))
        for p in passages.values()
    ))
    return passages, roots, uses


def parameter_point(out):
    saved = {r["parameter"]: r["exact_value"] for r in read(PRODUCTION / "interior_384/exact_parameters.tsv")}
    clearance = next(r["value"] for r in read(ARCHIVE / "best_params.tsv") if r["parameter"] == "k_clear")
    saved = {name: saved[name] for name in MECHANISMS} | {"k_clear": clearance}
    table(out / "mechanism_witness.tsv", (dict(parameter=k, exact_value=v) for k, v in saved.items()))
    return {k: exact(v) for k, v in saved.items()}


def operators(p, oxygens, out):
    model.rates = fixed_rates
    result = {}
    ledger = []
    dt = exact(DT)
    for oxygen in oxygens:
        o = exact(oxygen)
        op = model.operator(p, o)
        mu = op["mu"]
        division_death = []
        diagonal = []
        for i, chromosome in enumerate(model.GRID):
            v = op["values"][i]
            survive = v["survival"]
            # Full binomial expectation of biological deaths. Pruning and
            # grid-drop losses remain in the live operator only.
            biological = (2 * v["lam"] * (1 - p["p_wgd"])
                          * (1 - (1 - v["probability"] * (1 - survive)) ** chromosome))
            need(biological > 0 and (biological - op["ledger"][i]["biological"]).contains(0),
                 f"biological death rate cross-check failed at {oxygen}/{chromosome}")
            division_death.append(biological)
            diagonal.append(1 - dt * (v["lam"] + mu[i]))
        need(all(v > 0 for v in diagonal), f"Euler positivity at {oxygen}")
        need(all(1 - dt * v["lam"] > 0 for v in op["values"]), f"capacity branch at {oxygen}")
        ledger.append(dict(oxygen_pct=oxygen, mask_sha256=op["mask_sha256"],
                           minimum_prune_gap=str(min(r["minimum_margin"] for r in op["masks"])),
                           minimum_live_diagonal=str(min(v.lower() for v in diagonal))))
        result[oxygen] = (op["G"], arb_mat([mu]), arb_mat([division_death]))
    table(out / "operators.tsv", ledger)
    return result


def profile(cohort):
    center = 46 if cohort == "2N" else 88
    weights = [Fraction(9**4, 9**4 + (n - center)**4) for n in model.GRID]
    total = sum(weights)
    return arb_mat([[exact(w / total)] for w in weights]), model.GRID.index(center)


def sums(matrix):
    return [sum(matrix[i, j] for i in range(N)) for j in range(matrix.ncols())]


def compose(count, stock, jac, stock_jac):
    return count + stock, [a + b for a, b in zip(sums(jac), stock_jac)]


def propagate(pi, pi_j, a, observed_start, wanted, op, eta, decay):
    """Track root derivative columns plus the last, local live inoculum column.

    The visible initial count is held fixed during this propagation. Initial
    dead-class mixture is immaterial at fixed oxygen and common clearance.
    """
    dt = exact(DT)
    G, mu, mitotic = op
    ncols = pi_j.ncols() + 1
    last = ncols - 1
    x = pi * a
    jac = arb_mat([[pi_j[i, j] * a for j in range(last)] + [pi[i, 0]] for i in range(N)])
    d = observed_start - a
    dj = [arb(0)] * last + [arb(-1)]
    one_minus_death = [1 - dt * mu[0, i] for i in range(N)]
    G = G * dt
    mu = mu * (dt * eta[0])
    mitotic = mitotic * (dt * eta[1])
    minimum_crowding, minimum_visible = arb(1), None
    sampled = {}
    wanted = set(wanted)
    final = max(wanted)
    for step in range(final + 1):
        L = sum(x[i, 0] for i in range(N))
        crowd = 1 - L
        need(crowd > 0, f"crowding branch at tick {step}")
        # The positivity branch follows from the strictly positive root
        # distribution, reset scale, nonnegative off-diagonals and Euler diagonal.
        minimum_crowding = min(minimum_crowding, crowd.lower())
        if step in wanted:
            c, cj = compose(L, d, jac, dj)
            need(c * CAPACITY > 1, f"log floor at tick {step}")
            minimum_visible = c.lower() if minimum_visible is None else min(minimum_visible, c.lower())
            sampled[step] = (x, jac, c, cj, L)
        if step == final:
            break
        gx = G * x
        gj = G * jac
        dl = sums(jac)
        dmx = (mitotic * x)[0, 0]
        dmj = mitotic * jac
        dhj = mu * jac
        d = decay * d + (mu * x)[0, 0] + crowd * dmx
        dj = [decay * dj[j] + dhj[0, j] + crowd * dmj[0, j] - dl[j] * dmx for j in range(ncols)]
        jac = arb_mat([
            [one_minus_death[i] * jac[i, j] + crowd * gj[i, j] - gx[i, 0] * dl[j]
             for j in range(ncols)] for i in range(N)
        ])
        x = arb_mat([[one_minus_death[i] * x[i, 0] + crowd * gx[i, 0]] for i in range(N)])
    return sampled, minimum_crowding, minimum_visible


def normalized_state(x, jac, L, fix_last):
    """Differentiate donor composition subject to chosen local count pivots."""
    k = jac.ncols() - 1
    v = x / L
    J = arb_mat([[jac[i, j] + jac[i, k] * fix_last[j] for j in range(k)] for i in range(N)])
    dl = sums(J)
    return v, arb_mat([[(J[i, j] - v[i, 0] * dl[j]) / L for j in range(k)] for i in range(N)])


def etas(case, cohort):
    fractions = ("2/5", "3/5") if case == "split" and cohort == "4N" else ("1/3", "2/3")
    return [exact(k) for k in fractions]


def root_columns(passages, roots):
    sizes = Counter()
    for p in passages.values():
        if p.samples:
            sizes[roots[p.id]] += len(p.samples) - 2
    need(sum(sizes.values()) == 107, "wrong number of count rows after seeding/last-count elimination")
    cols = {}
    for root_id, k in sizes.items():
        ref = 46 if passages[root_id].cohort == "2N" else 88
        choices = [i for i, n in enumerate(model.GRID) if n != ref]
        # Spread independent coordinate directions over the entire allowed grid.
        cols[root_id] = [choices[(j * (len(choices) - 1)) // max(1, k - 1)] for j in range(k)]
    return cols


def midpoint_preconditioned_bound(matrix):
    n = matrix.nrows()
    mid = arb_mat([[matrix[i, j].mid() for j in range(n)] for i in range(n)])
    identity = arb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    error = identity - mid.inv() * matrix
    return max(sum(abs(error[i, j]) for j in range(n)).upper() for i in range(n))


def summarize_matrix(matrix, root_id, cols, labels, out):
    n = matrix.nrows()
    need(n == matrix.ncols() == len(cols) == len(labels), "incorrect residual minor")
    entries = [
        dict(root_id=root_id, record_id=labels[i], perturbed_chromosome=model.GRID[cols[j]],
             **encode(matrix[i, j]))
        for i in range(n) for j in range(n)
    ]
    residual = midpoint_preconditioned_bound(matrix)
    need(residual < 1, f"cannot certify residual minor {root_id}: {residual}")
    det = matrix.det()
    return entries, dict(root_id=root_id, dimension=n, neumann_error_upper=str(residual),
                         determinant_enclosure=str(det), determinant_directly_excludes_zero=nonzero(det))


def run(phase, case, precision, out, witness_source=None):
    ctx.prec = precision
    out.mkdir(exist_ok=True, parents=True)
    passages, roots, uses = load_design(out)
    p = parameter_point(out)
    decay = 1 - exact(DT) * p["k_clear"]
    need(decay > 0 and decay < 1, "invalid death clearance")
    ops = operators(p, sorted({p.oxygen for p in passages.values()}, key=Fraction), out)

    allocations = {} if witness_source is None else {
        r["passage_id"]: int(r["live_initial_witness_cells"]) for r in read(witness_source)
    }
    cols = root_columns(passages, roots) if phase == "counts" else {}
    output_states = {}
    pivots, states, witness, row_entries, checks = [], [], [], [], []
    residual_rows = defaultdict(list)
    residual_labels = defaultdict(list)
    for psg in passages.values():
        key = roots[psg.id]
        selected = cols.get(key, ())
        k = len(selected)
        if psg.parent is None:
            pi, ref = profile(psg.cohort)
            pi_j = arb_mat([[int(i == j) - int(i == ref) for j in selected] for i in range(N)])
        else:
            parent = output_states[psg.parent]
            px, pj, _, _, pL = parent["samples"][psg.donor_tick]
            pi, pi_j = normalized_state(px, pj, pL, parent["fix_last"])

        count = exact(Fraction(psg.seed_count, CAPACITY))
        a_cells = allocations.get(psg.id, psg.seed_count * 4 // 5)
        if psg.parent is not None:
            sibling_count = uses[psg.parent][psg.donor_tick]
            safe = pL * CAPACITY / (2 * sibling_count)
            if witness_source is None:
                a_cells = min(a_cells, int(float(safe.lower())))
            need(exact(a_cells) < safe, f"not enough parent live cells for {psg.id}")
        a = exact(Fraction(a_cells, CAPACITY))
        eta = etas(case, psg.cohort)
        visible_rest = count - a
        residual_true_dead = 2 * visible_rest / sum(eta)
        q = a / (a + residual_true_dead)
        need(a > 0 and visible_rest > 0 and q > 0 and q < 1, f"invalid inoculum witness {psg.id}")
        witness.append(dict(
            passage_id=psg.id, live_initial_witness_cells=a_cells,
            initial_visible_count_conditioning_cells=psg.seed_count,
            inoculum_viable_fraction_witness=str(q),
            initial_dead_hypoxia_fraction_witness="1/2",
            hypoxia_visible_fraction_witness=str(eta[0].mid()),
            mitotic_visible_fraction_witness=str(eta[1].mid())
        ))
        wanted = set(uses.get(psg.id, ())) | {r[0] for r in psg.samples}
        samples, crowd_margin, visible_margin = propagate(
            pi, pi_j, a, count, wanted, ops[psg.oxygen], eta, decay
        )
        fix_last = [arb(0)] * k
        if psg.samples:
            if phase == "ols":
                pivot = sum(exact(weight) * samples[tick][3][-1] / samples[tick][2]
                            for tick, _, weight in psg.samples)
                label = psg.id
            else:
                final_tick, label, _ = psg.samples[-1]
                derivatives = samples[final_tick][3]
                pivot = derivatives[-1]
                need(nonzero(pivot), f"zero last count pivot for {psg.id}")
                fix_last = [-derivatives[i] / pivot for i in range(k)]
                for tick, record_id, _ in psg.samples[1:-1]:
                    cj = samples[tick][3]
                    residual_rows[key].append([cj[i] + cj[-1] * fix_last[i] for i in range(k)])
                    residual_labels[key].append(record_id)
            need(nonzero(pivot), f"zero {phase} nuisance pivot for {psg.id}")
            pivots.append(dict(passage_id=psg.id, row_label=label, root_id=key,
                               witness_derivative=str(pivot), absolute_lower=str(pivot.abs_lower()),
                               **encode(pivot)))
        output_states[psg.id] = dict(samples=samples, fix_last=fix_last)
        for tick, record_id, _ in psg.samples:
            c = samples[tick][2]
            states.append(dict(record_id=record_id, passage_id=psg.id, euler_step=tick,
                               witness_visible_cells=str(c * CAPACITY)))
        checks.append(dict(passage_id=psg.id, min_crowding_margin=str(crowd_margin),
                           sampled_visible_cells_lower=str(visible_margin * CAPACITY)))
    need(len(pivots) == 114, "wrong seed/amplitude pivot count")
    table(out / "nuisance_pivots.tsv", pivots)
    table(out / "witness_inocula.tsv", witness)
    table(out / "witness_counts.tsv", states)
    table(out / "branch_checks.tsv", checks)
    if phase == "counts":
        checks = []
        for key, chosen in cols.items():
            entries, certificate = summarize_matrix(
                arb_mat(residual_rows[key]), key, chosen, residual_labels[key], out
            )
            row_entries.extend(entries)
            checks.append(certificate)
            print(key, certificate["dimension"], certificate["neumann_error_upper"], flush=True)
        need(sum(r["dimension"] for r in checks) == 107, "missing composition pivots")
        table(out / "root_residual_minors.tsv", row_entries)
        table(out / "root_rank_checks.tsv", checks)
    proof_rank = 114 if phase == "ols" else 2 * 114 + 107
    table(out / "rank.tsv", [dict(phase=phase, visibility=case, precision_bits=precision,
                                  observations=proof_rank, nuisance_rank_certified=proof_rank,
                                  joint_rank_certified=proof_rank,
                                  mechanisms_with_exact_local_compensation=14)])
    print(f"{phase} {case}: certified nuisance and joint row rank {proof_rank}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True, choices=("ols", "counts"))
    parser.add_argument("--visibility", required=True, choices=("shared", "split"))
    parser.add_argument("--bits", type=int, default=384)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--witness", type=Path)
    args = parser.parse_args()
    run(args.phase, args.visibility, args.bits, args.out, args.witness)
