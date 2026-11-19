"""Independently enclose the continuous culture observation Jacobian with mpmath.

The only old-study inputs are witness definitions, branch-hash references,
and the preregistered output rows and parameter columns. No previous matrix,
determinant, derivative implementation or numerical result is read.
"""
import csv
import hashlib
import math
import sys
from fractions import Fraction
from pathlib import Path

import mpmath
from mpmath import iv, mp

iv.dps = mp.dps = 90
GRID = range(22, 155)
N = len(GRID)
COHORTS = ("2N", "4N")
OXYGENS = ("1/5", "1", "41/2")
INITIAL_TOTAL = 10**6
CAPACITY = 10**12
EPSILON = iv.mpf("1e-8")
H = iv.mpf("1e-35")
RADIUS = iv.mpf("1e-10")
ZERO = iv.mpf(0)
LOG_RATIO = [iv.ln(iv.mpf(n) / 44) for n in GRID]
DISTANCE2 = [d * d for d in range(N)]
BINOMIAL = [tuple(math.comb(n, m) for m in range(n + 1)) for n in GRID]


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def rational(text):
    q = Fraction(text)
    return iv.mpf(q.numerator) / q.denominator


def positive_log(value):
    if not value.real > 0:
        raise ValueError("Logarithm argument intersects the principal branch cut")
    return iv.ln(value)


def positive_divide(numerator, denominator):
    if not denominator.real > 0:
        raise ValueError("Divisor not certified nonzero on the analytic neighborhood")
    return numerator / denominator


def power(base, exponent):
    return iv.exp(exponent * positive_log(base))


def witness_values(stage, witness):
    values = {}
    for row in read(stage / "inputs/witnesses.tsv"):
        if row["witness"] != witness:
            continue
        low, high, fraction = map(rational, (row["lower"], row["upper"], row["fraction"]))
        rule = row["rule"]
        if rule == "explicit":
            x = rational(row["value"])
        elif rule == "geometric":
            x = low * power(high / low, fraction)
        elif rule == "log1p":
            x = (1 + low) * power((1 + high) / (1 + low), fraction) - 1
        elif rule == "linear":
            x = low + fraction * (high - low)
        else:
            raise ValueError(rule)
        assert x > low and x < high
        values[row["parameter"]] = x
    assert len(values) == 20
    return values


def generator(p, oxygen, masks=None):
    """Build sparse gain/loss destinations from the physical event equations."""
    hypoxia = positive_divide(1, 1 + iv.exp(p["n_O"] * (iv.ln(rational(oxygen)) - positive_log(p["O2_crit"]))))
    gain_columns, death, new_masks = [], [], []
    for j, n in enumerate(GRID):
        log_ratio = LOG_RATIO[j]
        division = positive_divide(p["lam_max"], 1 + p["alpha_o2"] * hypoxia * iv.exp(p["gamma_growth"] * log_ratio))
        direct_death = p["mu_hp"] * hypoxia * iv.exp(p["gamma_mu"] * log_ratio)
        error = p["p_mis_base"] + positive_divide(p["p_misseg"] * direct_death, p["k_o_mis"] + direct_death)
        log_survival = positive_log(p["buffer_smax"]) - p["buffer_beta"] * iv.exp(-p["buffer_n_exp"] * log_ratio)
        ordinary = division * (1 - p["p_wgd"])
        moves = [(j, -division)]
        if 2 * n <= GRID[-1]:
            moves.append((2 * n - GRID[0], division * p["p_wgd"]))
        if masks is None:
            bits = []
            for m, binomial in enumerate(BINOMIAL[j]):
                probability = binomial * error**m * (1 - error)**(n - m)
                if probability > EPSILON:
                    bits.append(True)
                elif probability < EPSILON:
                    bits.append(False)
                else:
                    raise ValueError(f"Mask not strict at O={oxygen} N={n} m={m}")
            new_masks.append(tuple(m for m, keep in enumerate(bits) if keep))
            kept = new_masks[-1]
        else:
            kept = masks[j]
        for m in kept:
            viable = ordinary * BINOMIAL[j][m] * error**m * (1 - error)**(n - m)
            if m:
                viable *= iv.exp(m * log_survival)
            for dst in (j - m, j + m):
                if 0 <= dst < N:
                    moves.append((dst, viable))
        gain_columns.append(moves)
        death.append(direct_death)
    return (gain_columns, death), tuple(new_masks)


def branch_hash(masks):
    rows = []
    for n, kept in zip(GRID, masks):
        on = set(kept)
        rows.append(f"{n}:" + "".join("1" if m in on else "0" for m in range(n + 1)))
    return hashlib.sha256("\n".join(rows).encode()).hexdigest()


def gaussian_initial(p, cohort):
    mean, sd = p["init_mean_" + cohort], p["init_sd_" + cohort]
    width_factor = positive_divide(1, 2 * sd**2)
    weights = [iv.exp(-(n - mean)**2 * width_factor) for n in GRID]
    denominator = sum(weights)
    return [positive_divide(w, denominator) for w in weights]


def live_jets(operator, probability):
    """Directional chain rule for four time derivatives of a quadratic field."""
    columns, death = operator
    initial = [INITIAL_TOTAL * v for v in probability]

    def gain(state):
        answer = [ZERO] * N
        for j, moves in enumerate(columns):
            for i, rate in moves:
                answer[i] += rate * state[j]
        return answer

    crowding = 1 - iv.mpf(INITIAL_TOTAL) / CAPACITY
    initial_gain = gain(initial)

    def first_partial(state, state_gain, live):
        return [crowding * g - d * x - live * g0 / CAPACITY
                for g, d, x, g0 in zip(state_gain, death, state, initial_gain)]

    def second_partial(v_gain, v_live, w_gain, w_live):
        return [-(v_live * gw + w_live * gv) / CAPACITY for gv, gw in zip(v_gain, w_gain)]

    x1 = [crowding * g - d * x for g, d, x in zip(initial_gain, death, initial)]
    g1, l1 = gain(x1), sum(x1)
    x2 = first_partial(x1, g1, l1)
    g2, l2 = gain(x2), sum(x2)
    x3 = [v + w for v, w in zip(first_partial(x2, g2, l2), second_partial(g1, l1, g1, l1))]
    g3, l3 = gain(x3), sum(x3)
    x4 = [v + 3 * w for v, w in zip(first_partial(x3, g3, l3), second_partial(g1, l1, g2, l2))]
    return {1: l1, 2: l2 / 2, 3: l3 / 6, 4: sum(x4) / 24}


def karyotype_moment(width, initial_probability):
    width_factor = positive_divide(1, 2 * width**2)
    weights = [iv.exp(-d2 * width_factor) for d2 in DISTANCE2]
    center_distance2 = [(n - 44)**2 for n in GRID]
    expectation = ZERO
    for mother, probability in enumerate(initial_probability):
        column = [weights[abs(i - mother)] for i in range(N)]
        expectation += probability * positive_divide(sum(w * d2 for w, d2 in zip(column, center_distance2)), sum(column))
    return expectation


def observation(p, masks, base_operators=None, base_initials=None, base_moment=None, only=None):
    """Return the declared Taylor observations and, if requested, the 2N moment."""
    initials = base_initials or {cohort: gaussian_initial(p, cohort) for cohort in COHORTS}
    answer = {}
    if only != "sigma_kary":
        for oxygen in OXYGENS:
            operator = base_operators[oxygen] if base_operators else generator(p, oxygen, masks[oxygen])[0]
            for cohort in COHORTS:
                jets = live_jets(operator, initials[cohort])
                answer.update({(oxygen, cohort, order): value for order, value in jets.items()})
    if base_moment is not None:
        answer["moment"] = base_moment
    else:
        answer["moment"] = karyotype_moment(p["sigma_kary"], initials["2N"])
    return answer


def endpoint_number(interval):
    return mp.make_mpf(interval._mpi_[1])


def dyadic_bounds(interval):
    endpoints = {}
    for label, (sign, mantissa, exponent, _) in zip(("lower", "upper"), interval._mpi_):
        endpoints[f"{label}_mantissa"] = (-1)**sign * mantissa
        endpoints[f"{label}_exponent"] = exponent
    return endpoints


def full_matrix_rank(jacobian):
    """An arbitrary approximate inverse is checked with outward intervals."""
    n = len(jacobian)
    middle = mp.matrix([[mp.make_mpf(v.mid._mpi_[0]) for v in row] for row in jacobian])
    inverse = middle**-1
    product = iv.matrix(inverse.tolist()) * iv.matrix(jacobian)
    residual = iv.eye(n) - product
    norm_upper = max(endpoint_number(sum(abs(residual[i, j]) for j in range(n))) for i in range(n))
    assert mp.isfinite(norm_upper) and norm_upper < 1
    det = iv.det(iv.matrix(jacobian))
    assert 0 not in det and all(mp.isfinite(mp.make_mpf(endpoint)) for endpoint in det._mpi_)
    return norm_upper, det


def verify(stage, witness, output):
    case = next(row for row in read(stage / "inputs/cases.tsv") if row["case"] == "B04")
    assert [Fraction(x) for x in case["oxygens"].split(",")] == list(map(Fraction, OXYGENS))
    assert case["max_order"] == "4" and case["unknown_initial"] == "True"
    assert case["unknown_observer"] == "False"
    fixed = {row["setting"]: row["value"] for row in read(stage / "inputs/fixed_settings.tsv")
             if row["scope"] == "archived_cfg"}
    assert Fraction(fixed["init_total_size"]) == INITIAL_TOTAL
    assert Fraction(fixed["K"]) == CAPACITY
    assert int(fixed["N_MIN"]) == GRID[0] and int(fixed["N_MAX"]) == GRID[-1]
    p = witness_values(stage, witness)
    source_rows = stage / "verification" / ("B04_" + witness)
    candidates = read(source_rows / "candidate_registry.tsv")
    selected = next(row for row in candidates if row["candidate"] == "primary")
    column_names = selected["columns"].split(",")
    column_powers = [int(row["scale_power_two"]) for row in read(source_rows / "column_registry.tsv")]
    registry = {int(row["index"]): row for row in read(source_rows / "row_registry.tsv")}
    chosen_rows = [registry[int(i)] for i in selected["rows"].split(",")]
    names = [(row["oxygen"], row["cohort"], int(row["order"])) for row in chosen_rows]
    row_powers = [int(row["scale_power_two"]) for row in chosen_rows]
    assert len(names) == len(column_names) == 18 and len(set(names)) == 18
    assert all(name[2] in range(1, 5) for name in names)
    base_operators, masks = {}, {}
    for oxygen in OXYGENS:
        base_operators[oxygen], masks[oxygen] = generator(p, oxygen)
    expected = {row["oxygen"]: row["mask_sha256"] for row in read(source_rows / "branch_registry.tsv")}
    reconstructed = {oxygen: branch_hash(masks[oxygen]) for oxygen in OXYGENS}
    assert reconstructed == expected, (reconstructed, expected)
    print(witness, "3 strict branch hashes match", flush=True)

    base_initials = {cohort: gaussian_initial(p, cohort) for cohort in COHORTS}
    base_moment = karyotype_moment(p["sigma_kary"], base_initials["2N"])

    # In the fixed-mask analytic continuation, vary all coordinates over a
    # complex rectangle. It contains the radius-R disk in every coordinate.
    # All denominators have nonzero enclosures in this rectangle. Since the
    # center is real, the symmetric quotient omits even Taylor terms. Cauchy
    # bounds the discarded odd terms by M*h^2/(R^3*(1-(h/R)^2)).
    rectangle = iv.mpf([-RADIUS.b, RADIUS.b])
    complex_box = {name: iv.mpc(value + rectangle, rectangle) for name, value in p.items()}
    box_outputs = observation(complex_box, masks)
    print(witness, "analytic domain verified on complex boxes", flush=True)
    jacobian_columns = []
    for name in column_names + ["sigma_kary"]:
        replace_initial = name.startswith("init_")
        use_base = replace_initial or name == "sigma_kary"
        use_base_moment = not (name.endswith("_2N") or name == "sigma_kary")
        evaluations = []
        for step in (-H, H):
            variation = dict(p, **{name: p[name] + step})
            evaluations.append(observation(
                variation,
                masks,
                base_operators if use_base else None,
                None if replace_initial else base_initials,
                base_moment if use_base_moment else None,
                only=name,
            ))
        minus, plus = evaluations
        result = []
        for key in names + ["moment"]:
            if key not in plus or (key == "moment" and use_base_moment):
                result.append(ZERO)
                continue
            error = abs(box_outputs[key]) * H**2 / (RADIUS**3 * (1 - (H / RADIUS)**2))
            assert mp.isfinite(endpoint_number(error))
            result.append((plus[key] - minus[key]) / (2 * H) + iv.mpf([-error.b, error.b]))
        jacobian_columns.append(result)
        print(witness, name, "enclosed", flush=True)
    row_scales = [iv.mpf(2)**n for n in row_powers] + [iv.mpf(1)]
    column_scales = [iv.mpf(2)**n for n in column_powers] + [iv.mpf(1)]
    matrix = [[jacobian_columns[j][i] * row_scales[i] * column_scales[j] for j in range(19)]
              for i in range(19)]
    write(output / f"{witness}_matrix.tsv", [
        dict(row=i, column=j, observation=str((names + ["moment"])[i]),
             parameter=(column_names + ["sigma_kary"])[j], enclosure=str(matrix[i][j]),
             **dyadic_bounds(matrix[i][j]))
        for i in range(19) for j in range(19)
    ])
    report_rows = []
    for dimension, target in [(18, [row[:18] for row in matrix[:18]]), (19, matrix)]:
        norm, determinant = full_matrix_rank(target)
        report_rows.append(dict(witness=witness, dimension=dimension, mpmath=mpmath.__version__,
                                decimal_digits=iv.dps, mask_hashes_matched=3, status="FULL_RANK",
                                residual_norm_upper=mp.nstr(norm, 12),
                                determinant=str(determinant), sigma_partial=str(matrix[18][18]),
                                **dyadic_bounds(determinant)))
    write(output / f"{witness}_results.tsv", report_rows)
    print(witness, "FULL RANK 18 and 19", flush=True)


if __name__ == "__main__":
    verify(Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3]))
