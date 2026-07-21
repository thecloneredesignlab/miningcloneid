from __future__ import annotations

import csv
import json
import math
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from concurrent.futures import ProcessPoolExecutor, as_completed
from scipy.integrate import solve_ivp
from time import time

from mcmc_fit import (
    run_mcmc, run_mcmc_joint,
    compute_aic_bic_single, compute_aic_bic_joint,
)

# ── Simulation settings ────────────────────────────────────────────────────────
MIN_SIZE      = 1e5
MAX_SIZE      = 2e10
DEFAULT_LEN   = 7.0
ODE_STEP      = 0.05

BASE_BEAM_WIDTH    = 10
BASE_MAX_DEPTH     = 10
SAMPLED_BEAM_WIDTH = 10
SAMPLED_MAX_DEPTH  = 10
N_SAMPLED_RUNS     = 10

# ── Joint fitting mode ─────────────────────────────────────────────────────────
# When JOINT_FIT_MODE = True the __main__ block loads ALL matched harvests from
# HARVESTS_CSV_PATH and fits them simultaneously via run_mcmc_joint().  Global
# biological parameters are shared; each mouse gets per-parameter log-space
# epsilon offsets so individual variation is captured while all mice jointly
# inform the global parameters.
#
# Set JOINT_FIT_MODE = False to restore the original single-harvest behaviour
# (the script is then patched per-harvest by run_one_harvest.py).
JOINT_FIT_MODE    = True
HARVESTS_CSV_PATH = "../data/InVivoData_Gemcitabine/harvest_ploidy_mapping.csv"

# ── In-vivo live/dead model constants ─────────────────────────────────────────
N_TR                         = 5      # Fixed transit compartments (from in-vitro model selection)
CONFLUENCE_DEATH_EXPONENT    = 4.0    # Fixed exponent q for confluence-death term
HAPLOID_N:    int            = 23
OXYGEN_HAPLOID_N: int        = 22

# k_tr is not fitted; it is interpolated from the tumor's mean ploidy using the
# log-linear delay model from the paper (Eq. log τ(p) = b0 + b1·log(p/2)),
# calibrated on the in-vitro 2N and 4N MAP estimates.
K_TR_2N = 6.84   # day⁻¹ — transit rate for near-diploid (2N) cells
K_TR_4N = 2.82   # day⁻¹ — transit rate for near-tetraploid (4N) cells

_TAU_2N   = N_TR / K_TR_2N                             # ≈ 0.731 days (mean delay, 2N)
_TAU_4N   = N_TR / K_TR_4N                             # ≈ 1.773 days (mean delay, 4N)
_K_TR_B0  = math.log(_TAU_2N)                          # intercept at ploidy ratio = 2
_K_TR_B1  = math.log(_TAU_4N / _TAU_2N) / math.log(2) # slope per doubling of ploidy


def _interpolate_k_tr(initial_ploidy: dict, haploid_n: int = HAPLOID_N) -> float:
    """
    Compute k_tr by interpolating from the weighted-mean ploidy of *initial_ploidy*.

    Log-linear model:  log τ(p) = b0 + b1·log(p/2),   k_tr = N_TR / τ
    where p is the mean ploidy ratio (chromosome count / haploid_n).

    Calibration:
        p = 2 (2N):  k_tr = 6.84 day⁻¹
        p = 4 (4N):  k_tr = 2.82 day⁻¹
    """
    total = sum(initial_ploidy.values())
    if total <= 0:
        return K_TR_2N
    weighted_chr = sum(N * cnt for N, cnt in initial_ploidy.items())
    mean_ploidy  = (weighted_chr / total) / haploid_n
    mean_ploidy  = max(mean_ploidy, 0.5)   # guard against degenerate input
    log_tau      = _K_TR_B0 + _K_TR_B1 * math.log(mean_ploidy / 2.0)
    return N_TR / math.exp(log_tau)

GEMCITABINE_INTRACELLULAR_CPEAK_COEFF = 0.00971   # C_peak (au) per mg/kg dose

_CONFIG_PATH = "../config/drug_kinetics.json"


def _load_drug_kinetics_config(path: str = _CONFIG_PATH) -> dict:
    with open(path) as fh:
        raw = json.load(fh)
    pk_params: dict[str, dict[str, dict]] = {}
    for drug, params in raw["PK_PARAMS_TO_FIT"].items():
        pk_params[drug] = {}
        for param, cfg in params.items():
            pk_params[drug][param] = {
                "init":           cfg["init"],
                "prior_log_mean": math.log(cfg["init"]),
                "prior_log_std":  cfg["prior_log_std"],
                "step":           cfg["step"],
            }
    return {
        "DOSE_REFERENCE_MG_KG": float(raw["DOSE_REFERENCE_MG_KG"]),
        "CYCLE_LENGTHS":        {k: float(v) for k, v in raw["CYCLE_LENGTHS"].items()},
        "PK_PARAMS_TO_FIT":     pk_params,
    }


_DK = _load_drug_kinetics_config()

EXCEL_PATH = "../data/InVivoData_Gemcitabine/dt_Gem_VT_20260209_v5.xlsx"

DOSE_REFERENCE_MG_KG: float = _DK["DOSE_REFERENCE_MG_KG"]

# Reference C_peak used for beta-dose correction  (= 0.00971 × 120 mg/kg)
_REF_C_PEAK = GEMCITABINE_INTRACELLULAR_CPEAK_COEFF * DOSE_REFERENCE_MG_KG


def load_harvest_data(
        excel_path: str,
        harvest_name: str,
        verbose: bool = True,
) -> tuple[dict[int, float], list[tuple[float, float, str, float]], str]:
    df = pd.read_excel(excel_path)
    matches = df[df["harvest"] == harvest_name]
    if matches.empty:
        raise ValueError(f"No rows found with harvest == '{harvest_name}'")

    row     = matches.iloc[0]
    cols    = list(df.columns)
    start_idx = cols.index("Day_0")
    day_cols  = cols[start_idx:]
    values    = row[day_cols].dropna()

    dose_mg_kg = float(row["Dose"])
    first_tx   = float(row["Day of 1st treatment"])

    burdens_cm3: dict[int, float] = {}
    for col, val in values.items():
        day = int(col.split("_")[1])
        if day < first_tx:
            continue
        fval = float(val)
        if fval > 0.0:
            burdens_cm3[day] = fval

    last_day = float(max(burdens_cm3.keys())) if burdens_cm3 else first_tx

    schedule: list[tuple[float, float, str, float]] = [
        (first_tx, last_day, "gemcitabine", dose_mg_kg),
    ]

    ploidy_name = str(row["harvest"]) + ".sps.cbs"

    if verbose:
        print(f"  Loaded harvest  : {harvest_name}")
        print(f"  Dose (mg/kg)    : {dose_mg_kg}")
        print(f"  First tx day    : {first_tx}")
        print(f"  Days with data  : {sorted(burdens_cm3.keys())}")
        print(f"  Drug schedule   : {schedule}")

    return burdens_cm3, schedule, ploidy_name


SAMPLE_NAME = "SUM159-2N-120-RL_harvest"

# In joint mode the module-level default-harvest data is only needed to
# initialise module globals; the joint loader reloads everything properly.
# Suppress the per-harvest print statements so the output isn't misleading.
_VERBOSE_LOAD = __name__ == "__main__" and not JOINT_FIT_MODE

_OBSERVED_TUMOR_BURDENS_CM3, OBSERVED_DRUGS_ADMINISTERED, PLOIDY_SAMPLE_NAME = (
    load_harvest_data(EXCEL_PATH, SAMPLE_NAME, verbose=_VERBOSE_LOAD)
)

FIRST_TX_DAY: float = next(
    (start for start, _end, drug, _dose in OBSERVED_DRUGS_ADMINISTERED
     if drug != "none"),
    0.0,
)

# Only gemcitabine and none are modelled
DRUGS = [
    "gemcitabine",
    "none",
]

CYCLE_LENGTHS: dict[str, float] = _DK["CYCLE_LENGTHS"]

PLOIDY_CSV_PATH:       str = "../data/InVivoData_Gemcitabine/all_ploidy.csv"
OXYGEN_COUNT_TSV_PATH: str = (
    "../data/InVivoData_Gemcitabine/"
    "treatment_day_ploidy_counts_fit_invivo_o2_supply_demand_MAP_seed28.tsv"
)

_CELLS_PER_CM3 = 1e7


def load_treatment_day_ploidy_from_oxygen(
        count_tsv_path:  str,
        harvest_name:    str,
        oxygen_haploid_n: int  = OXYGEN_HAPLOID_N,
        label_filter:    str | None = "live",
        verbose:         bool = True,
) -> dict[int, float]:
    df   = pd.read_csv(count_tsv_path, sep="\t")
    mask = df["harvest"] == harvest_name
    if label_filter is not None:
        mask &= df["label"] == label_filter
    sub = df[mask]
    if sub.empty:
        raise ValueError(
            f"No rows found for harvest='{harvest_name}', label='{label_filter}' "
            f"in '{count_tsv_path}'."
        )
    result: dict[int, float] = {}
    for _, row in sub.iterrows():
        chr_count = int(round(float(row["ploidy"]) * oxygen_haploid_n))
        result[chr_count] = result.get(chr_count, 0.0) + float(row["cell_count"])

    if verbose:
        total  = sum(result.values())
        tx_day = sub["treatment_day"].iloc[0]
        print(f"  Oxygen model ploidy : {harvest_name}")
        print(f"  Treatment day       : {tx_day}")
        print(f"  Label filter        : {label_filter}")
        print(f"  Ploidy states       : {len(result)}")
        print(f"  Total cells (init)  : {total:.3e}")

    return result


def load_ploidy_distribution(
        csv_path:    str  = PLOIDY_CSV_PATH,
        sample_name: str  = PLOIDY_SAMPLE_NAME,
        verbose:     bool = True,
) -> np.ndarray:
    try:
        df = pd.read_csv(csv_path, sep="\t")
    except FileNotFoundError:
        if verbose:
            print(f"  WARNING: ploidy CSV not found at '{csv_path}'. "
                  "OBSERVED_END_PLOIDY_DISTRIBUTION will be empty.")
        return np.array([], dtype=float)

    mask   = df["file"].str.contains(sample_name, na=False)
    values = df.loc[mask, "total_chromosomes"].to_numpy(dtype=float)

    if values.size == 0 and verbose:
        print(f"  WARNING: no rows matched sample '{sample_name}'. "
              f"Biopsy likelihood term disabled.")
    elif verbose:
        print(f"  Loaded {values.size} ploidy values for '{sample_name}' "
              f"(mean={values.mean():.4f}, std={values.std():.4f})")

    return values


INITIAL_PLOIDY: dict[int, float] = load_treatment_day_ploidy_from_oxygen(
    OXYGEN_COUNT_TSV_PATH,
    SAMPLE_NAME,
    verbose=_VERBOSE_LOAD,
)

if INITIAL_PLOIDY and _OBSERVED_TUMOR_BURDENS_CM3:
    _anchor_day = int(FIRST_TX_DAY)
    if _anchor_day not in _OBSERVED_TUMOR_BURDENS_CM3:
        raise ValueError(
            f"No observed tumor burden at FIRST_TX_DAY={_anchor_day}. "
            f"Observed days: {sorted(_OBSERVED_TUMOR_BURDENS_CM3.keys())}."
        )
    _initial_cells   = sum(INITIAL_PLOIDY.values())
    _initial_vol_cm3 = _OBSERVED_TUMOR_BURDENS_CM3[_anchor_day]
    if _initial_cells > 0 and _initial_vol_cm3 > 0:
        _CELLS_PER_CM3 = _initial_cells / _initial_vol_cm3
        if _VERBOSE_LOAD:
            print(f"  Derived _CELLS_PER_CM3 = {_CELLS_PER_CM3:.3e} cells/cm^3")

OBSERVED_TUMOR_BURDENS: dict[int, float] = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}

OBSERVED_END_PLOIDY_DISTRIBUTION: np.ndarray = load_ploidy_distribution(
    verbose=_VERBOSE_LOAD
)

PK_PARAMS_TO_FIT: dict[str, dict[str, dict]] = _DK["PK_PARAMS_TO_FIT"]

START_BEAM_FROM_OBSERVED_END = False
PLOT_OBSERVED_DATA:     bool = False


def get_cycle_length(drug: str) -> float:
    return CYCLE_LENGTHS.get(drug, DEFAULT_LEN)


# ── Plasma PK helper ───────────────────────────────────────────────────────────

def _pulsed_dose_value(t: float, C_peak: float, half_life: float, period: float) -> float:
    """IV bolus: exponential decay between doses of fixed period."""
    if C_peak <= 0.0:
        return 0.0
    lam  = math.log(2.0) / max(half_life, 1e-12)
    tmod = math.fmod(t, period) if period > 0 else t
    return C_peak * math.exp(-lam * tmod)


# ── In-vivo live/dead ODE (from invitro_fitting.py model) ─────────────────────
# Fitted parameters: r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf
# Fixed/derived:     N_TR = 5, k_tr = interpolated from INITIAL_PLOIDY
# State layout:      [A, Z_1, …, Z_N_TR, D]

def _invivo_livedead_rhs(
    t,
    y,
    r, K, k_tr, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    C_peak, half_life, period,
):
    """
    In-vivo live/dead ODE driven by plasma PK.

    dA/dt = r * (1/(1+k_cyto*C_eff)) * A*(1 - A/K)
            - [mu_base + mu_conf*(A/K)^q + k_kill*Z_n] * A
    dZ_i/dt = k_tr*(C_eff - Z_1)           [i=1]
              k_tr*(Z_{i-1} - Z_i)          [i=2..N_TR]
    dD/dt   = [mu_base + mu_conf*(A/K)^q + k_kill*Z_n]*A - k_clear*D
    """
    n   = N_TR
    A   = max(float(y[0]), 0.0)
    Z   = np.asarray(y[1:1 + n], dtype=float)
    D   = max(float(y[-1]), 0.0)

    # Raw plasma drug signal
    C_raw = _pulsed_dose_value(t, C_peak, half_life, period)

    # Beta dose-scaling correction:
    #   C_eff = C_raw * (C_peak / C_peak_ref)^(beta_dose - 1)
    if C_raw <= 0.0 or C_peak <= 0.0 or _REF_C_PEAK <= 0.0:
        C_eff = 0.0
    else:
        C_eff = C_raw * (C_peak / _REF_C_PEAK) ** (beta_dose - 1.0)
        C_eff = max(C_eff, 0.0)

    # Transit chain
    dZ    = np.zeros(n)
    dZ[0] = k_tr * (C_eff - Z[0])
    for i in range(1, n):
        dZ[i] = k_tr * (Z[i - 1] - Z[i])
    Z_n = max(float(Z[-1]), 0.0)

    # Death hazards
    density   = max(A / max(K, 1e-30), 0.0)
    mu_ctrl   = mu_base + mu_conf * (density ** CONFLUENCE_DEATH_EXPONENT)
    kappa_gem = k_kill * Z_n
    lam_total = mu_ctrl + kappa_gem

    # Cytostasis (immediate, driven by C_eff)
    growth_mult = 1.0 / (1.0 + k_cyto * C_eff) if k_cyto > 0.0 else 1.0

    dA    = r * growth_mult * A * (1.0 - A / max(K, 1e-30)) - lam_total * A
    dD    = lam_total * A - k_clear * D

    return [dA] + dZ.tolist() + [dD]


def _run_live_dead_ode(
    alive0:    float,
    dead0:     float,
    drug:      str,
    r:         float,
    K:         float,
    k_kill:    float,
    k_clear:   float,
    k_cyto:    float,
    beta_dose: float,
    mu_base:   float,
    mu_conf:   float,
    pk_overrides: dict,
    cycle_len: float,
    dt:        float = ODE_STEP,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """
    Simulate one drug course segment.

    k_tr is not a parameter here — it is derived from the module-level
    INITIAL_PLOIDY via _interpolate_k_tr() and fixed for the duration of
    the simulation.

    Returns (alive_final, dead_final, alive_trajectory, t_points).
    """
    drug = drug.lower()

    # k_tr is derived from the tumor's mean ploidy, not fitted
    k_tr = _interpolate_k_tr(INITIAL_PLOIDY, HAPLOID_N)

    if drug == "none":
        C_peak   = 0.0
        half_life = 1.0
        period    = cycle_len if cycle_len > 0 else DEFAULT_LEN
    else:
        # gemcitabine — use PK overrides (may include MCMC-fitted params)
        C_peak    = float(pk_overrides.get("C_peak",    GEMCITABINE_INTRACELLULAR_CPEAK_COEFF * DOSE_REFERENCE_MG_KG))
        half_life = float(pk_overrides.get("half_life", 0.05))
        period    = float(pk_overrides.get("period",    7.0))

    n    = N_TR
    y0   = np.array([alive0] + [0.0] * n + [dead0], dtype=float)

    n_points = max(int(cycle_len / dt) + 1, 10)
    t_eval   = np.linspace(0.0, cycle_len, n_points)

    args = (r, K, k_tr, k_kill, k_clear, k_cyto, beta_dose,
            mu_base, mu_conf, C_peak, half_life, period)

    # The plasma PK driver _pulsed_dose_value() uses fmod(t, period): it decays
    # within a period and then jumps discontinuously back to C_peak at every
    # t = m*period.  LSODA is a multistep method and cannot integrate across a
    # discontinuity — stepping over one invalidates its solution history and the
    # implicit corrector fails to converge, producing
    #   "lsoda: Repeated convergence failures (perhaps bad Jacobian ...)".
    # We therefore integrate the interval [0, cycle_len] piecewise, restarting
    # the solver at every dose boundary so each sub-integration sees a smooth
    # RHS.  Results are re-sampled onto the original uniform grid so the return
    # contract (uniformly-spaced trajectory) is unchanged.
    if period and period > 0:
        n_bounds   = int(math.floor(cycle_len / period - 1e-9))
        boundaries = [(m + 1) * period for m in range(n_bounds)
                      if (m + 1) * period < cycle_len - 1e-9]
    else:
        boundaries = []
    seg_edges = [0.0] + boundaries + [cycle_len]

    ys      = [y0.reshape(-1, 1)]   # column vector so axis-1 concat works
    ts      = [np.array([0.0])]
    y_start = y0
    try:
        with warnings.catch_warnings():
            # Extreme MCMC proposals can still momentarily stress the solver;
            # suppress the noisy per-step LSODA UserWarning (we act on
            # result.success instead) so the fit log stays readable.
            warnings.filterwarnings(
                "ignore", message=".*[Ll]soda.*", category=UserWarning)
            warnings.filterwarnings(
                "ignore", message=".*[Cc]onvergence failures.*")
            for t0, t1 in zip(seg_edges[:-1], seg_edges[1:]):
                if t1 <= t0:
                    continue
                sub_eval = t_eval[(t_eval > t0 + 1e-12) & (t_eval <= t1 + 1e-12)]
                # Start integration just after the jump so the segment RHS is
                # smooth; the value at exactly t0 is the previous segment's end.
                result = solve_ivp(
                    _invivo_livedead_rhs,
                    (t0, t1),
                    y_start,
                    t_eval=sub_eval if sub_eval.size else None,
                    method="LSODA",
                    args=args,
                    rtol=1e-7, atol=1e-9, max_step=dt,
                )
                if not result.success or result.y.shape[1] < 1:
                    return (alive0, dead0,
                            np.array([alive0, alive0]),
                            np.array([0.0, cycle_len]))
                if sub_eval.size:
                    ys.append(result.y)
                    ts.append(result.t)
                # Carry the true end-of-segment state into the next segment.
                y_start = result.y[:, -1] if sub_eval.size else y_start
    except Exception:
        return alive0, dead0, np.array([alive0, alive0]), np.array([0.0, cycle_len])

    if len(ys) < 2:
        return alive0, dead0, np.array([alive0, alive0]), np.array([0.0, cycle_len])

    Y = np.concatenate(ys, axis=1)
    T = np.concatenate([np.atleast_1d(t) for t in ts])

    alive_traj  = np.clip(Y[0],  0.0, np.inf)
    dead_traj   = np.clip(Y[-1], 0.0, np.inf)
    alive_final = float(alive_traj[-1])
    dead_final  = float(dead_traj[-1])

    return alive_final, dead_final, alive_traj, T


# ── Schedule helpers ───────────────────────────────────────────────────────────

def _fill_gaps_with_none(
        schedule: list[tuple[float, float, str, float]],
) -> list[tuple[float, float, str, float]]:
    """Insert ('none', 0.0) segments for any gaps in *schedule*."""
    if not schedule:
        return []
    filled: list[tuple[float, float, str, float]] = []
    for start, end, drug, dose in schedule:
        if filled and filled[-1][1] < start:
            filled.append((filled[-1][1], start, "none", 0.0))
        filled.append((start, end, drug, dose))
    return filled


def _pk_overrides_for(
        drug:       str,
        pk_state:   dict[str, dict] | None,
        dose_mg_kg: float | None = None,
) -> dict:
    """Build PK override kwargs for *drug*.

    For gemcitabine, C_peak is computed from dose when *dose_mg_kg* is known;
    otherwise it is left as the MCMC-fitted value from pk_state.
    """
    overrides: dict = dict(pk_state.get(drug, {})) if pk_state else {}
    if drug == "gemcitabine" and dose_mg_kg is not None:
        overrides["C_peak"] = GEMCITABINE_INTRACELLULAR_CPEAK_COEFF * dose_mg_kg
    return overrides


# ── Core simulation functions ──────────────────────────────────────────────────

def simulate_next_state(
    alive0:       float,
    drug:         str,
    r:            float,
    K:            float,
    k_kill:       float,
    k_clear:      float,
    k_cyto:       float,
    beta_dose:    float,
    mu_base:      float,
    mu_conf:      float,
    pk_overrides: dict | None = None,
    cycle_len:    float | None = None,
) -> tuple[float, np.ndarray]:
    """Simulate one cycle; returns (alive_final, alive_trajectory)."""
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    alive_final, _dead_final, alive_traj, _t = _run_live_dead_ode(
        alive0=alive0, dead0=0.0, drug=drug,
        r=r, K=K, k_kill=k_kill, k_clear=k_clear,
        k_cyto=k_cyto, beta_dose=beta_dose, mu_base=mu_base, mu_conf=mu_conf,
        pk_overrides=overrides, cycle_len=T,
    )
    return alive_final, alive_traj


def get_observed_end_ploidy(
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    pk_state=None,
):
    """
    Simulate the observed drug schedule and return the final alive-cell count.

    Returns a dict {0: alive_final} for API compatibility.  Since the in-vivo
    live/dead model does not track karyotype, this is used only for beam-search
    initialisation when START_BEAM_FROM_OBSERVED_END is True.
    """
    if not OBSERVED_DRUGS_ADMINISTERED:
        return {0: float(sum(INITIAL_PLOIDY.values()))}

    alive = float(sum(INITIAL_PLOIDY.values()))
    schedule = _fill_gaps_with_none(OBSERVED_DRUGS_ADMINISTERED)

    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            alive, _ = simulate_next_state(
                alive, drug, r, K, k_kill, k_clear, k_cyto,
                beta_dose, mu_base, mu_conf,
                pk_overrides=overrides, cycle_len=cycle_len,
            )
        except Exception:
            return None
        if alive < MIN_SIZE or not math.isfinite(alive):
            return None

    return {0: alive}


def _simulate_burden_timeline(
    observed_schedule,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    pk_state=None,
):
    """Simulate entire observed schedule; return list of (day, alive_cells)."""
    alive    = float(sum(INITIAL_PLOIDY.values()))
    timeline = [(FIRST_TX_DAY, alive)]
    schedule = _fill_gaps_with_none(observed_schedule)

    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            new_alive, seg_traj = simulate_next_state(
                alive, drug, r, K, k_kill, k_clear, k_cyto,
                beta_dose, mu_base, mu_conf,
                pk_overrides=overrides, cycle_len=cycle_len,
            )
        except Exception:
            return None

        n_tp = len(seg_traj)
        if n_tp > 0:
            for t_idx, a_val in enumerate(seg_traj):
                day   = start_day + (t_idx + 1) * cycle_len / n_tp
                total = float(np.asarray(a_val).sum())
                timeline.append((day, total))

        alive = new_alive

    return timeline


def _simulate_burden_and_ploidy(
    observed_schedule,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    pk_state=None,
):
    """
    Simulate the full schedule ONCE and return (burden_timeline, ploidy_dict).

    The live/dead model does not track karyotype, so ploidy_dict is always
    returned as {} (an empty dict — the biopsy likelihood term is disabled).
    """
    alive    = float(sum(INITIAL_PLOIDY.values()))
    timeline = [(FIRST_TX_DAY, alive)]
    schedule = _fill_gaps_with_none(observed_schedule)

    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            new_alive, seg_traj = simulate_next_state(
                alive, drug, r, K, k_kill, k_clear, k_cyto,
                beta_dose, mu_base, mu_conf,
                pk_overrides=overrides, cycle_len=cycle_len,
            )
        except Exception:
            return None, None

        n_tp = len(seg_traj)
        if n_tp > 0:
            for t_idx, a_val in enumerate(seg_traj):
                day   = start_day + (t_idx + 1) * cycle_len / n_tp
                total = float(np.asarray(a_val).sum())
                timeline.append((day, total))

        if new_alive < MIN_SIZE or not math.isfinite(new_alive):
            return None, None

        alive = new_alive

    # No karyotype tracked — return empty ploidy dict
    return timeline, {}


# ── Beam search ────────────────────────────────────────────────────────────────

def _simulate_next_state_wrapper(
    alive0, drug, path, traj,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    pk_overrides, cycle_len=None,
):
    new_alive, seg_traj = simulate_next_state(
        alive0, drug, r, K, k_kill, k_clear, k_cyto,
        beta_dose, mu_base, mu_conf,
        pk_overrides=pk_overrides, cycle_len=cycle_len,
    )
    return new_alive, seg_traj, path, traj, drug, cycle_len


def _beam_search_step(
    current_beams,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    beam_width,
    pk_state=None,
    start_day=0.0,
):
    # Keep already-extinct beams unchanged
    next_candidates = [
        (burden, alive, path, traj, True)
        for burden, alive, path, traj, extinct in current_beams
        if extinct
    ]

    for burden, alive, path, traj, extinct in current_beams:
        if extinct:
            continue
        elapsed          = start_day + sum(d_len for _, _, d_len in path)
        available_drugs  = DRUGS if elapsed >= FIRST_TX_DAY else ["none"]

        for drug in available_drugs:
            beam_dose = (
                DOSE_REFERENCE_MG_KG if drug == "gemcitabine" else 0.0
            )
            overrides  = _pk_overrides_for(drug, pk_state,
                                            dose_mg_kg=beam_dose if drug != "none" else None)
            cycle_len  = (
                (FIRST_TX_DAY - elapsed)
                if drug == "none" and elapsed < FIRST_TX_DAY
                else None
            )

            new_alive, seg_traj, old_path, old_traj, drug_out, actual_cycle_len = (
                _simulate_next_state_wrapper(
                    alive, drug, path, traj,
                    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
                    overrides, cycle_len,
                )
            )
            new_burden    = new_alive
            resolved_len  = actual_cycle_len if actual_cycle_len is not None else get_cycle_length(drug_out)
            segment_info  = (drug_out, len(seg_traj), resolved_len)

            if new_burden < MIN_SIZE:
                next_candidates.append(
                    (new_burden, new_alive,
                     old_path + [segment_info], old_traj + list(seg_traj), True))
            elif new_burden <= MAX_SIZE:
                next_candidates.append(
                    (new_burden, new_alive,
                     old_path + [segment_info], old_traj + list(seg_traj), False))

    def _sort_key(x):
        burden, _alive, path, _traj, extinct = x
        if extinct:
            days = sum(d_len for _, _, d_len in path)
            return (0, days, burden)
        return (1, 0.0, burden)

    next_candidates.sort(key=_sort_key)
    return next_candidates[:beam_width]


def run_single_beam_search(
    run_idx,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    beam_width, max_depth,
    start_alive=None,
    pk_state=None,
    start_day=0.0,
):
    if start_alive is None:
        start_alive = float(sum(INITIAL_PLOIDY.values()))
    initial_burden = start_alive
    beam = [(initial_burden, start_alive, [],
             [np.array([start_alive])], False)]

    for _ in range(max_depth):
        beam = _beam_search_step(
            beam,
            r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
            beam_width, pk_state, start_day,
        )
        if not beam or all(b[4] for b in beam):
            break

    return beam[0] if beam else None


def _beam_search_worker(
    i,
    r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
    beam_width, max_depth,
    use_observed_end,
    pk_state,
):
    full_pk_state = {drug: dict(params) for drug, params in pk_state.items()}
    for drug, params in PK_PARAMS_TO_FIT.items():
        if drug not in full_pk_state:
            full_pk_state[drug] = {p: cfg["init"] for p, cfg in params.items()}

    if use_observed_end and OBSERVED_DRUGS_ADMINISTERED:
        end_state = get_observed_end_ploidy(
            r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
            pk_state=full_pk_state,
        )
        sa        = float(sum(end_state.values())) if end_state else float(sum(INITIAL_PLOIDY.values()))
        start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
    else:
        sa        = float(sum(INITIAL_PLOIDY.values()))
        start_day = FIRST_TX_DAY

    return run_single_beam_search(
        i,
        r, K, k_kill, k_clear, k_cyto, beta_dose, mu_base, mu_conf,
        beam_width, max_depth,
        start_alive=sa, pk_state=full_pk_state, start_day=start_day,
    )


# ── Joint-fitting helpers ──────────────────────────────────────────────────────

def _get_matching_harvests(csv_path: str) -> list[str]:
    """Return harvest names where has_match == True."""
    import csv as _csv
    harvests = []
    with open(csv_path, newline="") as fh:
        reader = _csv.DictReader(fh)
        for row in reader:
            if row["has_match"].strip() == "True":
                harvests.append(row["harvest"].strip())
    return harvests


def load_all_harvest_data_for_joint(harvest_names: list[str]) -> tuple[list[dict], list[str]]:
    """
    Load burden, ploidy and schedule data for every harvest in *harvest_names*.

    Returns (mice_data, valid_names) where mice_data is a list of per-mouse
    state dicts suitable for passing to run_mcmc_joint(), and valid_names is the
    subset of harvest_names that loaded without error.
    """
    mice_data:   list[dict] = []
    valid_names: list[str]  = []

    for harvest_name in harvest_names:
        try:
            burdens_cm3, schedule, ploidy_name = load_harvest_data(
                EXCEL_PATH, harvest_name, verbose=False,
            )
            first_tx = next(
                (start for start, _end, drug, _dose in schedule if drug != "none"),
                0.0,
            )
            initial_ploidy = load_treatment_day_ploidy_from_oxygen(
                OXYGEN_COUNT_TSV_PATH, harvest_name, verbose=False,
            )
            if not initial_ploidy:
                print(f"  [joint load] WARNING: no initial ploidy for {harvest_name} — skipped")
                continue

            anchor_day = int(first_tx)
            if anchor_day not in burdens_cm3:
                print(f"  [joint load] WARNING: no burden at day {anchor_day} "
                      f"for {harvest_name} — skipped")
                continue

            initial_cells = sum(initial_ploidy.values())
            initial_vol   = burdens_cm3[anchor_day]
            cells_per_cm3 = (initial_cells / initial_vol
                             if initial_cells > 0 and initial_vol > 0 else 1e7)

            obs = {day: vol * cells_per_cm3 for day, vol in burdens_cm3.items()}

            end_ploidy = load_ploidy_distribution(
                csv_path=PLOIDY_CSV_PATH, sample_name=ploidy_name, verbose=False,
            )

            mice_data.append({
                "harvest_name":   harvest_name,
                "initial_ploidy": initial_ploidy,
                "obs":            obs,
                "end_ploidy":     end_ploidy,
                "drugs":          schedule,
                "first_tx_day":   first_tx,
                "cells_per_cm3":  cells_per_cm3,
            })
            valid_names.append(harvest_name)
            print(f"  [joint load] {harvest_name}: {len(obs)} burden points, "
                  f"first_tx={first_tx}, cells/cm³={cells_per_cm3:.2e}")

        except Exception as exc:
            print(f"  [joint load] ERROR loading {harvest_name}: {exc}")

    return mice_data, valid_names


def plot_observed_fit(
    results_dir:       str,
    harvest_name:      str,
    map_params:        dict,
    pk_state:          dict,
    observed_schedule: list,
    observed_burdens:  dict,
    filename:          str = "observed_fit_tumor_burden.png",
    ylim:              tuple[float, float] | None = None,
) -> str | None:
    """
    Plot the FITTED (MAP) model simulated on the *actual* drug schedule the
    mouse received, overlaid on the observed tumor-burden data points.

    This is deliberately different from baseline_tumor_burden.png:
      • baseline_tumor_burden.png runs the fitted model forward on the
        beam-search *optimal* drug course (a counterfactual "what treatment
        should we have given"), so its curve is NOT expected to pass through
        the observed points.
      • observed_fit_tumor_burden.png (this plot) runs the fitted model on the
        *same drugs the mouse actually received*, so the predicted curve and
        the observed points are directly comparable — this is the plot to look
        at to judge whether the model fits the data.

    Requires the module globals used by _simulate_burden_timeline
    (INITIAL_PLOIDY, FIRST_TX_DAY) to already be set to this mouse's values.
    Returns the output path, or None if the simulation failed.

    If *ylim* (a (ymin, ymax) tuple) is given, the y-axis is fixed to that
    range instead of autoscaling. This is used to render a set of per-mouse
    plots that all share identical y-axis bounds, so they can be tiled into a
    grid whose panels are directly comparable in height (see
    regenerate_observed_fits.py --uniform and compile_fit_grid.py --uniform).
    """
    tl = _simulate_burden_timeline(
        observed_schedule,
        map_params["r"], map_params["K"], map_params["k_kill"],
        map_params["k_clear"], map_params["k_cyto"], map_params["beta_dose"],
        map_params["mu_base"], map_params["mu_conf"],
        pk_state=pk_state,
    )

    drug_colors = {"gemcitabine": "orange", "none": "yellow"}

    fig, ax = plt.subplots(figsize=(10, 5))

    # Shade the actual (observed) drug schedule, filling untreated gaps with
    # "none" so the whole observed window is accounted for.
    shaded: set = set()
    for start_day, end_day, drug, _dose in _fill_gaps_with_none(observed_schedule):
        ax.axvspan(start_day, end_day, color=drug_colors.get(drug, "gray"),
                   linewidth=0, alpha=0.15,
                   label=drug if drug not in shaded else None)
        shaded.add(drug)

    if tl is not None:
        days = np.array([t for t, _ in tl])
        vals = np.array([b for _, b in tl])
        ax.plot(days, vals, color="steelblue", linewidth=1.8,
                label="Predicted (fitted model, observed drugs)")
    else:
        print(f"  [{harvest_name}] WARNING: observed-drug simulation failed — "
              f"plotting data points only")

    obs_days = sorted(observed_burdens.keys())
    obs_vals = [observed_burdens[d] for d in obs_days]
    ax.scatter(obs_days, obs_vals, color="firebrick", zorder=5, s=50,
               label="Observed")

    ax.set_yscale("log")
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_xlabel("Day")
    ax.set_ylabel("Alive cells (log scale)")
    ax.set_title(f"Tumor Burden — Fitted Model on Observed Drugs\n({harvest_name})")
    ax.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1))
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path = f"{results_dir}/{filename}"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  [{harvest_name}] Saved: {out_path}")
    return out_path


def _save_mouse_results(
    harvest_name:    str,
    results_dir:     str,
    map_params:      dict,
    posterior_samps: list[dict],   # list of per-mouse param dicts (natural units)
    all_iter_params: list[dict],   # this mouse's own MAP trace (global + its epsilon), one entry per iteration
    mouse_data:      dict,
    pk_state_global: dict,
) -> None:
    """
    Run beam search and save all per-mouse outputs (CSVs + plots + animation)
    for one harvest.

    Before calling this function the caller must have set the module-level
    globals INITIAL_PLOIDY, FIRST_TX_DAY, OBSERVED_DRUGS_ADMINISTERED, and
    OBSERVED_TUMOR_BURDENS to this mouse's values so that the simulation
    functions use the correct data.
    """
    import csv as _csv

    os.makedirs(results_dir, exist_ok=True)

    r_map         = map_params["r"]
    K_map         = map_params["K"]
    k_kill_map    = map_params["k_kill"]
    k_clear_map   = map_params["k_clear"]
    k_cyto_map    = map_params["k_cyto"]
    beta_dose_map = map_params["beta_dose"]
    mu_base_map   = map_params["mu_base"]
    mu_conf_map   = map_params["mu_conf"]
    pk_state_map  = pk_state_global

    baseline_start_day = FIRST_TX_DAY
    map_start_alive    = float(sum(INITIAL_PLOIDY.values()))

    if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
        end_state_map = get_observed_end_ploidy(
            r_map, K_map, k_kill_map, k_clear_map,
            k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
            pk_state=pk_state_map,
        )
        if end_state_map:
            map_start_alive    = float(sum(end_state_map.values()))
            baseline_start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]

    # ── Baseline beam search ───────────────────────────────────────────────────
    baseline_res = run_single_beam_search(
        "baseline",
        r_map, K_map, k_kill_map, k_clear_map,
        k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
        BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
        start_alive=map_start_alive,
        pk_state=pk_state_map,
        start_day=baseline_start_day,
    )
    baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
    print(f"  [{harvest_name}] Baseline path: {baseline_path}")

    # ── Sampled beam searches ──────────────────────────────────────────────────
    rng          = np.random.default_rng()
    sel_idx      = rng.choice(len(posterior_samps), size=min(N_SAMPLED_RUNS, len(posterior_samps)),
                              replace=True)
    sel_samps    = [posterior_samps[i] for i in sel_idx]

    sampled_results: list = []
    sampled_params:  list = []
    n_full_maxout = 0
    use_obs_end   = START_BEAM_FROM_OBSERVED_END and bool(OBSERVED_DRUGS_ADMINISTERED)

    _env_workers  = int(os.environ.get("HARVEST_MAX_WORKERS", 0))
    _beam_workers = min(len(sel_samps), _env_workers or os.cpu_count())
    with ProcessPoolExecutor(max_workers=_beam_workers) as _pool:
        future_map = {
            _pool.submit(
                _beam_search_worker, i,
                float(sel_samps[i]["r"]),
                float(sel_samps[i]["K"]),
                float(sel_samps[i]["k_kill"]),
                float(sel_samps[i]["k_clear"]),
                float(sel_samps[i]["k_cyto"]),
                float(sel_samps[i]["beta_dose"]),
                float(sel_samps[i]["mu_base"]),
                float(sel_samps[i]["mu_conf"]),
                SAMPLED_BEAM_WIDTH, SAMPLED_MAX_DEPTH,
                use_obs_end,
                sel_samps[i]["pk"],
            ): i
            for i in range(len(sel_samps))
        }
        uniform_w = 1.0 / max(len(sel_samps), 1)
        run_weights: list[float] = []
        for future in as_completed(future_map):
            i   = future_map[future]
            res = future.result()
            if res is not None:
                sampled_results.append(res)
                run_weights.append(uniform_w)
                sampled_params.append(sel_samps[i])
            else:
                n_full_maxout += 1

    run_weights_arr = np.array(run_weights)
    if run_weights_arr.sum() > 0:
        run_weights_arr /= run_weights_arr.sum()

    print(f"  [{harvest_name}] Maxouts: {n_full_maxout}/{len(sel_samps)}")

    # ── Flip-rate table ────────────────────────────────────────────────────────
    flip_rate_rows = []
    _cycle_day     = baseline_start_day
    for i in range(len(baseline_path)):
        target_drug = baseline_path[i]
        cycle_len   = baseline_res[2][i][2]
        day_start   = int(round(_cycle_day))
        day_end     = int(round(_cycle_day + cycle_len))
        _cycle_day += cycle_len
        unweighted_flip = 0
        active_count    = 0
        for res in sampled_results:
            sp = [step[0] for step in res[2]]
            if i < len(sp):
                active_count += 1
                if sp[i] != target_drug or (i > 0 and sp[i - 1] != baseline_path[i - 1]):
                    unweighted_flip += 1
        raw_rate = (unweighted_flip / active_count) if active_count > 0 else 0.0
        flip_rate_rows.append({
            "cycle": i + 1, "baseline_drug": target_drug,
            "day_start": day_start, "day_end": day_end,
            "active_runs": active_count, "disagreement_rate": raw_rate,
        })

    _flip_csv = f"{results_dir}/flip_rate_table.csv"
    with open(_flip_csv, "w", newline="") as _fh:
        _writer = _csv.DictWriter(
            _fh,
            fieldnames=["cycle", "baseline_drug", "day_start", "day_end",
                        "active_runs", "disagreement_rate"],
        )
        _writer.writeheader()
        _writer.writerows(flip_rate_rows)
    print(f"  [{harvest_name}] Saved: {_flip_csv}")

    # ── Sampled run paths CSV ──────────────────────────────────────────────────
    if sampled_results:
        max_cycles = max(len(res[2]) for res in sampled_results)
        path_rows  = []
        for run_idx, (res, samp) in enumerate(zip(sampled_results, sampled_params)):
            final_burden, _alive, path, _traj, extinct = res
            sp  = [step[0] for step in path]
            row: dict = {"run": run_idx + 1}
            for c in range(max_cycles):
                row[f"cycle_{c + 1}"] = sp[c] if c < len(sp) else ""
            row.update({
                "r": samp["r"], "K": samp["K"], "beta_dose": samp["beta_dose"],
                "end_burden": final_burden, "outcome": "EXTINCT" if extinct else "alive",
            })
            path_rows.append(row)

        _path_csv   = f"{results_dir}/sampled_run_paths.csv"
        _cycle_cols = [f"cycle_{c + 1}" for c in range(max_cycles)]
        _fieldnames = ["run"] + _cycle_cols + ["r", "K", "beta_dose", "end_burden", "outcome"]
        with open(_path_csv, "w", newline="") as _fh:
            _writer = _csv.DictWriter(_fh, fieldnames=_fieldnames)
            _writer.writeheader()
            _writer.writerows(path_rows)
        print(f"  [{harvest_name}] Saved: {_path_csv}")

    # ── Baseline burden plot ───────────────────────────────────────────────────
    baseline_burden_timeline: list[tuple[float, float]] = []
    alive_state = map_start_alive
    day         = baseline_start_day
    baseline_burden_timeline.append((day, alive_state))

    if baseline_res is not None:
        for _drug, _n_seg, cycle_len in baseline_res[2]:
            overrides = _pk_overrides_for(
                _drug, pk_state_map,
                dose_mg_kg=DOSE_REFERENCE_MG_KG if _drug != "none" else None,
            )
            new_alive, seg_traj = simulate_next_state(
                alive_state, _drug,
                r_map, K_map, k_kill_map, k_clear_map,
                k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
                pk_overrides=overrides, cycle_len=cycle_len,
            )
            n_tp = len(seg_traj)
            for t_idx, a_val in enumerate(seg_traj):
                t = day + (t_idx + 1) * cycle_len / n_tp
                baseline_burden_timeline.append((t, float(a_val)))
            day += cycle_len
            alive_state = new_alive

    burden_days   = np.array([t for t, _ in baseline_burden_timeline])
    burden_values = np.array([b for _, b in baseline_burden_timeline])

    drug_colors = {"gemcitabine": "orange", "none": "yellow"}

    fig1, ax1 = plt.subplots(figsize=(10, 5))
    current_time = float(baseline_start_day)
    shaded: set = set()
    if baseline_res:
        for dn, _, dur in baseline_res[2]:
            ax1.axvspan(current_time, current_time + dur,
                        color=drug_colors.get(dn, "gray"), linewidth=0, alpha=0.15,
                        label=dn if dn not in shaded else None)
            shaded.add(dn)
            current_time += dur

    ax1.plot(burden_days, burden_values, color="steelblue",
             linewidth=1.8, label="Predicted (baseline path)")

    obs_days_plot = sorted(OBSERVED_TUMOR_BURDENS.keys())
    obs_vals_plot = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days_plot]
    ax1.scatter(obs_days_plot, obs_vals_plot, color="firebrick", zorder=5,
                label="Observed", s=50)

    ax1.set_yscale("log")
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Alive cells (log scale)")
    ax1.set_title(f"Tumor Burden — Baseline Path\n({harvest_name})")
    ax1.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1))
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(f"{results_dir}/baseline_tumor_burden.png", dpi=150, bbox_inches="tight")
    plt.close(fig1)
    print(f"  [{harvest_name}] Saved: {results_dir}/baseline_tumor_burden.png")

    # ── Observed-drug fit plot (fitted model on the drugs actually given) ──────
    # This is the plot to inspect to judge model fit against the data points,
    # since (unlike the baseline plot above) it uses the mouse's real schedule.
    plot_observed_fit(
        results_dir, harvest_name,
        map_params={
            "r": r_map, "K": K_map, "k_kill": k_kill_map, "k_clear": k_clear_map,
            "k_cyto": k_cyto_map, "beta_dose": beta_dose_map,
            "mu_base": mu_base_map, "mu_conf": mu_conf_map,
        },
        pk_state=pk_state_map,
        observed_schedule=OBSERVED_DRUGS_ADMINISTERED,
        observed_burdens=OBSERVED_TUMOR_BURDENS,
    )

    # ── Fitting animation (using this mouse's own effective-param trace) ───────
    N_ANIM_FRAMES = min(200, len(all_iter_params))
    frame_indices = np.linspace(0, len(all_iter_params) - 1, N_ANIM_FRAMES, dtype=int)

    anim_timelines = []
    for fi, idx in enumerate(frame_indices):
        p  = all_iter_params[idx]
        tl = _simulate_burden_timeline(
            OBSERVED_DRUGS_ADMINISTERED,
            p["r"], p["K"], p["k_kill"], p["k_clear"],
            p["k_cyto"], p["beta_dose"], p["mu_base"], p["mu_conf"],
            pk_state=p["pk"],
        )
        if tl is not None:
            anim_timelines.append((np.array([t for t, _ in tl]),
                                   np.array([b for _, b in tl]), p))
        else:
            anim_timelines.append((np.array([]), np.array([]), p))

    fig_anim, ax_anim = plt.subplots(figsize=(12, 6))
    fig_anim.subplots_adjust(top=0.85, right=0.82)
    ax_anim.set_yscale("log")
    ax_anim.set_xlabel("Day")
    ax_anim.set_ylabel("Alive cells (log scale)")
    ax_anim.grid(True, alpha=0.3)

    # Tune the y-axis to the CONVERGED fit (final ~10% of frames) plus the
    # observed points, instead of padding the union of every burn-in frame by
    # 0.1x-10x.  The old padding buried the per-dose delay oscillations
    # (~0.2 decade) under ~2 decades of empty space, so the fit looked flat.
    _Y_PAD = 10 ** 0.35   # ~0.35 decade of headroom on each side
    _tail  = anim_timelines[max(0, len(anim_timelines)
                                - max(1, len(anim_timelines) // 10)):]
    _conv  = np.array([v for _, vals, _ in _tail for v in vals if v > 0],
                      dtype=float)
    if _conv.size and obs_vals_plot:
        _lo  = min(float(np.percentile(_conv, 1)), min(obs_vals_plot))
        _hi  = max(float(np.percentile(_conv, 99)), max(obs_vals_plot))
        y_lo, y_hi = _lo / _Y_PAD, _hi * _Y_PAD
    else:
        y_lo, y_hi = 1e3, 1e12
    all_dv = [d for days, _, _ in anim_timelines for d in days]
    x_hi   = (max(max(all_dv), max(obs_days_plot)) * 1.05
               if all_dv else max(obs_days_plot, default=100) * 1.05)
    ax_anim.set_xlim(-1, x_hi)
    ax_anim.set_ylim(y_lo, y_hi)
    ax_anim.scatter(obs_days_plot, obs_vals_plot, color="firebrick", zorder=5,
                    s=50, label="Observed")

    fit_line, = ax_anim.plot([], [], color="steelblue", linewidth=2.0,
                              label="Predicted")
    title_text = ax_anim.set_title("", fontsize=11, pad=12)
    ax_anim.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
    ghost_artists: list = []

    def _anim_init():
        fit_line.set_data([], [])
        title_text.set_text("")
        return [fit_line, title_text]

    def _anim_update(frame_num):
        days, vals, p = anim_timelines[frame_num]
        phase = "burn-in" if p["burnin"] else "sampling"
        if frame_num > 0:
            pd2, pv2, _ = anim_timelines[frame_num - 1]
            if len(pd2) > 0:
                gl, = ax_anim.plot(pd2, pv2, color="steelblue",
                                   alpha=0.08, linewidth=0.6, zorder=1)
                ghost_artists.append(gl)
        fit_line.set_data(days if len(days) > 0 else [], vals if len(vals) > 0 else [])
        title_text.set_text(
            f"Fitting [{harvest_name}] — Iter {p['iter'] + 1} [{phase}]  "
            f"(mouse MAP)  logP={p['logP']:.1f}\n"
            f"r={p['r']:.4f}  K={p['K']:.2e}  k_kill={p['k_kill']:.3e}  β={p['beta_dose']:.4g}"
        )
        return [fit_line, title_text] + ghost_artists

    import matplotlib.animation as _animation
    anim = _animation.FuncAnimation(
        fig_anim, _anim_update, init_func=_anim_init,
        frames=len(anim_timelines), interval=80, blit=False,
    )
    anim.save(f"{results_dir}/gibbs_fitting_animation.gif", writer="pillow", dpi=100)
    plt.close(fig_anim)
    print(f"  [{harvest_name}] Saved: {results_dir}/gibbs_fitting_animation.gif")


# ── Main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    start_time = time()

    RESULTS_DIR = "results"
    os.makedirs(RESULTS_DIR, exist_ok=True)

    if JOINT_FIT_MODE:
        # ══════════════════════════════════════════════════════════════════════
        # JOINT MODE: load all matched mice, fit together, beam-search per mouse
        # ══════════════════════════════════════════════════════════════════════
        # Always write to results_joint/ in joint mode — ignore any RESULTS_DIR
        # value injected by run_one_harvest.py / run_all_harvests.py patching.
        RESULTS_DIR = "results_joint"
        os.makedirs(RESULTS_DIR, exist_ok=True)
        import json as _json

        print("=" * 70)
        print("Joint hierarchical fitting — all mice simultaneously")
        print("=" * 70)

        # 1. Discover and load all matched harvests
        print("\nDiscovering harvests...")
        all_harvest_names = _get_matching_harvests(HARVESTS_CSV_PATH)
        print(f"  Found {len(all_harvest_names)} matched harvests in CSV")

        print("\nLoading harvest data...")
        mice_data, valid_names = load_all_harvest_data_for_joint(all_harvest_names)
        print(f"\n  Successfully loaded {len(mice_data)} mice")

        if not mice_data:
            raise RuntimeError("No valid mice could be loaded — cannot run joint MCMC.")

        # 2. Joint MCMC
        print(f"\nStarting joint MCMC with {len(mice_data)} mice...")
        joint_result = run_mcmc_joint(
            mice_data        = mice_data,
            pk_params_to_fit = PK_PARAMS_TO_FIT,
            sample_names     = valid_names,
            caller_module_name = "beam_search_flip_rate_wgd",
            verbose          = True,
        )

        global_map      = joint_result["global_map"]
        mice_eps_map    = joint_result["mice_eps_map"]
        mice_map_params = joint_result["mice_map_params"]
        post_samples    = joint_result["post_samples"]
        all_trace       = joint_result["all_trace"]
        weights         = joint_result["weights"]

        print(f"\nJoint MCMC done in {time() - start_time:.1f}s")

        # 3a. AIC / BIC of the joint model vs. the ploidy-burden datapoints
        _ic = compute_aic_bic_joint(
            caller           = sys.modules[__name__],
            global_map       = global_map,
            mice_map_params  = mice_map_params,
            mice_data        = mice_data,
            pk_params_to_fit = PK_PARAMS_TO_FIT,
            sample_names     = valid_names,
        )
        print(f"\n{'='*60}\nModel selection (vs. ploidy-burden datapoints)\n{'='*60}")
        print(f"  Burden datapoints (n) : {_ic['n_datapoints']}  "
              f"across {_ic['n_mice']} mice")
        print(f"  Free parameters   (k) : {_ic['k_params']}  "
              f"({_ic['n_global_bio']} global bio + 1 sigma_B + "
              f"{_ic['n_epsilon']} epsilons + {_ic['n_pk_params']} PK)")
        print(f"  log-likelihood        : {_ic['loglik']:.4f}")
        print(f"  AIC                   : {_ic['aic']:.4f}")
        print(f"  BIC                   : {_ic['bic']:.4f}")

        # 3. Save global parameter summary
        _global_out = {
            "global_map": {k: v for k, v in global_map.items() if k != "pk"},
            "pk_map":     global_map["pk"],
            "model_selection": _ic,
            "mice": [
                {
                    "harvest": name,
                    "epsilons_log": eps,
                    "effective_params": {k: v for k, v in mpp.items()
                                         if k not in ("pk", "sigma_B", "energy")},
                }
                for name, eps, mpp in zip(valid_names, mice_eps_map, mice_map_params)
            ],
        }
        _global_json = f"{RESULTS_DIR}/joint_global_params.json"
        with open(_global_json, "w") as _fh:
            _json.dump(_global_out, _fh, indent=2)
        print(f"\nSaved global params: {_global_json}")

        # 4. Per-mouse beam search and output
        print(f"\n{'='*70}")
        print("Per-mouse beam search")
        print(f"{'='*70}")

        for m, (harvest_name, mouse_data, mouse_map, mouse_eps) in enumerate(
            zip(valid_names, mice_data, mice_map_params, mice_eps_map)
        ):
            print(f"\n[{m+1}/{len(valid_names)}] {harvest_name}")

            mouse_results_dir = os.path.join(RESULTS_DIR, harvest_name)
            os.makedirs(mouse_results_dir, exist_ok=True)

            # Set module-level globals to this mouse's data so that all
            # simulation functions (and ProcessPoolExecutor workers that
            # fork from this process) use the correct per-mouse data.
            INITIAL_PLOIDY              = mouse_data["initial_ploidy"]
            FIRST_TX_DAY                = mouse_data["first_tx_day"]
            OBSERVED_DRUGS_ADMINISTERED = mouse_data["drugs"]
            OBSERVED_TUMOR_BURDENS      = mouse_data["obs"]

            # Extract this mouse's posterior samples (effective natural-scale params)
            mouse_post_samps = [s["mice_params"][m] for s in post_samples]

            # Build this mouse's OWN iteration-by-iteration trace for the
            # fitting animation.  all_trace[i]["mice"][m] holds mouse m's
            # effective (global + epsilon) params at iteration i — using
            # those (instead of the shared global-only trace) is essential:
            # a mouse whose epsilon differs from the global value would
            # otherwise be animated with a curve that isn't even its own
            # fit, and per-mouse dosing/kill dynamics would be averaged away.
            mouse_iter_trace = [
                {
                    "iter":   t["iter"],
                    "burnin": t["burnin"],
                    "pk":     t["pk"],
                    "logP":   t["logP"],
                    **t["mice"][m],
                }
                for t in all_trace
            ]

            _save_mouse_results(
                harvest_name    = harvest_name,
                results_dir     = mouse_results_dir,
                map_params      = mouse_map,
                posterior_samps = mouse_post_samps,
                all_iter_params = mouse_iter_trace,   # this mouse's own MAP trace
                mouse_data      = mouse_data,
                pk_state_global = global_map["pk"],
            )

        print(f"\n{'='*70}")
        print(f"Joint fitting complete.  Total time: {time() - start_time:.1f}s")
        print(f"Global params : {_global_json}")
        print(f"Per-mouse dirs: {RESULTS_DIR}/<harvest_name>/")
        print(f"{'='*70}")

    else:
        # ══════════════════════════════════════════════════════════════════════
        # SINGLE-HARVEST MODE (original behaviour, used by run_one_harvest.py)
        # ══════════════════════════════════════════════════════════════════════
        _map_result, posterior_samples, selected_weights, all_iter_params = run_mcmc(
            initial_ploidy=INITIAL_PLOIDY,
            observed_burdens=OBSERVED_TUMOR_BURDENS,
            observed_end_ploidy=OBSERVED_END_PLOIDY_DISTRIBUTION,
            observed_drugs=OBSERVED_DRUGS_ADMINISTERED,
            pk_params_to_fit=PK_PARAMS_TO_FIT,
            haploid_n=HAPLOID_N,
            sample_name=SAMPLE_NAME,
            fn_get_end_ploidy=get_observed_end_ploidy,
            fn_simulate_burden=_simulate_burden_timeline,
            fn_fill_gaps=_fill_gaps_with_none,
            verbose=True,
        )

        r_map         = _map_result["r"]
        K_map         = _map_result["K"]
        k_kill_map    = _map_result["k_kill"]
        k_clear_map   = _map_result["k_clear"]
        k_cyto_map    = _map_result["k_cyto"]
        beta_dose_map = _map_result["beta_dose"]
        mu_base_map   = _map_result["mu_base"]
        mu_conf_map   = _map_result["mu_conf"]
        pk_state_map  = _map_result["pk"]
        k_tr_map      = _interpolate_k_tr(INITIAL_PLOIDY, HAPLOID_N)  # derived, not fitted

        print(f"\nMCMC fitting done in {time() - start_time:.1f}s")
        print(f"  MAP: r={r_map:.5f}  K={K_map:.3e}  k_tr={k_tr_map:.4f} (ploidy-derived)  "
              f"k_kill={k_kill_map:.4e}  β={beta_dose_map:.4f}")
        print(f"  mu_base={mu_base_map:.4f}  mu_conf={mu_conf_map:.4f}  "
              f"k_cyto={k_cyto_map:.4e}  k_clear={k_clear_map:.4f}")
        print(f"  Posterior samples: {len(posterior_samples)}")

        # AIC / BIC of the model vs. this mouse's ploidy-burden datapoints
        _ic = compute_aic_bic_single(
            caller           = sys.modules[__name__],
            map_params       = _map_result,
            initial_ploidy   = INITIAL_PLOIDY,
            obs              = OBSERVED_TUMOR_BURDENS,
            obs_drugs        = OBSERVED_DRUGS_ADMINISTERED,
            pk_params_to_fit = PK_PARAMS_TO_FIT,
            first_tx_day     = FIRST_TX_DAY,
        )
        print(f"\n{'='*60}\nModel selection (vs. ploidy-burden datapoints)\n{'='*60}")
        print(f"  Burden datapoints (n) : {_ic['n_datapoints']}")
        print(f"  Free parameters   (k) : {_ic['k_params']}  "
              f"({_ic['n_bio_params']} bio + 1 sigma_B + {_ic['n_pk_params']} PK)")
        print(f"  log-likelihood        : {_ic['loglik']:.4f}")
        print(f"  AIC                   : {_ic['aic']:.4f}")
        print(f"  BIC                   : {_ic['bic']:.4f}")

        os.makedirs(RESULTS_DIR, exist_ok=True)
        _ic_json = f"{RESULTS_DIR}/model_selection_aic_bic.json"
        with open(_ic_json, "w") as _fh:
            json.dump(_ic, _fh, indent=2)
        print(f"  Saved: {_ic_json}")

        if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
            end_state_map = get_observed_end_ploidy(
                r_map, K_map, k_kill_map, k_clear_map,
                k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
                pk_state=pk_state_map,
            )
            map_start_alive  = float(sum(end_state_map.values())) if end_state_map else float(sum(INITIAL_PLOIDY.values()))
            baseline_start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
        else:
            map_start_alive    = float(sum(INITIAL_PLOIDY.values()))
            baseline_start_day = FIRST_TX_DAY

        baseline_res = run_single_beam_search(
            "baseline",
            r_map, K_map, k_kill_map, k_clear_map,
            k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
            BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
            start_alive=map_start_alive,
            pk_state=pk_state_map,
            start_day=baseline_start_day,
        )
        baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
        print(f"Baseline path: {baseline_path}")

        rng              = np.random.default_rng()
        selected_idx     = rng.choice(len(posterior_samples), size=N_SAMPLED_RUNS, replace=True)
        selected_samples = [posterior_samples[i] for i in selected_idx]

        sampled_results, run_weights, sampled_params = [], [], []
        use_observed_end = START_BEAM_FROM_OBSERVED_END and bool(OBSERVED_DRUGS_ADMINISTERED)

        n_full_maxout = 0

        _env_workers  = int(os.environ.get("HARVEST_MAX_WORKERS", 0))
        _beam_workers = min(N_SAMPLED_RUNS, _env_workers or os.cpu_count())
        with ProcessPoolExecutor(max_workers=_beam_workers) as pool:
            future_map = {
                pool.submit(
                    _beam_search_worker, i,
                    float(selected_samples[i]["r"]),
                    float(selected_samples[i]["K"]),
                    float(selected_samples[i]["k_kill"]),
                    float(selected_samples[i]["k_clear"]),
                    float(selected_samples[i]["k_cyto"]),
                    float(selected_samples[i]["beta_dose"]),
                    float(selected_samples[i]["mu_base"]),
                    float(selected_samples[i]["mu_conf"]),
                    SAMPLED_BEAM_WIDTH, SAMPLED_MAX_DEPTH,
                    use_observed_end,
                    selected_samples[i]["pk"],
                ): i
                for i in range(N_SAMPLED_RUNS)
            }
            for future in as_completed(future_map):
                i   = future_map[future]
                res = future.result()
                if res is not None:
                    sampled_results.append(res)
                    run_weights.append(selected_weights[i])
                    sampled_params.append(selected_samples[i])
                else:
                    n_full_maxout += 1

        run_weights   = np.array(run_weights)
        run_weights  /= run_weights.sum()

        print(f"\n  Full maxouts: {n_full_maxout}/{N_SAMPLED_RUNS} runs "
              f"({100 * n_full_maxout / N_SAMPLED_RUNS:.1f}%)")

        print("\n" + "=" * 78)
        print(f"{'Cycle':<7} | {'Baseline Drug':<16} | {'Days':<16} | {'Disagreement Rate':>11}")
        print("-" * 78)
        _cycle_day    = baseline_start_day
        flip_rate_rows = []
        for i in range(len(baseline_path)):
            target_drug = baseline_path[i]
            cycle_len   = baseline_res[2][i][2]
            day_start   = int(round(_cycle_day))
            day_end     = int(round(_cycle_day + cycle_len))
            days_str    = f"Day {day_start}–{day_end}"
            _cycle_day += cycle_len
            unweighted_flip = 0
            active_count    = 0
            for res, w in zip(sampled_results, run_weights):
                sampled_path = [step[0] for step in res[2]]
                if i < len(sampled_path):
                    active_count += 1
                    end_mismatch   = sampled_path[i] != target_drug
                    start_mismatch = (i > 0 and sampled_path[i - 1] != baseline_path[i - 1])
                    if end_mismatch or start_mismatch:
                        unweighted_flip += 1
            raw_rate = (unweighted_flip / active_count) if active_count > 0 else 0.0
            print(f"{i + 1:<7} | {target_drug:<16} | {days_str:<16} | {raw_rate * 100:>9.2f}%")
            flip_rate_rows.append({
                "cycle": i + 1, "baseline_drug": target_drug,
                "day_start": day_start, "day_end": day_end,
                "active_runs": active_count, "disagreement_rate": raw_rate,
            })

        print("=" * 78)

        _flip_csv = f"{RESULTS_DIR}/flip_rate_table.csv"
        os.makedirs(RESULTS_DIR, exist_ok=True)
        with open(_flip_csv, "w", newline="") as _fh:
            _writer = csv.DictWriter(
                _fh,
                fieldnames=["cycle", "baseline_drug", "day_start", "day_end",
                            "active_runs", "disagreement_rate"],
            )
            _writer.writeheader()
            _writer.writerows(flip_rate_rows)
        print(f"  Saved: {_flip_csv}")

        if sampled_results:
            max_cycles    = max(len(res[2]) for res in sampled_results)
            drug_col_w    = 14
            param_col_w   = 10
            burden_col_w  = 12

            cycle_headers = " | ".join(
                f"{'Cycle ' + str(c + 1):<{drug_col_w}}" for c in range(max_cycles)
            )
            param_header = (
                f"{'r':<{param_col_w}} | "
                f"{'K':<{param_col_w}} | "
                f"{'beta_dose':<{param_col_w}} | "
                f"{'End Burden':<{burden_col_w}} | "
                f"{'Outcome':<16}"
            )
            header = f"{'Run':<6} | {cycle_headers} | {param_header}"
            sep    = "-" * len(header)

            print("\n" + "=" * len(header))
            print("Sampled Run Paths")
            print("=" * len(header))
            print(header)
            print(sep)

            path_rows = []
            for run_idx, (res, samp) in enumerate(zip(sampled_results, sampled_params)):
                final_burden, _alive, path, _traj, extinct = res
                sampled_path = [step[0] for step in path]

                cells = [f"{drug:<{drug_col_w}}" for drug in sampled_path]
                while len(cells) < max_cycles:
                    cells.append(f"{'—':<{drug_col_w}}")

                r_str   = f"{samp['r']:.4f}"
                k_str   = f"{samp['K']:.3e}"
                b_str   = f"{samp['beta_dose']:.4g}"
                outcome = "EXTINCT" if extinct else "alive"
                burden_str = f"{final_burden:.3e}"

                param_cells = (
                    f"{r_str:<{param_col_w}} | "
                    f"{k_str:<{param_col_w}} | "
                    f"{b_str:<{param_col_w}} | "
                    f"{burden_str:<{burden_col_w}} | "
                    f"{outcome:<16}"
                )
                print(f"{'Run ' + str(run_idx + 1):<6} | " + " | ".join(cells) + " | " + param_cells)

                row: dict = {"run": run_idx + 1}
                for c in range(max_cycles):
                    row[f"cycle_{c + 1}"] = sampled_path[c] if c < len(sampled_path) else ""
                row["r"]          = samp["r"]
                row["K"]          = samp["K"]
                row["beta_dose"]  = samp["beta_dose"]
                row["end_burden"] = final_burden
                row["outcome"]    = outcome
                path_rows.append(row)

            print("=" * len(header))

            _path_csv   = f"{RESULTS_DIR}/sampled_run_paths.csv"
            _cycle_cols = [f"cycle_{c + 1}" for c in range(max_cycles)]
            _fieldnames = ["run"] + _cycle_cols + ["r", "K", "beta_dose", "end_burden", "outcome"]
            with open(_path_csv, "w", newline="") as _fh:
                _writer = csv.DictWriter(_fh, fieldnames=_fieldnames)
                _writer.writeheader()
                _writer.writerows(path_rows)
            print(f"  Saved: {_path_csv}")

        print(f"Total time: {time() - start_time:.2f}s")

        # ── Baseline path simulation & plotting ────────────────────────────────────
        print("\nRunning baseline path simulation for plotting...")

        baseline_burden_timeline: list[tuple[float, float]] = []
        alive_state = map_start_alive
        day         = baseline_start_day
        baseline_burden_timeline.append((day, alive_state))

        if baseline_res is not None:
            for cycle_idx, (drug, n_seg_steps, cycle_len) in enumerate(baseline_res[2]):
                overrides = _pk_overrides_for(
                    drug, pk_state_map,
                    dose_mg_kg=DOSE_REFERENCE_MG_KG if drug != "none" else None,
                )
                new_alive, seg_traj = simulate_next_state(
                    alive_state, drug,
                    r_map, K_map, k_kill_map, k_clear_map,
                    k_cyto_map, beta_dose_map, mu_base_map, mu_conf_map,
                    pk_overrides=overrides, cycle_len=cycle_len,
                )
                n_tp = len(seg_traj)
                for t_idx, a_val in enumerate(seg_traj):
                    t = day + (t_idx + 1) * cycle_len / n_tp
                    baseline_burden_timeline.append((t, float(a_val)))
                day        += cycle_len
                alive_state = new_alive

        burden_days   = np.array([t for t, _ in baseline_burden_timeline])
        burden_values = np.array([b for _, b in baseline_burden_timeline])

        drug_colors: dict[str, str] = {
            "gemcitabine": "orange",
            "none":        "yellow",
        }

        def _add_drug_shading(ax, path_info, start_day: float = 0.0) -> None:
            if not path_info:
                return
            current_time  = float(start_day)
            shaded_labels: set[str] = set()
            for drug_name, _seg_len, duration in path_info:
                end_time = current_time + duration
                ax.axvspan(
                    current_time, end_time,
                    color=drug_colors.get(drug_name, "gray"),
                    linewidth=0, alpha=0.15,
                    label=drug_name if drug_name not in shaded_labels else None,
                )
                shaded_labels.add(drug_name)
                current_time = end_time

        _baseline_path_info = baseline_res[2] if baseline_res is not None else []

        fig1, ax1 = plt.subplots(figsize=(10, 5))
        _add_drug_shading(ax1, _baseline_path_info, start_day=baseline_start_day)
        ax1.plot(burden_days, burden_values, color="steelblue",
                 linewidth=1.8, label="Predicted (baseline path)")

        if PLOT_OBSERVED_DATA:
            obs_days = sorted(OBSERVED_TUMOR_BURDENS.keys())
            obs_vals = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days]
            ax1.scatter(obs_days, obs_vals, color="firebrick", zorder=5,
                        label="Observed", s=50)

        ax1.set_yscale("log")
        ax1.set_xlabel("Day")
        ax1.set_ylabel("Alive cells (log scale)")
        ax1.set_title(f"Tumor Burden — Baseline Path\n({SAMPLE_NAME})")
        ax1.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1))
        ax1.grid(True, alpha=0.3)
        fig1.tight_layout()
        fig1.savefig(f"{RESULTS_DIR}/baseline_tumor_burden.png", dpi=150, bbox_inches="tight")
        print(f"  Saved: {RESULTS_DIR}/baseline_tumor_burden.png")

        # ── Observed-drug fit plot (fitted model on the drugs actually given) ──────
        plot_observed_fit(
            RESULTS_DIR, SAMPLE_NAME,
            map_params={
                "r": r_map, "K": K_map, "k_kill": k_kill_map, "k_clear": k_clear_map,
                "k_cyto": k_cyto_map, "beta_dose": beta_dose_map,
                "mu_base": mu_base_map, "mu_conf": mu_conf_map,
            },
            pk_state=pk_state_map,
            observed_schedule=OBSERVED_DRUGS_ADMINISTERED,
            observed_burdens=OBSERVED_TUMOR_BURDENS,
        )

        # ── Gibbs fitting animation ────────────────────────────────────────────────
        print("\nGenerating Gibbs fitting animation...")

        N_ANIM_FRAMES = min(200, len(all_iter_params))
        frame_indices = np.linspace(0, len(all_iter_params) - 1,
                                    N_ANIM_FRAMES, dtype=int)

        anim_timelines: list[tuple[np.ndarray, np.ndarray, dict]] = []
        for fi, idx in enumerate(frame_indices):
            p  = all_iter_params[idx]
            tl = _simulate_burden_timeline(
                OBSERVED_DRUGS_ADMINISTERED,
                p["r"], p["K"], p["k_kill"], p["k_clear"],
                p["k_cyto"], p["beta_dose"], p["mu_base"], p["mu_conf"],
                pk_state=p["pk"],
            )
            if tl is not None:
                days = np.array([t for t, _ in tl])
                vals = np.array([b for _, b in tl])
                anim_timelines.append((days, vals, p))
            else:
                anim_timelines.append((np.array([]), np.array([]), p))
            if (fi + 1) % 50 == 0 or fi == len(frame_indices) - 1:
                print(f"  Simulated frame {fi + 1}/{N_ANIM_FRAMES}")

        fig_anim, ax_anim = plt.subplots(figsize=(12, 6))
        fig_anim.subplots_adjust(top=0.85, right=0.82)

        obs_days = sorted(OBSERVED_TUMOR_BURDENS.keys())
        obs_vals = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days]

        ax_anim.set_yscale("log")
        ax_anim.set_xlabel("Day")
        ax_anim.set_ylabel("Alive cells (log scale)")
        ax_anim.grid(True, alpha=0.3)

        # Tune the y-axis to the CONVERGED fit (final ~10% of frames) plus the
        # observed points, instead of padding the union of every burn-in frame
        # by 0.1x-10x.  The old padding buried the per-dose delay oscillations
        # (~0.2 decade) under ~2 decades of empty space, so the fit looked flat.
        _Y_PAD = 10 ** 0.35   # ~0.35 decade of headroom on each side
        _tail  = anim_timelines[max(0, len(anim_timelines)
                                    - max(1, len(anim_timelines) // 10)):]
        _conv  = np.array([v for _, vals, _ in _tail for v in vals if v > 0],
                          dtype=float)
        if _conv.size and obs_vals:
            _lo  = min(float(np.percentile(_conv, 1)), min(obs_vals))
            _hi  = max(float(np.percentile(_conv, 99)), max(obs_vals))
            y_lo, y_hi = _lo / _Y_PAD, _hi * _Y_PAD
        else:
            y_lo, y_hi = 1e3, 1e12
        all_day_vals = [d for days, _, _ in anim_timelines for d in days]
        x_hi = (max(max(all_day_vals), max(obs_days)) * 1.05
                if all_day_vals else max(obs_days) * 1.05)
        ax_anim.set_xlim(-1, x_hi)
        ax_anim.set_ylim(y_lo, y_hi)
        ax_anim.scatter(obs_days, obs_vals, color="firebrick", zorder=5,
                        s=50, label="Observed")

        fit_line, = ax_anim.plot([], [], color="steelblue", linewidth=2.0,
                                  label="Predicted (current iter)")
        title_text  = ax_anim.set_title("", fontsize=11, pad=12)
        ax_anim.legend(loc="upper left", bbox_to_anchor=(1.02, 1))
        ghost_artists = []

        def _anim_init():
            fit_line.set_data([], [])
            title_text.set_text("")
            return [fit_line, title_text]

        def _anim_update(frame_num):
            days, vals, p = anim_timelines[frame_num]
            phase = "burn-in" if p["burnin"] else "sampling"
            if frame_num > 0:
                prev_days, prev_vals, _ = anim_timelines[frame_num - 1]
                if len(prev_days) > 0:
                    gl, = ax_anim.plot(prev_days, prev_vals,
                                       color="steelblue", alpha=0.08,
                                       linewidth=0.6, zorder=1)
                    ghost_artists.append(gl)
            if len(days) > 0:
                fit_line.set_data(days, vals)
            else:
                fit_line.set_data([], [])
            title_text.set_text(
                f"Gibbs Fitting — Iteration {p['iter'] + 1} [{phase}]\n"
                f"r={p['r']:.4f}  K={p['K']:.2e}  "
                f"k_kill={p['k_kill']:.3e}  β={p['beta_dose']:.4g}  logP={p['logP']:.1f}"
            )
            return [fit_line, title_text] + ghost_artists

        anim = animation.FuncAnimation(
            fig_anim, _anim_update, init_func=_anim_init,
            frames=len(anim_timelines), interval=80, blit=False,
        )
        anim.save(f"{RESULTS_DIR}/gibbs_fitting_animation.gif", writer="pillow", dpi=100)
        print(f"  Saved: {RESULTS_DIR}/gibbs_fitting_animation.gif")
        plt.close(fig_anim)

        plt.show()