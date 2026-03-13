from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import os
from time import time

from Missegregation_Model import (
    simulate_karyotype_ode_piecewise,
    get_concentration_curve,
    build_times_with_doses,
    f,
)

# =============================================================================
# Excel data loader
# =============================================================================

EXCEL_PATH = "../data/InVivoData_Gemcitabine/dt_Gem_VT_20260209_v5.xlsx"

# Reference dose for C_peak scaling.
# C_peak in PK_PARAMS_TO_FIT is calibrated at this dose (mg/kg).
# For any other dose d:  C_peak_effective = C_peak_ref * (d / DOSE_REFERENCE_MG_KG).
# Dose 0  →  C_peak = 0  →  no drug effect (functionally equivalent to "none").
DOSE_REFERENCE_MG_KG: float = 120.0


def load_harvest_data(
    excel_path: str,
    harvest_name: str,
) -> tuple[dict[int, float], list[tuple[float, float, str, float]], str]:
    """Load one harvest row from the Excel sheet.

    Returns
    -------
    burdens_cm3
        Mapping of {day: tumour_volume_cm3}.  Day 0 is always set to 0.1
        (the implanted volume); subsequent zero measurements (below the
        detection limit) are skipped until the first truly detectable value.
    schedule
        List of 4-tuples (start_day, end_day, drug, dose_mg_kg).
        All treatment windows use drug="gemcitabine"; the pre-treatment window
        uses drug="none" with dose 0.0.  Because every mouse receives
        gemcitabine (just at different concentrations), dose=0 is handled by
        scaling C_peak to zero rather than switching drug identity.
    ploidy_name
        CBS filename stem for this harvest.
    """
    df = pd.read_excel(excel_path)
    matches = df[df["harvest"] == harvest_name]
    if matches.empty:
        raise ValueError(f"No rows found with harvest == '{harvest_name}'")

    row  = matches.iloc[0]
    cols = list(df.columns)

    start_idx = cols.index("Day_0")
    day_cols  = cols[start_idx:]
    values    = row[day_cols].dropna()

    # Day 0 is the implantation day.  The sheet records 0 (below detection
    # limit) but the injected volume is always 0.1 cm³.
    # For subsequent days: skip any leading zeros (undetectable tumour) and
    # only begin recording once the first positive measurement appears.
    burdens_cm3: dict[int, float] = {0: 0.1}
    found_nonzero = False
    for col, val in values.items():
        day = int(col.split("_")[1])
        if day == 0:
            continue                       # already set above
        fval = float(val)
        if not found_nonzero and fval == 0.0:
            continue                       # skip undetectable leading zeros
        found_nonzero = True
        burdens_cm3[day] = fval

    dose_mg_kg = float(row["Dose"])
    last_day   = float(max(burdens_cm3.keys()))
    first_tx   = float(row["Day of 1st treatment"])

    # All mice are treated with gemcitabine; dose_mg_kg encodes the
    # concentration.  Dose 0 means zero peak plasma concentration (no effect).
    schedule: list[tuple[float, float, str, float]] = [
        (0.0,      first_tx, "none",        0.0),
        (first_tx, last_day, "gemcitabine", dose_mg_kg),
    ]

    ploidy_name = str(row["harvest"]) + ".sps.cbs"

    print(f"  Loaded harvest  : {harvest_name}")
    print(f"  Dose (mg/kg)    : {dose_mg_kg}")
    print(f"  Days with data  : {sorted(burdens_cm3.keys())}")
    print(f"  Drug schedule   : {schedule}")
    print(f"  Ploidy CBS name : {ploidy_name}")

    return burdens_cm3, schedule, ploidy_name


SAMPLE_NAME = "SUM159-2N-0-O_harvest"

_OBSERVED_TUMOR_BURDENS_CM3, OBSERVED_DRUGS_ADMINISTERED, PLOIDY_SAMPLE_NAME = (
    load_harvest_data(EXCEL_PATH, SAMPLE_NAME)
)

# =============================================================================
# Drug / model constants
# =============================================================================

DRUGS = [
    "gemcitabine",
    "bay1895344",
    "alisertib",
    "ispinesib",
    "none",
]

MIN_SIZE    = 1e5
MAX_SIZE    = 2e10
DEFAULT_LEN = 7.0

CYCLE_LENGTHS = {
    "gemcitabine":    28.0,
    "bay1895344":      7.0,
    "alisertib":      21.0,
    "ispinesib":       7.0,
    "volasertib":      7.0,
    "cytarabine":      7.0,
    "umi-77":          7.0,
    "navitoclax":      7.0,
    "none":            7.0,
    "abt-199":        28.0,
    "abt-263":        28.0,
    "capecitabine":   21.0,
    "ceralasertib":   28.0,
    "osi-027":        28.0,
    "adavosertib":    28.0,
    "tegafur":        28.0,
    "tas":            28.0,
    "5-azacytidine":  28.0,
}

BASE_BEAM_WIDTH    = 10
BASE_MAX_DEPTH     = 100
SAMPLED_BEAM_WIDTH = 10
SAMPLED_MAX_DEPTH  = 100
N_SAMPLED_RUNS     = 20

ODE_STEP_FINE   = 0.05
ODE_STEP_COARSE = 0.20

_CELLS_PER_CM3 = 1e7

OBSERVED_TUMOR_BURDENS = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}

HAPLOID_N: int = 23

PLOIDY_TSV_PATH: str = "../data/InVivoData_Gemcitabine/all_ploidy.tsv"


def load_ploidy_distribution(
        tsv_path: str = PLOIDY_TSV_PATH,
        sample_name: str = PLOIDY_SAMPLE_NAME,
) -> np.ndarray:
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except FileNotFoundError:
        print(f"  WARNING: ploidy TSV not found at '{tsv_path}'. "
              "OBSERVED_END_PLOIDY_DISTRIBUTION will be empty — "
              "biopsy likelihood term disabled.")
        return np.array([], dtype=float)

    mask   = df["file"].str.contains(sample_name, na=False)
    values = df.loc[mask, "ploidy"].to_numpy(dtype=float)

    if values.size == 0:
        print(f"  WARNING: no rows matched sample '{sample_name}' in "
              f"'{tsv_path}'. Biopsy likelihood term disabled.")
    else:
        print(f"  Loaded {values.size} ploidy values for '{sample_name}' "
              f"(mean={values.mean():.4f}, std={values.std():.4f})")

    return values


OBSERVED_END_PLOIDY_DISTRIBUTION: np.ndarray = load_ploidy_distribution()

PLOIDY_SIGMA_CHR:         float = 0.5
PLOIDY_LIKELIHOOD_WEIGHT: float = 1.0

START_BEAM_FROM_OBSERVED_END = False

INITIAL_PLOIDY = {46: 0.94e6, 69: 0.05e6, 92: 0.01e6}

R_BASE_PRIOR_MEAN    = 0.28
R_BASE_PRIOR_STD     = 0.05
K_CAP_PRIOR_LOG_MEAN = np.log(6e10)
K_CAP_PRIOR_LOG_STD  = 0.8
BETA_INIT            = 0.05
BETA_PRIOR_LOG_MEAN  = np.log(BETA_INIT)
BETA_PRIOR_LOG_STD   = 0.8

# ── Initial (pre-adaptation) step sizes ──────────────────────────────────────
GIBBS_R_STEP        = 0.05
GIBBS_K_LOG_STEP    = 1.8
GIBBS_BETA_LOG_STEP = 0.30

# ── Adaptive Metropolis hyperparameters ──────────────────────────────────────
# AM kicks in after AM_ADAPT_START iterations and re-tunes every
# AM_ADAPT_INTERVAL iterations thereafter (during burn-in only).
# Proposal std  ←  AM_SCALE * std(chain_history) + AM_EPSILON
# AM_SCALE = 2.38 is the Gelman-Roberts-Gilks (1996) optimum for d = 1.
AM_ADAPT_START    = 100    # iterations before AM begins
AM_ADAPT_INTERVAL = 50     # how often (in iters) to recompute step sizes
AM_SCALE          = 2.38   # optimal for scalar MH (d = 1)
AM_EPSILON        = 1e-6   # floor: prevents step from collapsing to zero

LIKELIHOOD_SIGMA = 0.35
N_GIBBS_SAMPLES  = 200
GIBBS_BURNIN     = 100

PK_PARAMS_TO_FIT: dict[str, dict[str, dict]] = {
    "gemcitabine": {
        # C_peak here is the reference peak plasma concentration at
        # DOSE_REFERENCE_MG_KG (120 mg/kg).  For any other dose the
        # effective C_peak is rescaled automatically in _pk_overrides_for.
        "C_peak":    {"init": 0.032, "prior_log_mean": np.log(0.032), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 0.05,  "prior_log_mean": np.log(0.05),  "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 7.0,   "prior_log_mean": np.log(7.0),   "prior_log_std": 0.3, "step": 0.15},
    },
    "bay1895344": {
        "C_peak":    {"init": 0.5,  "prior_log_mean": np.log(0.5),  "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 0.50, "prior_log_mean": np.log(0.50), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 0.5,  "prior_log_mean": np.log(0.5),  "prior_log_std": 0.3, "step": 0.15},
    },
    "alisertib": {
        "C_peak":    {"init": 1.53, "prior_log_mean": np.log(1.53), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 19.0, "prior_log_mean": np.log(19.0), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 7.0,  "prior_log_mean": np.log(7.0),  "prior_log_std": 0.3, "step": 0.15},
    },
    "ispinesib": {
        "C_peak":    {"init": 0.09, "prior_log_mean": np.log(0.09), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 1.04, "prior_log_mean": np.log(1.04), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 7.0,  "prior_log_mean": np.log(7.0),  "prior_log_std": 0.3, "step": 0.15},
    },
    "volasertib": {
        "C_peak":    {"init": 1.0, "prior_log_mean": np.log(1.0), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 4.0, "prior_log_mean": np.log(4.0), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 7.0, "prior_log_mean": np.log(7.0), "prior_log_std": 0.3, "step": 0.15},
    },
    "cytarabine": {
        "C_peak":    {"init": 1.0, "prior_log_mean": np.log(1.0), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 0.2, "prior_log_mean": np.log(0.2), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 3.5, "prior_log_mean": np.log(3.5), "prior_log_std": 0.3, "step": 0.15},
    },
    "umi-77": {
        "C_peak":    {"init": 1.0, "prior_log_mean": np.log(1.0), "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 0.8, "prior_log_mean": np.log(0.8), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 7.0, "prior_log_mean": np.log(7.0), "prior_log_std": 0.3, "step": 0.15},
    },
    "navitoclax": {
        "C_peak":    {"init": 1.0,  "prior_log_mean": np.log(1.0),  "prior_log_std": 1.0, "step": 0.30},
        "half_life": {"init": 0.73, "prior_log_mean": np.log(0.73), "prior_log_std": 1.0, "step": 0.30},
        "period":    {"init": 1.0,  "prior_log_mean": np.log(1.0),  "prior_log_std": 0.3, "step": 0.15},
    },
    "abt-199":      {"dose":  {"init": 100.0, "prior_log_mean": np.log(100.0), "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.3,   "prior_log_mean": np.log(0.3),   "prior_log_std": 1.0, "step": 0.30}},
    "abt-263":      {"dose":  {"init": 100.0, "prior_log_mean": np.log(100.0), "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.5,   "prior_log_mean": np.log(0.5),   "prior_log_std": 1.0, "step": 0.30}},
    "capecitabine": {"dose":  {"init": 100.0, "prior_log_mean": np.log(100.0), "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.6,   "prior_log_mean": np.log(0.6),   "prior_log_std": 1.0, "step": 0.30}},
    "ceralasertib": {"dose":  {"init": 80.0,  "prior_log_mean": np.log(80.0),  "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.5,   "prior_log_mean": np.log(0.5),   "prior_log_std": 1.0, "step": 0.30}},
    "osi-027":      {"dose":  {"init": 50.0,  "prior_log_mean": np.log(50.0),  "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.5,   "prior_log_mean": np.log(0.5),   "prior_log_std": 1.0, "step": 0.30}},
    "adavosertib":  {"dose":  {"init": 100.0, "prior_log_mean": np.log(100.0), "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.6,   "prior_log_mean": np.log(0.6),   "prior_log_std": 1.0, "step": 0.30}},
    "tegafur":      {"dose":  {"init": 40.0,  "prior_log_mean": np.log(40.0),  "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.5,   "prior_log_mean": np.log(0.5),   "prior_log_std": 1.0, "step": 0.30}},
    "tas":          {"dose":  {"init": 60.0,  "prior_log_mean": np.log(60.0),  "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 0.7,   "prior_log_mean": np.log(0.7),   "prior_log_std": 1.0, "step": 0.30}},
    "5-azacytidine":{"dose":  {"init": 100.0, "prior_log_mean": np.log(100.0), "prior_log_std": 1.0, "step": 0.30},
                     "ke_day":{"init": 2.0,   "prior_log_mean": np.log(2.0),   "prior_log_std": 1.0, "step": 0.30}},
}


# =============================================================================
# Model helpers
# =============================================================================

def get_cycle_length(drug: str) -> float:
    return CYCLE_LENGTHS.get(drug, DEFAULT_LEN)


def _pk_overrides_for(
    drug: str,
    pk_state: dict[str, dict] | None,
    dose_mg_kg: float | None = None,
) -> dict:
    """Build PK override kwargs for *drug*.

    If *dose_mg_kg* is provided and the drug is "gemcitabine", C_peak is
    scaled linearly by dose_mg_kg / DOSE_REFERENCE_MG_KG.  This means:
      - dose 120 mg/kg  →  C_peak unchanged   (full effect)
      - dose  60 mg/kg  →  C_peak × 0.5
      - dose  30 mg/kg  →  C_peak × 0.25
      - dose   0 mg/kg  →  C_peak = 0          (no drug effect)
    """
    overrides: dict = dict(pk_state.get(drug, {})) if pk_state else {}
    if dose_mg_kg is not None and drug == "gemcitabine" and DOSE_REFERENCE_MG_KG > 0:
        scale    = dose_mg_kg / DOSE_REFERENCE_MG_KG
        base     = overrides.get(
            "C_peak", PK_PARAMS_TO_FIT["gemcitabine"]["C_peak"]["init"]
        )
        overrides["C_peak"] = base * scale
    return overrides


def _fill_gaps_with_none(
        schedule: list[tuple[float, float, str, float]],
) -> list[tuple[float, float, str, float]]:
    """Insert ("none", 0.0) segments to cover any gaps in *schedule*."""
    if not schedule:
        return []
    filled: list[tuple[float, float, str, float]] = []
    for start, end, drug, dose in schedule:
        if filled and filled[-1][1] < start:
            filled.append((filled[-1][1], start, "none", 0.0))
        filled.append((start, end, drug, dose))
    return filled


def _observed_drug_names() -> list[str]:
    return [drug for _, _, drug, _ in OBSERVED_DRUGS_ADMINISTERED]


def _observed_drug_set() -> set[str]:
    return {drug.lower() for _, _, drug, _ in OBSERVED_DRUGS_ADMINISTERED}


def _run_ode(ploidy_status: dict, drug: str,
             r_base: float, k_cap: float,
             pk_overrides: dict,
             cycle_len: float,
             dt: float,
             beta: float = BETA_INIT) -> tuple[dict, np.ndarray]:
    C_fn  = get_concentration_curve(drug, **pk_overrides)
    TIMES = build_times_with_doses(
        (0.0, cycle_len), dt,
        drug_name=drug, include_days=True, eps=1e-8,
    )
    _t, Ns, T_mat, _T_total, _M = simulate_karyotype_ode_piecewise(
        ploidy_status, drug,
        t_span=(0.0, cycle_len),
        r=r_base, Kcap=k_cap, beta=beta,
        N_min=10, N_max=90,
        C_fn=C_fn, f_param_fn=f, t_eval=TIMES, max_step=dt,
        renormalize_M=True,
    )
    final_counts = T_mat[:, -1]
    new_status   = {int(N): float(c) for N, c in zip(Ns, final_counts) if c > 0}
    trajectory   = T_mat.T[1:]
    return new_status, trajectory


def simulate_next_state(ploidy_status, drug, r_base, k_cap=6e10,
                        pk_overrides=None, cycle_len=None, beta=BETA_INIT):
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    return _run_ode(ploidy_status, drug, r_base, k_cap, overrides, T,
                    dt=ODE_STEP_FINE, beta=beta)


def simulate_next_state_cheap(ploidy_status, drug, r_base, k_cap,
                               pk_overrides=None, cycle_len=None, beta=BETA_INIT):
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    return _run_ode(ploidy_status, drug, r_base, k_cap, overrides, T,
                    dt=ODE_STEP_COARSE, beta=beta)


def get_observed_end_ploidy(r_base, k_cap, pk_state=None, beta=BETA_INIT):
    if not OBSERVED_DRUGS_ADMINISTERED:
        return dict(INITIAL_PLOIDY)
    ploidy   = dict(INITIAL_PLOIDY)
    schedule = _fill_gaps_with_none(OBSERVED_DRUGS_ADMINISTERED)
    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            ploidy, _ = simulate_next_state(ploidy, drug, r_base, k_cap,
                                            pk_overrides=overrides,
                                            cycle_len=cycle_len, beta=beta)
        except Exception:
            return None
        total = sum(ploidy.values())
        if total < MIN_SIZE or not np.isfinite(total):
            return None
    return ploidy


def _log_likelihood_ploidy(predicted_ploidy, observed_ploidy_values,
                            sigma=PLOIDY_SIGMA_CHR,
                            weight=PLOIDY_LIKELIHOOD_WEIGHT):
    total = sum(predicted_ploidy.values())
    if total <= 0:
        return -np.inf
    chr_bins  = np.array(list(predicted_ploidy.keys()),   dtype=float)
    log_props = np.log(np.array(list(predicted_ploidy.values()), dtype=float) / total)
    obs_chr   = observed_ploidy_values * HAPLOID_N
    diff      = obs_chr[:, None] - chr_bins[None, :]
    log_phi   = -0.5 * (diff / sigma) ** 2 - np.log(sigma * np.sqrt(2 * np.pi))
    log_mix   = log_props[None, :] + log_phi
    log_probs = log_mix.max(axis=1) + np.log(
        np.exp(log_mix - log_mix.max(axis=1, keepdims=True)).sum(axis=1))
    return float(weight * log_probs.sum())


def _simulate_burden_timeline(observed_schedule, r_base, k_cap,
                               pk_state=None, beta=BETA_INIT):
    ploidy   = dict(INITIAL_PLOIDY)
    timeline = [(0.0, float(sum(ploidy.values())))]
    schedule = _fill_gaps_with_none(observed_schedule)
    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            new_ploidy, seg_traj = simulate_next_state_cheap(
                ploidy, drug, r_base, k_cap,
                pk_overrides=overrides, cycle_len=cycle_len, beta=beta)
        except Exception:
            return None
        n_tp = len(seg_traj)
        if n_tp > 0:
            for t_idx, burden_by_ploidy in enumerate(seg_traj):
                day   = start_day + (t_idx + 1) * cycle_len / n_tp
                total = float(np.asarray(burden_by_ploidy).sum())
                timeline.append((day, total))
        ploidy = new_ploidy
    return timeline


def _log_likelihood(r_base, k_cap, observed_schedule, observed_burdens,
                    pk_state=None, beta=BETA_INIT, observed_end_ploidy=None):
    if not (0.0 < r_base < 1.0) or k_cap <= 0 or beta <= 0:
        return -np.inf
    timeline = _simulate_burden_timeline(observed_schedule, r_base, k_cap,
                                         pk_state, beta)
    if timeline is None:
        return -np.inf
    sim_days    = np.array([t[0] for t in timeline])
    sim_burdens = np.array([t[1] for t in timeline])
    log_lik = 0.0
    for obs_day, obs_burden in observed_burdens.items():
        if obs_burden <= 0:
            continue
        idx         = int(np.argmin(np.abs(sim_days - obs_day)))
        pred_burden = sim_burdens[idx]
        if pred_burden <= 0:
            return -np.inf
        log_ratio = np.log(obs_burden) - np.log(pred_burden)
        log_lik  += -0.5 * (log_ratio / LIKELIHOOD_SIGMA) ** 2
    if observed_end_ploidy is not None and len(observed_end_ploidy) > 0:
        predicted_end_ploidy = get_observed_end_ploidy(r_base, k_cap,
                                                       pk_state, beta)
        if predicted_end_ploidy is None:
            return -np.inf
        log_lik += _log_likelihood_ploidy(predicted_end_ploidy, observed_end_ploidy)
    return log_lik


def _log_prior(r_base, k_cap, beta, pk_state=None):
    if beta <= 0:
        return -np.inf
    log_p  = -0.5 * ((r_base - R_BASE_PRIOR_MEAN) / R_BASE_PRIOR_STD) ** 2
    log_p += -0.5 * ((np.log(k_cap) - K_CAP_PRIOR_LOG_MEAN) / K_CAP_PRIOR_LOG_STD) ** 2
    log_p += -0.5 * ((np.log(beta)  - BETA_PRIOR_LOG_MEAN)  / BETA_PRIOR_LOG_STD)  ** 2
    if pk_state:
        for drug, params in pk_state.items():
            spec = PK_PARAMS_TO_FIT.get(drug, {})
            for param, val in params.items():
                if val <= 0:
                    return -np.inf
                cfg   = spec[param]
                log_p += -0.5 * ((np.log(val) - cfg["prior_log_mean"]) / cfg["prior_log_std"]) ** 2
    return log_p


def _log_posterior(r_base, k_cap, beta, observed_schedule, observed_burdens,
                   pk_state=None, observed_end_ploidy=None):
    lp = _log_prior(r_base, k_cap, beta, pk_state)
    if not np.isfinite(lp):
        return -np.inf
    return lp + _log_likelihood(r_base, k_cap, observed_schedule,
                                observed_burdens, pk_state, beta,
                                observed_end_ploidy)


# =============================================================================
# Adaptive Metropolis helper
# =============================================================================

def _am_step(history: list[float], fallback: float) -> float:
    """Return an AM-tuned proposal std from scalar chain history.

    Uses the Gelman-Roberts-Gilks (1996) formula for d = 1:
        σ_prop = AM_SCALE * std(history) + AM_EPSILON

    Falls back to *fallback* when the history is too short or has zero variance.
    """
    if len(history) < AM_ADAPT_START:
        return fallback
    std = float(np.std(history))
    proposed = AM_SCALE * std + AM_EPSILON
    # Sanity-clamp: never let the step shrink below 1 % of the fallback value
    return max(proposed, 0.01 * fallback)


# =============================================================================
# Gibbs sampler with Adaptive Metropolis proposals
# =============================================================================

def run_gibbs_sampler(observed_schedule, observed_burdens,
                      observed_end_ploidy=None,
                      n_samples=N_GIBBS_SAMPLES, burnin=GIBBS_BURNIN):
    """Metropolis-within-Gibbs with Adaptive Metropolis (AM) proposals.

    How AM works here
    -----------------
    Each scalar parameter keeps a running history of its chain values (on the
    natural scale for r_base, on the log scale for k_cap / beta / PK params).
    Every AM_ADAPT_INTERVAL iterations *during burn-in*, the proposal std for
    that parameter is recomputed as:

        σ = AM_SCALE * std(history) + AM_EPSILON          (Haario et al. 2001)

    AM_SCALE = 2.38 is the Gelman-Roberts-Gilks (1996) optimal constant for
    d = 1 scalar updates. Adaptation is frozen once sampling begins so that
    the post-burn-in chain is a valid MCMC sample from the fixed posterior.
    """
    observed_drug_set = {d.lower() for _, _, d, _ in observed_schedule}
    active_pk = {
        drug: {p: cfg["init"] for p, cfg in params.items()}
        for drug, params in PK_PARAMS_TO_FIT.items()
        if drug in observed_drug_set
    }
    n_pk = sum(len(v) for v in active_pk.values())
    use_ploidy_data = (observed_end_ploidy is not None
                       and len(observed_end_ploidy) > 0)
    print(f"--- Gibbs Sampler + AM: R_BASE + K_CAP + beta"
          + (f" + {n_pk} PK param(s)" if n_pk else "")
          + (" + biopsy ploidy distribution" if use_ploidy_data else "") + " ---")
    print(f"    AM_ADAPT_START={AM_ADAPT_START}  "
          f"AM_ADAPT_INTERVAL={AM_ADAPT_INTERVAL}  "
          f"AM_SCALE={AM_SCALE}  AM_EPSILON={AM_EPSILON}")

    rng  = np.random.default_rng()
    r    = R_BASE_PRIOR_MEAN
    k    = np.exp(K_CAP_PRIOR_LOG_MEAN)
    beta = BETA_INIT
    pk   = {drug: dict(params) for drug, params in active_pk.items()}

    cur_lp = _log_posterior(r, k, beta, observed_schedule, observed_burdens,
                            pk, observed_end_ploidy)

    # ── Current adaptive step sizes (updated during burn-in) ─────────────────
    r_step    = GIBBS_R_STEP
    k_step    = GIBBS_K_LOG_STEP
    beta_step = GIBBS_BETA_LOG_STEP
    pk_steps  = {
        drug: {p: PK_PARAMS_TO_FIT[drug][p]["step"] for p in params}
        for drug, params in active_pk.items()
    }

    # ── Per-parameter chain histories (for AM; stored on proposal scale) ─────
    # r_base: natural scale; k/beta/PK: log scale
    r_hist:    list[float] = [r]
    k_hist:    list[float] = [np.log(k)]
    beta_hist: list[float] = [np.log(beta)]
    pk_hist: dict[str, dict[str, list[float]]] = {
        drug: {p: [np.log(params[p])] for p in params}
        for drug, params in active_pk.items()
    }

    samples, log_posts = [], []
    accept: dict[str, int] = {"r": 0, "k": 0, "beta": 0}
    for drug, params in pk.items():
        for p in params:
            accept[f"{drug}.{p}"] = 0

    total = n_samples + burnin

    for i in range(total):
        in_burnin = i < burnin

        # ── Adapt step sizes every AM_ADAPT_INTERVAL iters during burn-in ────
        if in_burnin and i > 0 and (i % AM_ADAPT_INTERVAL == 0):
            r_step    = _am_step(r_hist,    GIBBS_R_STEP)
            k_step    = _am_step(k_hist,    GIBBS_K_LOG_STEP)
            beta_step = _am_step(beta_hist, GIBBS_BETA_LOG_STEP)
            for drug, params in pk_steps.items():
                for p in params:
                    pk_steps[drug][p] = _am_step(
                        pk_hist[drug][p],
                        PK_PARAMS_TO_FIT[drug][p]["step"],
                    )

        if i % 20 == 0:
            phase = "burn-in" if in_burnin else "sampling"
            print(f"  iter {i:>4}/{total}  [{phase}]  "
                  f"R={r:.4f}(σ={r_step:.4f})  "
                  f"K={k:.3e}(σ={k_step:.3f})  "
                  f"beta={beta:.4g}(σ={beta_step:.3f})  "
                  f"logP={cur_lp:.2f}")

        # ── R_BASE ────────────────────────────────────────────────────────────
        r_prop  = r + rng.normal(0.0, r_step)
        prop_lp = _log_posterior(r_prop, k, beta, observed_schedule,
                                 observed_burdens, pk, observed_end_ploidy)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            r, cur_lp = r_prop, prop_lp
            if not in_burnin:
                accept["r"] += 1
        r_hist.append(r)

        # ── K_CAP (log scale) ─────────────────────────────────────────────────
        k_prop  = np.exp(np.log(k) + rng.normal(0.0, k_step))
        prop_lp = _log_posterior(r, k_prop, beta, observed_schedule,
                                 observed_burdens, pk, observed_end_ploidy)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            k, cur_lp = k_prop, prop_lp
            if not in_burnin:
                accept["k"] += 1
        k_hist.append(np.log(k))

        # ── beta (log scale) ──────────────────────────────────────────────────
        beta_prop = np.exp(np.log(beta) + rng.normal(0.0, beta_step))
        prop_lp   = _log_posterior(r, k, beta_prop, observed_schedule,
                                   observed_burdens, pk, observed_end_ploidy)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            beta, cur_lp = beta_prop, prop_lp
            if not in_burnin:
                accept["beta"] += 1
        beta_hist.append(np.log(beta))

        # ── PK params (log scale) ─────────────────────────────────────────────
        for drug, params in pk.items():
            for param in list(params.keys()):
                step     = pk_steps[drug][param]
                val_prop = np.exp(np.log(params[param]) + rng.normal(0.0, step))
                pk_prop  = {d: dict(ps) for d, ps in pk.items()}
                pk_prop[drug][param] = val_prop
                prop_lp = _log_posterior(r, k, beta, observed_schedule,
                                         observed_burdens, pk_prop,
                                         observed_end_ploidy)
                if np.log(rng.uniform()) < prop_lp - cur_lp:
                    pk[drug][param] = val_prop
                    cur_lp = prop_lp
                    if not in_burnin:
                        accept[f"{drug}.{param}"] += 1
                pk_hist[drug][param].append(np.log(pk[drug][param]))

        if not in_burnin:
            samples.append({"r_base": r, "k_cap": k, "beta": beta,
                             "pk": {d: dict(ps) for d, ps in pk.items()}})
            log_posts.append(cur_lp)

    log_posts = np.array(log_posts)

    print(f"\n  Final adaptive step sizes:")
    print(f"    R_BASE : {r_step:.5f}   K_CAP(log) : {k_step:.4f}   "
          f"beta(log) : {beta_step:.4f}")
    for drug, params in pk_steps.items():
        for p, s in params.items():
            print(f"    {drug}.{p} : {s:.4f}")

    print(f"\n  Acceptance rates (sampling phase):")
    print(f"    R_BASE : {accept['r']    / n_samples:.3f}")
    print(f"    K_CAP  : {accept['k']    / n_samples:.3f}")
    print(f"    beta   : {accept['beta'] / n_samples:.3f}")
    for drug, params in pk.items():
        for p in params:
            key = f"{drug}.{p}"
            print(f"    {key:<22} : {accept[key] / n_samples:.3f}")

    r_mean    = np.mean([s["r_base"] for s in samples])
    k_mean    = np.mean([s["k_cap"]  for s in samples])
    beta_mean = np.mean([s["beta"]   for s in samples])
    print(f"\n  Posterior means  R={r_mean:.4f}  K={k_mean:.3e}  "
          f"beta={beta_mean:.4g}")

    return samples, log_posts


def compute_bma_weights(log_posteriors: np.ndarray) -> np.ndarray:
    lp = np.asarray(log_posteriors, dtype=float)
    w  = np.exp(lp - np.max(lp))
    return w / w.sum()


def _simulate_next_state_wrapper(ploidy, drug, path, traj, r_base, k_cap,
                                  pk_overrides, beta):
    new_status, seg_traj = simulate_next_state(ploidy, drug, r_base, k_cap,
                                               pk_overrides=pk_overrides,
                                               beta=beta)
    return new_status, seg_traj, path, traj, drug


def _beam_search_step(current_beams, executor, r_base, k_cap, beam_width,
                      pk_state=None, beta=BETA_INIT):
    futures = []
    for burden, ploidy, path, traj, extinct in current_beams:
        if extinct:
            continue
        for drug in DRUGS:
            # Beam search explores future treatments at full reference dose
            overrides = _pk_overrides_for(drug, pk_state,
                                          dose_mg_kg=DOSE_REFERENCE_MG_KG)
            futures.append(executor.submit(
                _simulate_next_state_wrapper,
                ploidy, drug, path, traj, r_base, k_cap, overrides, beta))
    next_candidates = []
    for future in as_completed(futures):
        next_ploidy, seg_traj, old_path, old_traj, drug = future.result()
        new_burden   = sum(next_ploidy.values())
        segment_info = (drug, len(seg_traj))
        if new_burden < MIN_SIZE:
            next_candidates.append(
                (new_burden, next_ploidy,
                 old_path + [segment_info], old_traj + list(seg_traj), True))
        elif new_burden <= MAX_SIZE:
            next_candidates.append(
                (new_burden, next_ploidy,
                 old_path + [segment_info], old_traj + list(seg_traj), False))
    next_candidates.sort(key=lambda x: x[0])
    return next_candidates[:beam_width]


def run_single_beam_search(run_idx, r_base, k_cap, beam_width, max_depth,
                            start_ploidy=None, pk_state=None, beta=BETA_INIT):
    if start_ploidy is None:
        start_ploidy = dict(INITIAL_PLOIDY)
    initial_burden = sum(start_ploidy.values())
    beam = [(initial_burden, start_ploidy, [],
             [np.array(list(start_ploidy.values()))], False)]
    with ThreadPoolExecutor(max_workers=len(DRUGS) * beam_width) as executor:
        for _ in range(max_depth):
            beam = _beam_search_step(beam, executor, r_base, k_cap,
                                     beam_width, pk_state, beta)
            if not beam or all(b[4] for b in beam):
                break
    return beam[0] if beam else None


def _beam_search_worker(i, r_i, k_i, beam_width, max_depth,
                         use_observed_end, pk_state, beta_i):
    if use_observed_end and OBSERVED_DRUGS_ADMINISTERED:
        sp = get_observed_end_ploidy(r_i, k_i, pk_state, beta=beta_i)
        if sp is None:
            sp = dict(INITIAL_PLOIDY)
    else:
        sp = dict(INITIAL_PLOIDY)
    return run_single_beam_search(i, r_i, k_i, beam_width, max_depth,
                                  start_ploidy=sp, pk_state=pk_state,
                                  beta=beta_i)


if __name__ == "__main__":
    start_time = time()

    posterior_samples, log_posteriors = run_gibbs_sampler(
        OBSERVED_DRUGS_ADMINISTERED,
        OBSERVED_TUMOR_BURDENS,
        observed_end_ploidy=OBSERVED_END_PLOIDY_DISTRIBUTION,
        n_samples=N_GIBBS_SAMPLES,
        burnin=GIBBS_BURNIN,
    )

    bma_weights = compute_bma_weights(log_posteriors)

    map_idx      = int(np.argmax(log_posteriors))
    map_sample   = posterior_samples[map_idx]
    r_base_map   = map_sample["r_base"]
    k_cap_map    = map_sample["k_cap"]
    beta_map     = map_sample["beta"]
    pk_state_map = map_sample["pk"]

    print(f"\nMAP  R={r_base_map:.4f}  K={k_cap_map:.3e}  beta={beta_map:.4g}")

    if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
        map_start_ploidy = get_observed_end_ploidy(r_base_map, k_cap_map,
                                                    pk_state_map, beta=beta_map)
        if map_start_ploidy is None:
            map_start_ploidy = dict(INITIAL_PLOIDY)
    else:
        map_start_ploidy = dict(INITIAL_PLOIDY)

    baseline_res  = run_single_beam_search(
        "baseline", r_base_map, k_cap_map,
        BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
        start_ploidy=map_start_ploidy,
        pk_state=pk_state_map,
        beta=beta_map)
    baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
    print(f"Baseline path: {baseline_path}")

    rng = np.random.default_rng()
    selected_idx     = rng.choice(len(posterior_samples),
                                  size=N_SAMPLED_RUNS, p=bma_weights, replace=True)
    selected_samples = [posterior_samples[i] for i in selected_idx]
    selected_weights = bma_weights[selected_idx]
    selected_weights = selected_weights / selected_weights.sum()

    sampled_results, run_weights = [], []
    use_observed_end = START_BEAM_FROM_OBSERVED_END and bool(OBSERVED_DRUGS_ADMINISTERED)

    # Counts runs where every candidate on the beam exceeded MAX_SIZE at some
    # depth, causing the beam to empty and the worker to return None.  These
    # runs are excluded from the flip-rate denominator and their BMA weight is
    # silently redistributed among the surviving runs.
    n_full_maxout = 0

    with ProcessPoolExecutor(max_workers=min(N_SAMPLED_RUNS, os.cpu_count())) as pool:
        future_map = {
            pool.submit(
                _beam_search_worker, i,
                float(selected_samples[i]["r_base"]),
                float(selected_samples[i]["k_cap"]),
                SAMPLED_BEAM_WIDTH, SAMPLED_MAX_DEPTH,
                use_observed_end,
                selected_samples[i]["pk"],
                float(selected_samples[i]["beta"]),
            ): i
            for i in range(N_SAMPLED_RUNS)
        }
        for future in as_completed(future_map):
            i   = future_map[future]
            res = future.result()
            if res is not None:
                sampled_results.append(res)
                run_weights.append(selected_weights[i])
            else:
                n_full_maxout += 1

    run_weights = np.array(run_weights)
    run_weights /= run_weights.sum()

    print(f"\n  Full maxouts (all candidates exceeded MAX_SIZE={MAX_SIZE:.0e}): "
          f"{n_full_maxout}/{N_SAMPLED_RUNS} runs "
          f"({100 * n_full_maxout / N_SAMPLED_RUNS:.1f}%) — "
          f"excluded from flip rate, their BMA weight redistributed.")

    print("\n" + "=" * 65)
    print(f"{'Cycle':<7} | {'Baseline Drug':<16} | {'Unweighted':>11}")
    print("-" * 65)
    for i in range(len(baseline_path)):
        target_drug     = baseline_path[i]
        unweighted_flip = 0
        active_count    = 0
        for res, w in zip(sampled_results, run_weights):
            sampled_path = [step[0] for step in res[2]]
            if i < len(sampled_path):
                active_count += 1
                if sampled_path[i] != target_drug:
                    unweighted_flip += 1
        raw_rate = (unweighted_flip / active_count) if active_count > 0 else 0.0
        print(f"{i + 1:<7} | {target_drug:<16} | {raw_rate * 100:>9.2f}%")

    print("=" * 65)
    print(f"Total time: {time() - start_time:.2f}s")