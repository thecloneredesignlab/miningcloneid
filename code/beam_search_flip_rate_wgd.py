from __future__ import annotations

import csv
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from time import time
from ploidy_model_wgd_missegg_transit import (
    simulate_karyotype_ode_piecewise,
    get_concentration_curve,
    build_times_with_doses,
    f,
)

# PARAMETERS:
MIN_SIZE    = 1e5
MAX_SIZE    = 2e10
DEFAULT_LEN = 7.0

BASE_BEAM_WIDTH    = 100
BASE_MAX_DEPTH     = 100
SAMPLED_BEAM_WIDTH = 100
SAMPLED_MAX_DEPTH  = 100
N_SAMPLED_RUNS     = 100

ODE_STEP   = 0.05

_CELLS_PER_CM3 = 1e7

R_BASE_FIRST_GUESS    = 0.28
R_BASE_PRIOR_STD     = 0.05
K_CAP_FIRST_LOG_GUESS = np.log(6e10)
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
N_GIBBS_SAMPLES  = 2000
GIBBS_BURNIN     = 1000

HAPLOID_N: int = 23

# =============================================================================
# Drug kinetics config loader
# =============================================================================

_CONFIG_PATH = "../config/drug_kinetics.json"


def _load_drug_kinetics_config(path: str = _CONFIG_PATH) -> dict:
    """Load drug kinetics parameters from the JSON config file.

    The config stores `init`, `prior_log_std`, and `step` for each PK
    parameter.  `prior_log_mean` is derived here as `log(init)` so the
    config stays human-readable (no raw log values to maintain by hand).
    """
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

# =============================================================================
# Excel data loader
# =============================================================================

EXCEL_PATH = "../data/InVivoData_Gemcitabine/dt_Gem_VT_20260209_v5.xlsx"

# Reference dose for C_peak scaling (loaded from config).
# C_peak in PK_PARAMS_TO_FIT is calibrated at this dose (mg/kg).
# For any other dose d:  C_peak_effective = C_peak_ref * (d / DOSE_REFERENCE_MG_KG).
# Dose 0  →  C_peak = 0  →  no drug effect (functionally equivalent to "none").
DOSE_REFERENCE_MG_KG: float = _DK["DOSE_REFERENCE_MG_KG"]


def load_harvest_data(
    excel_path: str,
    harvest_name: str,
    verbose: bool = True,
) -> tuple[dict[int, float], list[tuple[float, float, str, float]], str]:
    """Load one harvest row from the Excel sheet.

    Returns
    -------
    burdens_cm3
        Mapping of {day: tumor_volume_cm3}.  Day 0 is always set to 0.1
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
    # For subsequent days: skip any leading zeros (undetectable tumor) and
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

    if verbose:
        print(f"  Loaded harvest  : {harvest_name}")
        print(f"  Dose (mg/kg)    : {dose_mg_kg}")
        print(f"  Days with data  : {sorted(burdens_cm3.keys())}")
        print(f"  Drug schedule   : {schedule}")
        print(f"  Ploidy CBS name : {ploidy_name}")

    return burdens_cm3, schedule, ploidy_name


SAMPLE_NAME = "SUM159-4N-120-RL_harvest"

_OBSERVED_TUMOR_BURDENS_CM3, OBSERVED_DRUGS_ADMINISTERED, PLOIDY_SAMPLE_NAME = (
    load_harvest_data(EXCEL_PATH, SAMPLE_NAME, verbose=(__name__ == "__main__"))
)

# Day on which the first active treatment begins.  Before this point the beam
# search is restricted to "none" (no drug), mirroring the observed schedule.
FIRST_TX_DAY: float = next(
    (start for start, _end, drug, _dose in OBSERVED_DRUGS_ADMINISTERED
     if drug != "none"),
    0.0,
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

# Loaded from config/drug_kinetics.json
CYCLE_LENGTHS: dict[str, float] = _DK["CYCLE_LENGTHS"]

OBSERVED_TUMOR_BURDENS = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}

PLOIDY_TSV_PATH: str = "../data/InVivoData_Gemcitabine/all_ploidy.tsv"

# =============================================================================
# Injected-cell CBS paths — one file per starting ploidy condition
# =============================================================================

CBS_2N_PATH: str = (
    "../data/InVivoData_Gemcitabine/"
    "2N_InjectedCells_SUM-159_NLS_2N_A7M_K_harvest.sps.cbs"
)
CBS_4N_PATH: str = (
    "../data/InVivoData_Gemcitabine/"
    "4N_InjectedCells_SUM-159_NLS_4N_A5M_K_harvest.sps.cbs"
)


def load_initial_ploidy_from_cbs(
    cbs_path: str,
    total_injected: float = 1e6,
    verbose: bool = True,
) -> dict[int, float]:
    try:
        df = pd.read_csv(cbs_path, sep="\t")
    except FileNotFoundError:
        raise FileNotFoundError(
            f"Injected-cell CBS file not found: '{cbs_path}'. "
            "Update CBS_2N_PATH / CBS_4N_PATH to point to the correct location."
        )

    df_chr = df[df["chr"] != 999]
    additional_dna = df[df["chr"] == 999]

    cell_cols = [c for c in df.columns if c.startswith("SP_")]
    n_cells = len(cell_cols)
    if n_cells == 0:
        raise ValueError(f"No SP_* cell columns found in '{cbs_path}'.")

    # Sum copy numbers across all autosomes for each cell
    chr_sums = df_chr[cell_cols].sum(axis=0)

    # Add extra DNA as a fractional multiplier (e.g., additional_dna=0.1 → multiply by 1.1)
    if len(additional_dna) > 0:
        extra_dna_values = additional_dna[cell_cols].iloc[0]
        chr_sums = chr_sums * (1.0 + extra_dna_values)

    # Build {rounded_chr_count: cell_count} and scale to total_injected
    from collections import Counter
    counts = Counter(int(round(v)) for v in chr_sums.values)
    initial_ploidy = {
        chr_n: (n_cells_at / n_cells) * total_injected
        for chr_n, n_cells_at in sorted(counts.items())
    }

    if verbose:
        print(f"  CBS file        : {cbs_path}")
        print(f"  Cells measured  : {n_cells}")
        for chr_n, cells in initial_ploidy.items():
            print(f"    chr_count={chr_n:3d} → {cells:.0f} injected cells "
                  f"({counts[chr_n]}/{n_cells} measured)")

    return initial_ploidy


def _select_cbs_path(harvest_name: str) -> str:
    """Return the CBS path that matches *harvest_name*'s ploidy condition.

    Harvest names follow the pattern ``SUM159-{2N|4N}-…``.  The string "4N"
    takes precedence over "2N" to avoid false matches.
    """
    name_upper = harvest_name.upper()
    if "4N" in name_upper:
        return CBS_4N_PATH
    if "2N" in name_upper:
        return CBS_2N_PATH
    raise ValueError(
        f"Cannot infer ploidy condition (2N / 4N) from harvest name "
        f"'{harvest_name}'.  Rename the harvest or set INITIAL_PLOIDY manually."
    )


def load_ploidy_distribution(
        tsv_path: str = PLOIDY_TSV_PATH,
        sample_name: str = PLOIDY_SAMPLE_NAME,
        verbose: bool = True,
) -> np.ndarray:
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except FileNotFoundError:
        if verbose:
            print(f"  WARNING: ploidy TSV not found at '{tsv_path}'. "
                  "OBSERVED_END_PLOIDY_DISTRIBUTION will be empty — "
                  "biopsy likelihood term disabled.")
        return np.array([], dtype=float)

    mask   = df["file"].str.contains(sample_name, na=False)
    values = df.loc[mask, "ploidy"].to_numpy(dtype=float)

    if values.size == 0:
        if verbose:
            print(f"  WARNING: no rows matched sample '{sample_name}' in "
                  f"'{tsv_path}'. Biopsy likelihood term disabled.")
    else:
        if verbose:
            print(f"  Loaded {values.size} ploidy values for '{sample_name}' "
                  f"(mean={values.mean():.4f}, std={values.std():.4f})")

    return values


OBSERVED_END_PLOIDY_DISTRIBUTION: np.ndarray = load_ploidy_distribution(
    verbose=(__name__ == "__main__")
)

PLOIDY_SIGMA_CHR:         float = 0.5
PLOIDY_LIKELIHOOD_WEIGHT: float = 1.0

START_BEAM_FROM_OBSERVED_END = False

# =============================================================================
# Plotting flags
# =============================================================================

# Set to True to overlay observed tumor burden data points on Plot 1 and
# the combined burden+ploidy plot (Plot 1b).
PLOT_OBSERVED_DATA: bool = False

# Derive the starting ploidy from the injected-cell CBS file that corresponds
# to this harvest's ploidy condition (2N or 4N).
# Formula: initial_cells[chr_n] = (cells_at_chr_n / n_measured) * 1e6
# where 1e6 is the number of cells injected per mouse.
_CBS_PATH_FOR_HARVEST: str = _select_cbs_path(SAMPLE_NAME)
INITIAL_PLOIDY: dict[int, float] = load_initial_ploidy_from_cbs(
    _CBS_PATH_FOR_HARVEST,
    total_injected=1e6,
    verbose=(__name__ == "__main__"),
)

# Loaded from config/drug_kinetics.json (prior_log_mean derived as log(init))
PK_PARAMS_TO_FIT: dict[str, dict[str, dict]] = _DK["PK_PARAMS_TO_FIT"]


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
                    dt=ODE_STEP, beta=beta)


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
            new_ploidy, seg_traj = simulate_next_state(
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
    # Flat priors: return 0 for all valid parameters
    if not (0.0 < r_base < 1.0) or k_cap <= 0:
        return -np.inf
    if pk_state:
        for drug, params in pk_state.items():
            for param, val in params.items():
                if val <= 0:
                    return -np.inf
    return 0.0


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

    FIX: non-finite values (inf/nan) that can enter the history via exp()
    overflow are filtered out before computing std.  Without this, np.std()
    computes mean=inf then tries inf-inf=nan, triggering a secondary
    RuntimeWarning and returning nan, which collapses the step size to
    AM_EPSILON — effectively freezing the chain.
    """
    if len(history) < AM_ADAPT_START:
        return fallback
    # Filter out non-finite values (inf/nan) before computing std.
    finite_history = [v for v in history if np.isfinite(v)]
    if len(finite_history) < AM_ADAPT_START:
        return fallback
    std = float(np.std(finite_history))
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
    r    = R_BASE_FIRST_GUESS
    k    = np.exp(K_CAP_FIRST_LOG_GUESS)
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

    # Maximum allowed value for log-scale proposals.
    # np.exp(30) ≈ 1.07e13, np.exp(-30) ≈ 9.4e-14 — generous bounds that
    # are physically unreachable for any PK / beta / K parameter while
    # keeping np.exp() safely below the overflow threshold (~710).
    _LOG_PROP_MAX = 30.0

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
        # FIX: clamp the log-proposal before exp() to prevent overflow.
        _k_log_prop = np.clip(np.log(k) + rng.normal(0.0, k_step),
                              -_LOG_PROP_MAX, _LOG_PROP_MAX)
        k_prop  = np.exp(_k_log_prop)
        prop_lp = _log_posterior(r, k_prop, beta, observed_schedule,
                                 observed_burdens, pk, observed_end_ploidy)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            k, cur_lp = k_prop, prop_lp
            if not in_burnin:
                accept["k"] += 1
        k_hist.append(np.log(k))

        # ── beta (log scale) ──────────────────────────────────────────────────
        # FIX: clamp the log-proposal before exp() to prevent overflow.
        _beta_log_prop = np.clip(np.log(beta) + rng.normal(0.0, beta_step),
                                 -_LOG_PROP_MAX, _LOG_PROP_MAX)
        beta_prop = np.exp(_beta_log_prop)
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
                step = pk_steps[drug][param]
                # FIX: clamp the log-proposal before exp() to prevent overflow.
                # The previous code had a typo: assigned to `val_prop_log_prop`
                # but then read from `_log_prop` (NameError / wrong variable).
                _log_prop = np.clip(
                    np.log(params[param]) + rng.normal(0.0, step),
                    -_LOG_PROP_MAX, _LOG_PROP_MAX,
                )
                val_prop = np.exp(_log_prop)
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
                                  pk_overrides, beta, cycle_len=None):
    new_status, seg_traj = simulate_next_state(ploidy, drug, r_base, k_cap,
                                               pk_overrides=pk_overrides,
                                               cycle_len=cycle_len,
                                               beta=beta)
    return new_status, seg_traj, path, traj, drug, cycle_len


def _beam_search_step(current_beams, executor, r_base, k_cap, beam_width,
                      pk_state=None, beta=BETA_INIT, start_day=0.0):
    futures = []
    # Carry extinct beams forward unchanged so they are never silently dropped.
    # Without this, a beam that first reaches MIN_SIZE at depth D is removed
    # from next_candidates at depth D+1 (the extinct guard skips expansion and
    # never re-adds it), making it impossible for that path to win.
    next_candidates = [
        (burden, ploidy, path, traj, True)
        for burden, ploidy, path, traj, extinct in current_beams
        if extinct
    ]
    for burden, ploidy, path, traj, extinct in current_beams:
        if extinct:
            continue
        # Only allow active treatments after FIRST_TX_DAY; before that,
        # restrict to "none" to mirror the pre-treatment period.
        elapsed = start_day + sum(d_len for _, _, d_len in path)
        available_drugs = DRUGS if elapsed >= FIRST_TX_DAY else ["none"]
        for drug in available_drugs:
            # Beam search explores future treatments at full reference dose
            overrides = _pk_overrides_for(drug, pk_state,
                                          dose_mg_kg=DOSE_REFERENCE_MG_KG)
            # Snap the last pre-treatment "none" cycle to end exactly at
            # FIRST_TX_DAY so treatment begins on the correct calendar day.
            if drug == "none" and elapsed < FIRST_TX_DAY:
                cycle_len = FIRST_TX_DAY - elapsed
            else:
                cycle_len = None  # use drug default
            futures.append(executor.submit(
                _simulate_next_state_wrapper,
                ploidy, drug, path, traj, r_base, k_cap, overrides, beta,
                cycle_len))
    for future in as_completed(futures):
        next_ploidy, seg_traj, old_path, old_traj, drug, actual_cycle_len = future.result()
        new_burden   = sum(next_ploidy.values())
        resolved_len = actual_cycle_len if actual_cycle_len is not None else get_cycle_length(drug)
        segment_info = (drug, len(seg_traj), resolved_len)
        if new_burden < MIN_SIZE:
            next_candidates.append(
                (new_burden, next_ploidy,
                 old_path + [segment_info], old_traj + list(seg_traj), True))
        elif new_burden <= MAX_SIZE:
            next_candidates.append(
                (new_burden, next_ploidy,
                 old_path + [segment_info], old_traj + list(seg_traj), False))
    # Sort key: extinct beams always rank above alive beams.
    # Among extinct beams, prefer fewest days to extinction, breaking ties by burden.
    # Among alive beams, prefer lower burden.
    def _sort_key(x):
        burden, _ploidy, path, _traj, extinct = x
        if extinct:
            days_to_extinction = sum(d_len for _, _, d_len in path)
            return (0, days_to_extinction, burden)   # extinct: fewest days, then lowest burden
        else:
            return (1, 0.0, burden)                  # alive: always below all extinct beams

    next_candidates.sort(key=_sort_key)
    return next_candidates[:beam_width]


def run_single_beam_search(run_idx, r_base, k_cap, beam_width, max_depth,
                            start_ploidy=None, pk_state=None, beta=BETA_INIT,
                            start_day=0.0):
    if start_ploidy is None:
        start_ploidy = dict(INITIAL_PLOIDY)
    initial_burden = sum(start_ploidy.values())
    beam = [(initial_burden, start_ploidy, [],
             [np.array(list(start_ploidy.values()))], False)]
    with ThreadPoolExecutor(max_workers=len(DRUGS) * beam_width) as executor:
        for _ in range(max_depth):
            beam = _beam_search_step(beam, executor, r_base, k_cap,
                                     beam_width, pk_state, beta, start_day)
            if not beam or all(b[4] for b in beam):
                break
    return beam[0] if beam else None


def _beam_search_worker(i, r_i, k_i, beam_width, max_depth,
                         use_observed_end, pk_state, beta_i):
    # The Gibbs sampler only fits PK for drugs present in the observed schedule
    # (typically just gemcitabine).  Fill in init values for every other drug so
    # that the beam search explores all candidate treatments under consistent PK
    # assumptions — exactly as the baseline beam search does (cf. lines 826-828).
    full_pk_state = {drug: dict(params) for drug, params in pk_state.items()}
    for drug, params in PK_PARAMS_TO_FIT.items():
        if drug not in full_pk_state:
            full_pk_state[drug] = {p: cfg["init"] for p, cfg in params.items()}

    if use_observed_end and OBSERVED_DRUGS_ADMINISTERED:
        sp = get_observed_end_ploidy(r_i, k_i, full_pk_state, beta=beta_i)
        if sp is None:
            sp = dict(INITIAL_PLOIDY)
        # Beam starts at the end of the observed schedule — already past FIRST_TX_DAY.
        start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
    else:
        sp = dict(INITIAL_PLOIDY)
        start_day = 0.0
    return run_single_beam_search(i, r_i, k_i, beam_width, max_depth,
                                  start_ploidy=sp, pk_state=full_pk_state,
                                  beta=beta_i, start_day=start_day)


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

    pk_state_map = map_sample["pk"]  # only has gemcitabine
    # Fill in known PK for all other drugs from their init values
    for drug, params in PK_PARAMS_TO_FIT.items():
        if drug not in pk_state_map:
            pk_state_map[drug] = {p: cfg["init"] for p, cfg in params.items()}

    print(f"\nMAP  R={r_base_map:.4f}  K={k_cap_map:.3e}  beta={beta_map:.4g}")



    if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
        map_start_ploidy = get_observed_end_ploidy(r_base_map, k_cap_map,
                                                    pk_state_map, beta=beta_map)
        if map_start_ploidy is None:
            map_start_ploidy = dict(INITIAL_PLOIDY)
        # Beam starts at the end of the observed schedule — already past FIRST_TX_DAY.
        baseline_start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
    else:
        map_start_ploidy = dict(INITIAL_PLOIDY)
        baseline_start_day = 0.0

    baseline_res  = run_single_beam_search(
        "baseline", r_base_map, k_cap_map,
        BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
        start_ploidy=map_start_ploidy,
        pk_state=pk_state_map,
        beta=beta_map,
        start_day=baseline_start_day)
    baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
    print(f"Baseline path: {baseline_path}")

    rng = np.random.default_rng()
    # selected_idx     = rng.choice(len(posterior_samples),
    #                               size=N_SAMPLED_RUNS, p=bma_weights, replace=True)
    selected_idx     = rng.choice(len(posterior_samples),
                                  size=N_SAMPLED_RUNS, replace=True)
    selected_samples = [posterior_samples[i] for i in selected_idx]
    selected_weights = bma_weights[selected_idx]
    selected_weights = selected_weights / selected_weights.sum()

    sampled_results, run_weights, sampled_params = [], [], []
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
                sampled_params.append(selected_samples[i])
            else:
                n_full_maxout += 1

    run_weights = np.array(run_weights)
    run_weights /= run_weights.sum()

    print(f"\n  Full maxouts (all candidates exceeded MAX_SIZE={MAX_SIZE:.0e}): "
          f"{n_full_maxout}/{N_SAMPLED_RUNS} runs "
          f"({100 * n_full_maxout / N_SAMPLED_RUNS:.1f}%) — "
          f"excluded from flip rate, their BMA weight redistributed.")

    print("\n" + "=" * 78)
    print(f"{'Cycle':<7} | {'Baseline Drug':<16} | {'Days':<16} | {'Disagreement Rate':>11}")
    print("-" * 78)
    _cycle_day       = baseline_start_day
    flip_rate_rows   = []
    for i in range(len(baseline_path)):
        target_drug     = baseline_path[i]
        cycle_len       = baseline_res[2][i][2]
        day_start       = int(round(_cycle_day))
        day_end         = int(round(_cycle_day + cycle_len))
        days_str        = f"Day {day_start}–{day_end}"
        _cycle_day     += cycle_len
        unweighted_flip = 0
        active_count    = 0
        for res, w in zip(sampled_results, run_weights):
            sampled_path = [step[0] for step in res[2]]
            if i < len(sampled_path):
                active_count += 1
                if sampled_path[i] != target_drug:
                    unweighted_flip += 1
        raw_rate = (unweighted_flip / active_count) if active_count > 0 else 0.0
        print(f"{i + 1:<7} | {target_drug:<16} | {days_str:<16} | {raw_rate * 100:>9.2f}%")
        flip_rate_rows.append({
            "cycle":             i + 1,
            "baseline_drug":     target_drug,
            "day_start":         day_start,
            "day_end":           day_end,
            "active_runs":       active_count,
            "disagreement_rate": raw_rate,
        })

    print("=" * 78)

    _flip_csv = "flip_rate_table.csv"
    with open(_flip_csv, "w", newline="") as _fh:
        _writer = csv.DictWriter(
            _fh,
            fieldnames=["cycle", "baseline_drug", "day_start", "day_end",
                        "active_runs", "disagreement_rate"],
        )
        _writer.writeheader()
        _writer.writerows(flip_rate_rows)
    print(f"  Saved: {_flip_csv}")

    # =========================================================================
    # Table 2 — Path of each sampled run
    # =========================================================================

    if sampled_results:
        max_cycles = max(len(res[2]) for res in sampled_results)
        drug_col_w = 14   # width of each per-cycle drug column

        # Fixed-width columns appended after the drug sequence
        param_col_w   = 10
        outcome_col_w = 16
        burden_col_w  = 12

        cycle_headers = " | ".join(
            f"{'Cycle ' + str(c + 1):<{drug_col_w}}" for c in range(max_cycles)
        )
        param_header  = (
            f"{'r_base':<{param_col_w}} | "
            f"{'k_cap':<{param_col_w}} | "
            f"{'beta':<{param_col_w}} | "
            f"{'End Burden':<{burden_col_w}} | "
            f"{'Outcome':<{outcome_col_w}}"
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
            final_burden, _ploidy, path, _traj, extinct = res
            sampled_path = [step[0] for step in path]

            # Drug-sequence cells (console)
            cells = [f"{drug:<{drug_col_w}}" for drug in sampled_path]
            while len(cells) < max_cycles:
                cells.append(f"{'—':<{drug_col_w}}")

            # Parameter values
            r_str = f"{samp['r_base']:.4f}"
            k_str = f"{samp['k_cap']:.3e}"
            b_str = f"{samp['beta']:.4g}"

            # End-of-path outcome
            outcome    = "EXTINCT" if extinct else "alive"
            burden_str = f"{final_burden:.3e}"

            param_cells = (
                f"{r_str:<{param_col_w}} | "
                f"{k_str:<{param_col_w}} | "
                f"{b_str:<{param_col_w}} | "
                f"{burden_str:<{burden_col_w}} | "
                f"{outcome:<{outcome_col_w}}"
            )
            print(f"{'Run ' + str(run_idx + 1):<6} | " + " | ".join(cells) + " | " + param_cells)

            # CSV row — one column per cycle, then fixed params/outcome columns
            row: dict = {"run": run_idx + 1}
            for c in range(max_cycles):
                row[f"cycle_{c + 1}"] = sampled_path[c] if c < len(sampled_path) else ""
            row["r_base"]      = samp["r_base"]
            row["k_cap"]       = samp["k_cap"]
            row["beta"]        = samp["beta"]
            row["end_burden"]  = final_burden
            row["outcome"]     = outcome
            path_rows.append(row)

        print("=" * len(header))

        _path_csv    = "sampled_run_paths.csv"
        _cycle_cols  = [f"cycle_{c + 1}" for c in range(max_cycles)]
        _fieldnames  = ["run"] + _cycle_cols + ["r_base", "k_cap", "beta",
                                                 "end_burden", "outcome"]
        with open(_path_csv, "w", newline="") as _fh:
            _writer = csv.DictWriter(_fh, fieldnames=_fieldnames)
            _writer.writeheader()
            _writer.writerows(path_rows)
        print(f"  Saved: {_path_csv}")

    print(f"Total time: {time() - start_time:.2f}s")

    print("\nRunning baseline path simulation for plotting...")

    # Build a schedule from the baseline path decisions.
    baseline_burden_timeline: list[tuple[float, float]] = []
    baseline_ploidy_snapshots: list[tuple[float, dict]] = []

    ploidy_state = dict(map_start_ploidy)
    day = 0.0

    # Record initial state
    baseline_burden_timeline.append((day, float(sum(ploidy_state.values()))))
    baseline_ploidy_snapshots.append((day, dict(ploidy_state)))

    if baseline_res is not None:
        for cycle_idx, (drug, n_seg_steps, cycle_len) in enumerate(baseline_res[2]):
            overrides = _pk_overrides_for(drug, pk_state_map,
                                          dose_mg_kg=DOSE_REFERENCE_MG_KG)
            try:
                new_ploidy, seg_traj = simulate_next_state(
                    ploidy_state, drug, r_base_map, k_cap_map,
                    pk_overrides=overrides, cycle_len=cycle_len, beta=beta_map)
            except Exception as exc:
                print(f"  WARNING: cycle {cycle_idx} ({drug}) failed: {exc}")
                break

            C_fn  = get_concentration_curve(drug, **overrides)
            TIMES = build_times_with_doses(
                (0.0, cycle_len), ODE_STEP,
                drug_name=drug, include_days=True, eps=1e-8,
            )
            _t, Ns_full, T_mat_full, _T_total, _M = simulate_karyotype_ode_piecewise(
                ploidy_state, drug,
                t_span=(0.0, cycle_len),
                r=r_base_map, Kcap=k_cap_map, beta=beta_map,
                N_min=10, N_max=90,
                C_fn=C_fn, f_param_fn=f, t_eval=TIMES, max_step=ODE_STEP,
                renormalize_M=True,
            )
            # T_mat_full: shape (len(Ns_full), len(TIMES))
            # Skip the first column (t=0, already recorded as initial state).
            n_tp = T_mat_full.shape[1]
            for t_idx in range(1, n_tp):
                t = day + TIMES[t_idx]
                col = T_mat_full[:, t_idx]
                total = float(col.sum())
                baseline_burden_timeline.append((t, total))
                ploidy_dict_t = {int(n): float(c) for n, c in zip(Ns_full, col) if c > 0}
                baseline_ploidy_snapshots.append((t, ploidy_dict_t))

            day += cycle_len
            ploidy_state = new_ploidy

    # =========================================================================
    # Pre-compute ploidy bins (shared by Plot 1b and Plot 2)
    #
    # For each cycle snapshot:
    #   • divide each chromosome count (key) by HAPLOID_N (23) → ploidy ratio
    #   • round to the nearest 0.1
    #   • accumulate raw cell counts and fractions per bin
    # =========================================================================

    all_bins: set[float] = set()
    # Each entry: (day, {bin: raw_cell_count}, {bin: fraction})
    ploidy_bin_data: list[tuple[float, dict[float, float], dict[float, float]]] = []

    for snap_day, ploidy_dict in baseline_ploidy_snapshots:
        total_cells = sum(ploidy_dict.values())
        if total_cells <= 0:
            ploidy_bin_data.append((snap_day, {}, {}))
            continue
        binned_counts: dict[float, float] = {}
        for chr_count, cell_count in ploidy_dict.items():
            ratio       = chr_count / HAPLOID_N                      # divide by 23
            rounded_bin = round(round(ratio * 10) / 10, 1)           # nearest 0.1
            binned_counts[rounded_bin] = binned_counts.get(rounded_bin, 0.0) + cell_count
        frac_binned = {b: v / total_cells for b, v in binned_counts.items()}
        ploidy_bin_data.append((snap_day, binned_counts, frac_binned))
        all_bins.update(binned_counts.keys())

    sorted_bins   = sorted(all_bins)
    snap_days_arr = np.array([d[0] for d in ploidy_bin_data])
    cmap          = plt.get_cmap("tab20")

    # Build count and fraction matrices  (rows = bins, cols = time-points)
    count_matrix = np.zeros((len(sorted_bins), len(ploidy_bin_data)))
    frac_matrix  = np.zeros((len(sorted_bins), len(ploidy_bin_data)))
    for col_idx, (_, counts, fracs) in enumerate(ploidy_bin_data):
        for row_idx, b in enumerate(sorted_bins):
            count_matrix[row_idx, col_idx] = counts.get(b, 0.0)
            frac_matrix[row_idx,  col_idx] = fracs.get(b,  0.0)

    # =========================================================================
    # Drug shading helper — shared by all three plots
    # =========================================================================

    drug_colors: dict[str, str] = {
        "gemcitabine": "orange",
        "bay1895344":  "red",
        "alisertib":   "green",
        "ispinesib":   "blue",
        "none":        "yellow",
    }

    def _add_drug_shading(ax, path_info: list[tuple[str, int]]) -> None:
        """Shade the x-axis of *ax* by drug cycle using axvspan."""
        if not path_info:
            return
        current_time = 0.0
        shaded_labels: set[str] = set()
        for drug_name, _seg_len, duration in path_info:
            end_time = current_time + duration
            ax.axvspan(
                current_time, end_time,
                color=drug_colors.get(drug_name, "gray"),
                linewidth=0,
                alpha=0.15,
                label=drug_name if drug_name not in shaded_labels else None,
            )
            shaded_labels.add(drug_name)
            current_time = end_time

    # Collect the baseline path drug sequence (may be None if search failed)
    _baseline_path_info: list[tuple[str, int]] = (
        baseline_res[2] if baseline_res is not None else []
    )

    # =========================================================================
    # Plot 1 — Total tumor burden over the baseline path
    # =========================================================================

    burden_days   = np.array([t for t, _ in baseline_burden_timeline])
    burden_values = np.array([b for _, b in baseline_burden_timeline])

    fig1, ax1 = plt.subplots(figsize=(10, 5))

    # Drug-cycle shading (drawn first so data lines sit on top)
    _add_drug_shading(ax1, _baseline_path_info)

    ax1.plot(burden_days, burden_values, color="steelblue",
             linewidth=1.8, label="Predicted (baseline path)")

    if PLOT_OBSERVED_DATA:
        obs_days = sorted(OBSERVED_TUMOR_BURDENS.keys())
        obs_vals = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days]
        ax1.scatter(obs_days, obs_vals, color="firebrick", zorder=5,
                    label="Observed", s=50)

    ax1.set_xlabel("Day")
    ax1.set_ylabel("Number of cells")
    ax1.set_title(f"Tumor Burden — Baseline Path\n({SAMPLE_NAME})")
    ax1.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1, 1))
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig("baseline_tumor_burden.png", dpi=150)
    print("  Saved: baseline_tumor_burden.png")

    fig1b, ax_combined = plt.subplots(figsize=(12, 6))

    # Drug-cycle shading
    _add_drug_shading(ax_combined, _baseline_path_info)

    # Total burden — cells
    ax_combined.plot(burden_days, burden_values,
                     color="steelblue", linewidth=2.5, zorder=3,
                     label="Total burden (predicted)")

    if PLOT_OBSERVED_DATA:
        obs_days = sorted(OBSERVED_TUMOR_BURDENS.keys())
        obs_vals = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days]
        ax_combined.scatter(obs_days, obs_vals, color="firebrick", zorder=5,
                            label="Observed", s=50)

    # Per-ploidy cell counts — same unit, same axis
    for row_idx, b in enumerate(sorted_bins):
        counts = count_matrix[row_idx, :]
        if counts.max() > 1e-6:
            ax_combined.plot(snap_days_arr, counts,
                             color=cmap(row_idx % 20),
                             linewidth=2.5,
                             label=f"Ploidy {b:.1f}×")

    ax_combined.set_xlabel("Day")
    ax_combined.set_ylabel("Number of cells")
    ax_combined.legend(title="Treatment / Ploidy", loc="upper left", fontsize=7, ncol=2)
    ax_combined.grid(True, alpha=0.3)

    fig1b.suptitle(
        f"Tumor Burden & Ploidy Cell Counts — Baseline Path\n"
        f"({SAMPLE_NAME}  |  ploidy = chromosomes / {HAPLOID_N}, rounded to 0.1)",
        fontsize=10)
    fig1b.tight_layout()
    fig1b.savefig("baseline_burden_and_ploidy_counts.png", dpi=150)
    print("  Saved: baseline_burden_and_ploidy_counts.png")

    # =========================================================================
    # Plot 2 — Ploidy fraction distribution over the baseline path
    # =========================================================================

    fig2, ax2 = plt.subplots(figsize=(11, 6))

    # Drug-cycle shading
    _add_drug_shading(ax2, _baseline_path_info)

    for row_idx, b in enumerate(sorted_bins):
        fracs = frac_matrix[row_idx, :]
        if fracs.max() > 1e-6:
            ax2.plot(snap_days_arr, fracs,
                     color=cmap(row_idx % 20),
                     linewidth=2.5,
                     label=f"Ploidy {b:.1f}×")

    ax2.set_xlabel("Day (end of each cycle)")
    ax2.set_ylabel("Fraction of cells")
    ax2.set_title(f"Ploidy Distribution Over Time — Baseline Path\n"
                  f"(predicted ploidy / {HAPLOID_N}, rounded to nearest 0.1)")
    ax2.legend(title="Treatment / Ploidy", loc="upper right", fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig("baseline_ploidy_distribution.png", dpi=150)
    print("  Saved: baseline_ploidy_distribution.png")

    # =========================================================================
    # Plot 3 — Average ploidy over the baseline path
    #
    # For each time-point snapshot the weighted mean ploidy is:
    #   avg_ploidy(t) = [Σ (chr_count × cell_count)] / total_cells / HAPLOID_N
    # =========================================================================

    avg_ploidy_values: list[float] = []
    for _snap_day, ploidy_dict in baseline_ploidy_snapshots:
        total_cells = sum(ploidy_dict.values())
        if total_cells <= 0:
            avg_ploidy_values.append(float("nan"))
        else:
            weighted_sum = sum(chr_n * cnt for chr_n, cnt in ploidy_dict.items())
            avg_ploidy_values.append(weighted_sum / total_cells / HAPLOID_N)

    avg_ploidy_arr = np.array(avg_ploidy_values)

    fig3, ax3 = plt.subplots(figsize=(10, 5))

    # Drug-cycle shading (drawn first so the line sits on top)
    _add_drug_shading(ax3, _baseline_path_info)

    ax3.plot(snap_days_arr, avg_ploidy_arr,
             color="darkorchid", linewidth=2.0, label="Mean ploidy")

    ax3.set_xlabel("Day")
    ax3.set_ylabel(f"Average ploidy (chromosomes / {HAPLOID_N})")
    ax3.set_title(
        f"Average Ploidy Over Treatment — Baseline Path\n({SAMPLE_NAME})"
    )
    ax3.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1, 1))
    ax3.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig("baseline_average_ploidy.png", dpi=150)
    print("  Saved: baseline_average_ploidy.png")

    plt.show()
