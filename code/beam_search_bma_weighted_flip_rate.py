from __future__ import annotations

import sys

import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import os
from time import time
from PujanEarlyVersionModel import ploidy_forcast

# =============================================================================
# DRUG / MODEL CONFIGURATION
# =============================================================================

# --- ACTIVE (PK + PD fully defined) -----------------------------------------
DRUGS = [
    "gemcitabine",    # IV  — C_peak=0.032, half_life=0.05,  period=7
    "bay1895344",     # IV  — C_peak=0.5,   half_life=0.5,   period=0.5
    "alisertib",      # IV  — C_peak=1.53,  half_life=19.0,  period=7
    "ispinesib",      # IV  — C_peak=0.09,  half_life=1.04,  period=7
    "none",           # IV  — C_peak=0  (no active drug)
    # --- PK defined. Uncomment after adding the corresponding f() entry ------
    # "volasertib",   # IV  — C_peak=1.0,   half_life=4.0,   period=7
    # "cytarabine",   # IV  — C_peak=1.0,   half_life=0.2,   period=3.5
    # "umi-77",       # IV  — C_peak=1.0,   half_life=0.8,   period=7
    # "navitoclax",   # IV  — C_peak=1.0,   half_life=0.73,  period=1
    # "abt-199",      # Oral — dose=100, F=0.6, Vd=250, ka=1.0, ke=0.3/day
    # "abt-263",      # Oral — dose=100, F=0.6, Vd=120, ka=2.0, ke=0.5/day
    # "capecitabine", # Oral — dose=100, F=0.8, Vd=40,  ka=3.0, ke=0.6/day
    # "ceralasertib", # Oral — dose=80,  F=0.5, Vd=100, ka=2.0, ke=0.5/day
    # "osi-027",      # Oral — dose=50,  F=0.5, Vd=80,  ka=1.8, ke=0.5/day
    # "adavosertib",  # Oral — dose=100, F=0.6, Vd=65,  ka=2.4, ke=0.6/day
    # "tegafur",      # Oral — dose=40,  F=0.5, Vd=45,  ka=1.6, ke=0.5/day
    # "tas",          # Oral — dose=60,  F=0.5, Vd=40,  ka=2.4, ke=0.7/day
    # "5-azacytidine",# Oral — dose=100, F=0.2, Vd=40,  ka=3.0, ke=2.0/day
]

MIN_SIZE    = 1e5
MAX_SIZE    = 2e10
DEFAULT_LEN = 7.0

# Treatment-cycle duration (days) used by the BEAM SEARCH (future cycles only).
# For observed cycles the duration is always (end_day - start_day) from the
# schedule below, so CYCLE_LENGTHS is NOT consulted for observed data.
CYCLE_LENGTHS = {
    # IV
    "gemcitabine":    28.0,   # 4-week cycle (weekly infusions x3 + 1 rest week)
    "bay1895344":      7.0,   # twice-daily dosing, 1-week evaluation window
    "alisertib":      21.0,   # 3-week cycle
    "ispinesib":       7.0,
    "volasertib":      7.0,
    "cytarabine":      7.0,
    "umi-77":          7.0,
    "navitoclax":      7.0,   # daily dosing, weekly evaluation window
    "none":            7.0,
    # Oral
    "abt-199":        28.0,
    "abt-263":        28.0,
    "capecitabine":   21.0,   # 14 days on / 7 days off
    "ceralasertib":   28.0,
    "osi-027":        28.0,
    "adavosertib":    28.0,
    "tegafur":        28.0,
    "tas":            28.0,
    "5-azacytidine":  28.0,
}

# Beam-search hyperparameters
BASE_BEAM_WIDTH    = 40
BASE_MAX_DEPTH     = 100
SAMPLED_BEAM_WIDTH = 40
SAMPLED_MAX_DEPTH  = 100
N_SAMPLED_RUNS     = 20

# =============================================================================
# PATIENT DATA
# =============================================================================
# Enter every drug cycle that has already been administered, in order.
#
# Format: list of (start_day, end_day, drug_name) tuples.
#   - start_day / end_day : days since treatment start (day 0).
#   - The model uses (end_day - start_day) as the exact cycle duration.
#     CYCLE_LENGTHS is NOT used for observed cycles.
#   - Any gap between consecutive entries is automatically filled with "none".
#   - Entries must be sorted by start_day (ascending) and must not overlap.

OBSERVED_DRUGS_ADMINISTERED: list[tuple[float, float, str]] = [
    ( 0,  7, "none"),
    ( 7, 14, "none"),
    (14, 21, "none"),
    (21, 28, "none"),
    (28, 31, "none")
]

# Observed tumor-burden measurements: { day_since_tx_start : burden_in_cm3 }
# Day 0 is the first day of treatment.
# Values are in cm3; multiplied by 1e7 cells/cm3 to get cell counts.
_CELLS_PER_CM3 = 1e7

_OBSERVED_TUMOR_BURDENS_CM3 = {
     0:  40.00,
     3: 125.44,
     7: 274.44,
    10: 689.70,
    14: 778.53,
    17: 1056.30,
    21: 1245.46,
    24: 1916.60,
    28: 1767.87,
    31: 1729.65,
}

# Converted to cells automatically - do not edit this line.
OBSERVED_TUMOR_BURDENS = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}
# sys.exit()

# Beam-search starting point
# True  -> beam search begins from the simulated ploidy state at the END of
#          OBSERVED_DRUGS_ADMINISTERED (i.e., forecasting future treatment).
# False -> beam search begins from INITIAL_PLOIDY (full simulation from scratch).
START_BEAM_FROM_OBSERVED_END = False

# Initial ploidy distribution (cells per ploidy class at t = 0)
INITIAL_PLOIDY = {2.0: 3.8e8, 3.0: 0.1e8, 4.0: 0.1e8}

if _OBSERVED_TUMOR_BURDENS_CM3:
    # Make initial ploidy consistent with the first observed burden if not already set.
    initial = _OBSERVED_TUMOR_BURDENS_CM3[0] * _CELLS_PER_CM3
    if sum(INITIAL_PLOIDY.values()) != initial:
        print(f"Warning: Initial ploidy burden ({sum(INITIAL_PLOIDY.values()):.2e} cells) "
              f"does not match first observed burden ({initial:.2e} cells). "
              f"Please adjust INITIAL_PLOIDY or the first entry in _OBSERVED_TUMOR_BURDENS_CM3.")

# =============================================================================
# PRIOR AND GIBBS MCMC SETTINGS
# =============================================================================

# R_BASE prior: Normal(mean, std)
R_BASE_PRIOR_MEAN = 0.575
R_BASE_PRIOR_STD  = 0.10

# K_CAP prior: Log-Normal(log_mean, log_std)
K_CAP_PRIOR_LOG_MEAN = np.log(6e10)
K_CAP_PRIOR_LOG_STD  = 0.8

# Metropolis-within-Gibbs proposal widths for R_BASE and K_CAP
GIBBS_R_STEP     = 0.02   # std of Normal proposal for R_BASE
GIBBS_K_LOG_STEP = 0.30   # std of Normal proposal for log(K_CAP)

# Observation-noise std on log scale (tune to data variability)
LIKELIHOOD_SIGMA = 1.5

# Number of posterior samples / burn-in
N_GIBBS_SAMPLES = 200
GIBBS_BURNIN    = 100

# Cheaper N_SIMS during MCMC (full 1000 used in beam search)
N_SIMS_LIKELIHOOD = 100

# =============================================================================
# PK PARAMETERS TO FIT
# =============================================================================
# For each observed drug, list the PK parameters you want to infer alongside
# R_BASE and K_CAP.  The Gibbs sampler only touches entries for drugs that
# actually appear in OBSERVED_DRUGS_ADMINISTERED — the rest are ignored.
#
# Each parameter entry has:
#   "init"           - starting value (matches PujanEarlyVersionModel default)
#   "prior_log_mean" - mean of the log-Normal prior  (= log of the default)
#   "prior_log_std"  - std  of the log-Normal prior  (1.0 ≈ one order of mag)
#   "step"           - Metropolis proposal std on the log scale
#
# All PK params are strictly positive → sampled on the log scale.
# Comment-out / delete a param entry to hold it fixed at its model default.
#
# IV route params  : C_peak  (peak plasma conc, same units as EC50 in f())
#                    half_life (hours → days in the model)
# Oral route params: dose    (mg, same normalisation as F·dose/Vd in model)
#                    ke_day  (elimination rate constant, 1/day)
# Both routes      : period  (dosing interval in days) — fixed by default

PK_PARAMS_TO_FIT: dict[str, dict[str, dict]] = {

    # ── IV drugs ──────────────────────────────────────────────────────────────
    "gemcitabine": {
        "C_peak": {
            "init": 0.032, "prior_log_mean": np.log(0.032),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 0.05, "prior_log_mean": np.log(0.05),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 7.0, "prior_log_mean": np.log(7.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "bay1895344": {
        "C_peak": {
            "init": 0.5, "prior_log_mean": np.log(0.5),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 0.50, "prior_log_mean": np.log(0.50),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 0.5, "prior_log_mean": np.log(0.5),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "alisertib": {
        "C_peak": {
            "init": 1.53, "prior_log_mean": np.log(1.53),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 19.0, "prior_log_mean": np.log(19.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 7.0, "prior_log_mean": np.log(7.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "ispinesib": {
        "C_peak": {
            "init": 0.09, "prior_log_mean": np.log(0.09),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 1.04, "prior_log_mean": np.log(1.04),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 7.0, "prior_log_mean": np.log(7.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "volasertib": {
        "C_peak": {
            "init": 1.0, "prior_log_mean": np.log(1.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 4.0, "prior_log_mean": np.log(4.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 7.0, "prior_log_mean": np.log(7.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "cytarabine": {
        "C_peak": {
            "init": 1.0, "prior_log_mean": np.log(1.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 0.2, "prior_log_mean": np.log(0.2),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 3.5, "prior_log_mean": np.log(3.5),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "umi-77": {
        "C_peak": {
            "init": 1.0, "prior_log_mean": np.log(1.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 0.8, "prior_log_mean": np.log(0.8),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 7.0, "prior_log_mean": np.log(7.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    "navitoclax": {
        # Listed as IV in drug_dosing_schedules; uses pulsed_dose PK.
        "C_peak": {
            "init": 1.0, "prior_log_mean": np.log(1.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "half_life": {
            "init": 0.73, "prior_log_mean": np.log(0.73),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "period": {"init": 1.0, "prior_log_mean": np.log(1.0),
                   "prior_log_std": 0.3, "step": 0.15},
    },

    # ── Oral drugs ────────────────────────────────────────────────────────────
    # Oral AUC ∝ F·dose / (Vd·ke_day).  F and Vd are strongly correlated with
    # dose and ke_day and hard to identify separately from tumor-burden data,
    # so only dose and ke_day are fitted by default.  Uncomment F / Vd / ka_day
    # entries to add them to the inference if you have supporting PK data.

    "abt-199": {
        "dose": {
            "init": 100.0, "prior_log_mean": np.log(100.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.3, "prior_log_mean": np.log(0.3),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "abt-263": {
        "dose": {
            "init": 100.0, "prior_log_mean": np.log(100.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.5, "prior_log_mean": np.log(0.5),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "capecitabine": {
        "dose": {
            "init": 100.0, "prior_log_mean": np.log(100.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.6, "prior_log_mean": np.log(0.6),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "ceralasertib": {
        "dose": {
            "init": 80.0, "prior_log_mean": np.log(80.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.5, "prior_log_mean": np.log(0.5),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "osi-027": {
        "dose": {
            "init": 50.0, "prior_log_mean": np.log(50.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.5, "prior_log_mean": np.log(0.5),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "adavosertib": {
        "dose": {
            "init": 100.0, "prior_log_mean": np.log(100.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.6, "prior_log_mean": np.log(0.6),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "tegafur": {
        "dose": {
            "init": 40.0, "prior_log_mean": np.log(40.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.5, "prior_log_mean": np.log(0.5),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "tas": {
        "dose": {
            "init": 60.0, "prior_log_mean": np.log(60.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 0.7, "prior_log_mean": np.log(0.7),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    "5-azacytidine": {
        "dose": {
            "init": 100.0, "prior_log_mean": np.log(100.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
        "ke_day": {
            "init": 2.0, "prior_log_mean": np.log(2.0),
            "prior_log_std": 1.0, "step": 0.30,
        },
    },

    # "none" has C_peak = 0 — nothing to fit.
}
# =============================================================================


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_cycle_length(drug: str) -> float:
    """Default cycle length for beam-search (future) cycles only."""
    return CYCLE_LENGTHS.get(drug, DEFAULT_LEN)


def _pk_overrides_for(drug: str, pk_state: dict[str, dict] | None) -> dict:
    """Extract the PK override dict for a single drug from the full pk_state."""
    if pk_state is None:
        return {}
    return dict(pk_state.get(drug, {}))


def _fill_gaps_with_none(
        schedule: list[tuple[float, float, str]],
) -> list[tuple[float, float, str]]:
    """Return schedule with "none" cycles inserted for any uncovered gaps.

    A gap exists when one entry's end_day < the next entry's start_day.
    The inserted "none" entry spans exactly that gap: (prev_end, next_start).
    """
    if not schedule:
        return []
    filled: list[tuple[float, float, str]] = []
    for start, end, drug in schedule:
        if filled and filled[-1][1] < start:
            filled.append((filled[-1][1], start, "none"))
        filled.append((start, end, drug))
    return filled


def _observed_drug_names() -> list[str]:
    """Drug names from the raw schedule (gaps not filled)."""
    return [drug for _, _, drug in OBSERVED_DRUGS_ADMINISTERED]


def _observed_drug_set() -> set[str]:
    """Unique drug names present in OBSERVED_DRUGS_ADMINISTERED."""
    return {drug.lower() for _, _, drug in OBSERVED_DRUGS_ADMINISTERED}


# ---------------------------------------------------------------------------
# Core simulation helpers
# ---------------------------------------------------------------------------

def simulate_next_state(ploidy_status: dict, drug: str,
                        r_base: float, k_cap: float = 6e10,
                        pk_overrides: dict | None = None,
                        cycle_len: float | None = None) -> tuple:
    """Run one treatment cycle; return (new_ploidy_dict, segment_trajectory).

    cycle_len overrides CYCLE_LENGTHS[drug] when provided — always pass
    (end_day - start_day) for observed cycles so exact durations are respected.
    segment_trajectory shape: (n_timepoints, n_ploidies).
    """
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    ploidies, _, _, _, Tpaths = ploidy_forcast(
        ploidy_status, drug, T,
        N_SIMS=1000, R_BASE=r_base, K_CAP=k_cap, **overrides
    )
    final_per_ploidy    = Tpaths[:, :, -1]
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)
    mean_trajectory     = np.mean(Tpaths, axis=0).T   # (n_tp, n_ploidies)
    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k])
                  for k in range(len(ploidies))}
    return new_status, mean_trajectory[1:]


def simulate_next_state_cheap(ploidy_status: dict, drug: str,
                               r_base: float, k_cap: float,
                               pk_overrides: dict | None = None,
                               cycle_len: float | None = None) -> tuple:
    """Same as simulate_next_state but fewer simulations (for MCMC)."""
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    ploidies, _, _, _, Tpaths = ploidy_forcast(
        ploidy_status, drug, T,
        N_SIMS=N_SIMS_LIKELIHOOD, R_BASE=r_base, K_CAP=k_cap, **overrides
    )
    final_per_ploidy    = Tpaths[:, :, -1]
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)
    mean_trajectory     = np.mean(Tpaths, axis=0).T
    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k])
                  for k in range(len(ploidies))}
    return new_status, mean_trajectory[1:]


def get_observed_end_ploidy(r_base: float, k_cap: float,
                             pk_state: dict[str, dict] | None = None) -> dict | None:
    """Simulate through the full observed schedule; return final ploidy state.

    Gaps are auto-filled with "none".  Each cycle's duration is taken from
    the explicit (start_day, end_day) values — CYCLE_LENGTHS is not used.
    Returns None if simulation fails or tumour leaves [MIN_SIZE, MAX_SIZE].
    """
    if not OBSERVED_DRUGS_ADMINISTERED:
        return dict(INITIAL_PLOIDY)

    ploidy   = dict(INITIAL_PLOIDY)
    schedule = _fill_gaps_with_none(OBSERVED_DRUGS_ADMINISTERED)

    for start_day, end_day, drug in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state)
        try:
            ploidy, _ = simulate_next_state(ploidy, drug, r_base, k_cap,
                                            pk_overrides=overrides,
                                            cycle_len=cycle_len)
        except Exception:
            return None
        total = sum(ploidy.values())
        if total < MIN_SIZE or total > MAX_SIZE:
            return None
    return ploidy


# =============================================================================
# GIBBS SAMPLER  (Metropolis-within-Gibbs: R_BASE, K_CAP, PK params)
# =============================================================================

def _simulate_burden_timeline(
        observed_schedule: list[tuple[float, float, str]],
        r_base: float,
        k_cap: float,
        pk_state: dict[str, dict] | None = None,
) -> list | None:
    """Return [(day, total_burden), ...] for the full observed schedule.

    Gaps between entries are auto-filled with "none" cycles.
    Each cycle's duration is (end_day - start_day); timeline points are placed
    at uniformly-spaced internal steps anchored to start_day.
    Uses the cheap SDE so MCMC stays fast.
    """
    ploidy   = dict(INITIAL_PLOIDY)
    timeline = [(0.0, float(sum(ploidy.values())))]
    schedule = _fill_gaps_with_none(observed_schedule)

    for start_day, end_day, drug in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state)
        try:
            new_ploidy, seg_traj = simulate_next_state_cheap(
                ploidy, drug, r_base, k_cap,
                pk_overrides=overrides, cycle_len=cycle_len)
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


def _log_likelihood(r_base: float, k_cap: float,
                    observed_schedule: list[tuple[float, float, str]],
                    observed_burdens: dict,
                    pk_state: dict[str, dict] | None = None) -> float:
    """Log-Normal observation model: -0.5 * [log(obs/pred) / sigma]^2 per point."""
    if not (0.0 < r_base < 1.0) or k_cap <= 0:
        return -np.inf

    timeline = _simulate_burden_timeline(observed_schedule, r_base, k_cap, pk_state)
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

    return log_lik


def _log_prior(r_base: float, k_cap: float,
               pk_state: dict[str, dict] | None = None) -> float:
    """Log-prior: Normal for R_BASE, log-Normal for K_CAP and each PK param."""
    log_p = -0.5 * ((r_base - R_BASE_PRIOR_MEAN) / R_BASE_PRIOR_STD) ** 2

    log_k  = np.log(k_cap)
    log_p += -0.5 * ((log_k - K_CAP_PRIOR_LOG_MEAN) / K_CAP_PRIOR_LOG_STD) ** 2

    if pk_state:
        for drug, params in pk_state.items():
            spec = PK_PARAMS_TO_FIT.get(drug, {})
            for param, val in params.items():
                if val <= 0:
                    return -np.inf
                cfg   = spec[param]
                lv    = np.log(val)
                log_p += -0.5 * ((lv - cfg["prior_log_mean"]) / cfg["prior_log_std"]) ** 2

    return log_p


def _log_posterior(r_base: float, k_cap: float,
                   observed_schedule: list[tuple[float, float, str]],
                   observed_burdens: dict,
                   pk_state: dict[str, dict] | None = None) -> float:
    lp = _log_prior(r_base, k_cap, pk_state)
    if not np.isfinite(lp):
        return -np.inf
    return lp + _log_likelihood(r_base, k_cap, observed_schedule,
                                observed_burdens, pk_state)


def run_gibbs_sampler(observed_schedule: list[tuple[float, float, str]],
                      observed_burdens: dict,
                      n_samples: int = N_GIBBS_SAMPLES,
                      burnin: int    = GIBBS_BURNIN) -> tuple:
    """Metropolis-within-Gibbs for R_BASE, K_CAP, and PK parameters."""
    # Derive active drug set from the raw schedule (gaps not yet filled)
    observed_drug_set = {d.lower() for _, _, d in observed_schedule}

    # Only fit PK params for drugs the patient has actually received.
    active_pk: dict[str, dict] = {
        drug: {p: cfg["init"] for p, cfg in params.items()}
        for drug, params in PK_PARAMS_TO_FIT.items()
        if drug in observed_drug_set
    }
    n_pk = sum(len(v) for v in active_pk.values())
    print(f"--- Gibbs Sampler: R_BASE + K_CAP"
          + (f" + {n_pk} PK param(s)" if n_pk else "") + " ---")
    for drug, params in active_pk.items():
        for p, v in params.items():
            cfg = PK_PARAMS_TO_FIT[drug][p]
            print(f"    {drug}.{p:<12}  init={v:.4g}  "
                  f"prior=logN({cfg['prior_log_mean']:.3f}, {cfg['prior_log_std']:.2f})")

    rng = np.random.default_rng(42)
    r   = R_BASE_PRIOR_MEAN
    k   = np.exp(K_CAP_PRIOR_LOG_MEAN)
    pk  = {drug: dict(params) for drug, params in active_pk.items()}

    cur_lp = _log_posterior(r, k, observed_schedule, observed_burdens, pk)

    samples   = []
    log_posts = []
    accept: dict[str, int] = {"r": 0, "k": 0}
    for drug, params in pk.items():
        for p in params:
            accept[f"{drug}.{p}"] = 0

    total = n_samples + burnin

    for i in range(total):
        if i % 20 == 0:
            phase  = "burn-in" if i < burnin else "sampling"
            pk_str = "  ".join(
                f"{d}.{p}={v:.3g}"
                for d, ps in pk.items() for p, v in ps.items()
            )
            print(f"  iter {i:>4}/{total}  [{phase}]  "
                  f"R={r:.4f}  K={k:.3e}"
                  + (f"  {pk_str}" if pk_str else "")
                  + f"  logP={cur_lp:.2f}")

        # Step 1: update R_BASE
        r_prop  = r + rng.normal(0.0, GIBBS_R_STEP)
        prop_lp = _log_posterior(r_prop, k, observed_schedule, observed_burdens, pk)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            r, cur_lp = r_prop, prop_lp
            if i >= burnin:
                accept["r"] += 1

        # Step 2: update K_CAP on log scale
        lk_prop = np.log(k) + rng.normal(0.0, GIBBS_K_LOG_STEP)
        k_prop  = np.exp(lk_prop)
        prop_lp = _log_posterior(r, k_prop, observed_schedule, observed_burdens, pk)
        if np.log(rng.uniform()) < prop_lp - cur_lp:
            k, cur_lp = k_prop, prop_lp
            if i >= burnin:
                accept["k"] += 1

        # Step 3: update each PK parameter on log scale
        for drug, params in pk.items():
            for param in list(params.keys()):
                step     = PK_PARAMS_TO_FIT[drug][param]["step"]
                cur_val  = params[param]
                lv_prop  = np.log(cur_val) + rng.normal(0.0, step)
                val_prop = np.exp(lv_prop)

                pk_prop               = {d: dict(ps) for d, ps in pk.items()}
                pk_prop[drug][param]  = val_prop

                prop_lp = _log_posterior(r, k, observed_schedule,
                                         observed_burdens, pk_prop)
                if np.log(rng.uniform()) < prop_lp - cur_lp:
                    pk[drug][param] = val_prop
                    cur_lp = prop_lp
                    if i >= burnin:
                        accept[f"{drug}.{param}"] += 1

        if i >= burnin:
            samples.append({
                "r_base": r,
                "k_cap":  k,
                "pk":     {d: dict(ps) for d, ps in pk.items()},
            })
            log_posts.append(cur_lp)

    log_posts = np.array(log_posts)

    print("\n  Acceptance rates:")
    print(f"    R_BASE   : {accept['r'] / n_samples:.2f}")
    print(f"    K_CAP    : {accept['k'] / n_samples:.2f}")
    for drug, params in pk.items():
        for p in params:
            key = f"{drug}.{p}"
            print(f"    {key:<22}: {accept[key] / n_samples:.2f}")

    r_mean = np.mean([s["r_base"] for s in samples])
    k_mean = np.mean([s["k_cap"]  for s in samples])
    print(f"\n  Posterior mean  R_BASE: {r_mean:.4f}  K_CAP: {k_mean:.3e}")
    for drug, params in pk.items():
        for p in params:
            p_mean = np.mean([s["pk"][drug][p] for s in samples])
            print(f"    {drug}.{p} posterior mean: {p_mean:.4g}")

    return samples, log_posts


# =============================================================================
# BMA WEIGHTS
# =============================================================================

def compute_bma_weights(log_posteriors: np.ndarray) -> np.ndarray:
    """Softmax of log-posteriors -> BMA weights (Raftery et al. 2005)."""
    lp = np.asarray(log_posteriors, dtype=float)
    w  = np.exp(lp - np.max(lp))
    return w / w.sum()


# =============================================================================
# BEAM SEARCH
# =============================================================================

def _simulate_next_state_wrapper(ploidy, drug, path, traj, r_base, k_cap,
                                  pk_overrides):
    # Beam search always uses default CYCLE_LENGTHS (no explicit end_day).
    new_status, seg_traj = simulate_next_state(ploidy, drug, r_base, k_cap,
                                               pk_overrides=pk_overrides)
    return new_status, seg_traj, path, traj, drug


def _beam_search_step(current_beams, executor, r_base, k_cap, beam_width,
                      pk_state: dict[str, dict] | None):
    futures = []
    for burden, ploidy, path, traj, extinct in current_beams:
        if extinct:
            continue
        for drug in DRUGS:
            overrides = _pk_overrides_for(drug, pk_state)
            futures.append(executor.submit(
                _simulate_next_state_wrapper,
                ploidy, drug, path, traj, r_base, k_cap, overrides))

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


def run_single_beam_search(run_idx, r_base: float, k_cap: float,
                            beam_width: int, max_depth: int,
                            start_ploidy: dict | None = None,
                            pk_state: dict[str, dict] | None = None):
    if start_ploidy is None:
        start_ploidy = dict(INITIAL_PLOIDY)

    initial_burden = sum(start_ploidy.values())
    beam = [(initial_burden, start_ploidy, [],
             [np.array(list(start_ploidy.values()))], False)]

    with ThreadPoolExecutor(max_workers=len(DRUGS) * beam_width) as executor:
        for _ in range(max_depth):
            beam = _beam_search_step(beam, executor, r_base, k_cap,
                                     beam_width, pk_state)
            if not beam or all(b[4] for b in beam):
                break

    return beam[0] if beam else None


def _beam_search_worker(i: int, r_i: float, k_i: float,
                         beam_width: int, max_depth: int,
                         use_observed_end: bool,
                         pk_state: dict[str, dict] | None) -> tuple | None:
    if use_observed_end and OBSERVED_DRUGS_ADMINISTERED:
        sp = get_observed_end_ploidy(r_i, k_i, pk_state)
        if sp is None:
            sp = dict(INITIAL_PLOIDY)
    else:
        sp = dict(INITIAL_PLOIDY)
    return run_single_beam_search(i, r_i, k_i, beam_width, max_depth,
                                  start_ploidy=sp, pk_state=pk_state)


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    start_time = time()

    # 1. Gibbs sampler: fit R_BASE, K_CAP, and PK params to patient data
    posterior_samples, log_posteriors = run_gibbs_sampler(
        OBSERVED_DRUGS_ADMINISTERED,
        OBSERVED_TUMOR_BURDENS,
        n_samples=N_GIBBS_SAMPLES,
        burnin=GIBBS_BURNIN,
    )

    # BMA weights proportional to posterior probability of each parameter draw
    bma_weights = compute_bma_weights(log_posteriors)

    # MAP sample
    map_idx      = int(np.argmax(log_posteriors))
    map_sample   = posterior_samples[map_idx]
    r_base_map   = map_sample["r_base"]
    k_cap_map    = map_sample["k_cap"]
    pk_state_map = map_sample["pk"]

    print(f"\nMAP estimates  R_BASE: {r_base_map:.4f}  K_CAP: {k_cap_map:.3e}")
    for drug, params in pk_state_map.items():
        for p, v in params.items():
            print(f"  MAP {drug}.{p} = {v:.4g}")

    # 2. Compute beam-search starting ploidy
    if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
        print("\n--- Computing end-of-observation ploidy state (MAP params) ---")
        map_start_ploidy = get_observed_end_ploidy(r_base_map, k_cap_map,
                                                    pk_state_map)
        if map_start_ploidy is None:
            print("  WARNING: simulation failed at MAP params - "
                  "falling back to INITIAL_PLOIDY")
            map_start_ploidy = dict(INITIAL_PLOIDY)
        else:
            total  = sum(map_start_ploidy.values())
            filled = _fill_gaps_with_none(OBSERVED_DRUGS_ADMINISTERED)
            print(f"  Start ploidy total burden: {total:.3e} cells  "
                  f"(after {len(filled)} cycle(s) incl. gap fills)")
    else:
        map_start_ploidy = dict(INITIAL_PLOIDY)
        print("\n  Using INITIAL_PLOIDY as start (START_BEAM_FROM_OBSERVED_END=False)")

    # 3. Baseline beam search at MAP parameters
    print(f"\n--- Baseline beam search (MAP: R={r_base_map:.4f}, K={k_cap_map:.3e}) ---")
    baseline_res  = run_single_beam_search(
        "baseline", r_base_map, k_cap_map,
        BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
        start_ploidy=map_start_ploidy,
        pk_state=pk_state_map)
    baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
    print(f"Baseline path: {baseline_path}")

    # 4. BMA ensemble: draw N_SAMPLED_RUNS from the posterior
    print(f"\n--- BMA ensemble ({N_SAMPLED_RUNS} runs drawn from posterior) ---")
    rng = np.random.default_rng()

    selected_idx     = rng.choice(len(posterior_samples),
                                  size=N_SAMPLED_RUNS, p=bma_weights, replace=True)
    selected_samples = [posterior_samples[i] for i in selected_idx]
    selected_weights = bma_weights[selected_idx]
    selected_weights = selected_weights / selected_weights.sum()

    sampled_results  = []
    run_weights      = []
    use_observed_end = START_BEAM_FROM_OBSERVED_END and bool(OBSERVED_DRUGS_ADMINISTERED)

    with ProcessPoolExecutor(max_workers=min(N_SAMPLED_RUNS, os.cpu_count())) as pool:
        future_map = {
            pool.submit(
                _beam_search_worker,
                i,
                float(selected_samples[i]["r_base"]),
                float(selected_samples[i]["k_cap"]),
                SAMPLED_BEAM_WIDTH,
                SAMPLED_MAX_DEPTH,
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

    run_weights = np.array(run_weights)
    run_weights /= run_weights.sum()

    # 5. BMA-weighted cycle-by-cycle flip rate
    print("\n" + "=" * 65)
    print(f"{'Cycle':<7} | {'Baseline Drug':<16} | "
          f"{'BMA Flip Rate':>14} | {'Unweighted':>11}")
    print("-" * 65)

    for i in range(len(baseline_path)):
        target_drug     = baseline_path[i]
        weighted_flip   = 0.0
        active_weight   = 0.0
        unweighted_flip = 0
        active_count    = 0

        for res, w in zip(sampled_results, run_weights):
            sampled_path = [step[0] for step in res[2]]
            if i < len(sampled_path):
                active_weight += w
                active_count  += 1
                if sampled_path[i] != target_drug:
                    weighted_flip   += w
                    unweighted_flip += 1

        bma_rate = (weighted_flip   / active_weight) if active_weight > 0 else 0.0
        raw_rate = (unweighted_flip / active_count)  if active_count  > 0 else 0.0

        print(f"{i + 1:<7} | {target_drug:<16} | "
              f"{bma_rate * 100:>12.2f}%  | {raw_rate * 100:>9.2f}%")

    print("=" * 65)
    print(f"Total time: {time() - start_time:.2f}s")