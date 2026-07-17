from __future__ import annotations

import csv
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for batch/headless runs
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from concurrent.futures import ProcessPoolExecutor, as_completed
from time import time
from ploidy_model_wgd_missegg_transit import (
    simulate_karyotype_ode_piecewise,
    get_concentration_curve,
    build_times_with_doses,
    f,
)

_f = f
from mcmc_fit import run_mcmc

# PARAMETERS:
MIN_SIZE = 1e5
MAX_SIZE = 2e10
DEFAULT_LEN = 7.0

BASE_BEAM_WIDTH =    100
BASE_MAX_DEPTH =     100
SAMPLED_BEAM_WIDTH = 100
SAMPLED_MAX_DEPTH =  100
N_SAMPLED_RUNS =     100

ODE_STEP = 0.05

# Fallback cells-per-cm^3
_CELLS_PER_CM3 = 1e7

r_base = 0.28
k_cap = 6e10
beta = 0.05
BETA_INIT = beta

HAPLOID_N: int = 23

OXYGEN_HAPLOID_N: int = 22

_CONFIG_PATH = "../config/drug_kinetics.json"


def _load_drug_kinetics_config(path: str = _CONFIG_PATH) -> dict:
    """Load drug kinetics parameters from the JSON config file.
    """
    with open(path) as fh:
        raw = json.load(fh)

    pk_params: dict[str, dict[str, dict]] = {}
    for drug, params in raw["PK_PARAMS_TO_FIT"].items():
        pk_params[drug] = {}
        for param, cfg in params.items():
            pk_params[drug][param] = {
                "init": cfg["init"],
                "prior_log_mean": math.log(cfg["init"]),
                "prior_log_std": cfg["prior_log_std"],
                "step": cfg["step"],
            }

    return {
        "DOSE_REFERENCE_MG_KG": float(raw["DOSE_REFERENCE_MG_KG"]),
        "CYCLE_LENGTHS": {k: float(v) for k, v in raw["CYCLE_LENGTHS"].items()},
        "PK_PARAMS_TO_FIT": pk_params,
    }


_DK = _load_drug_kinetics_config()

EXCEL_PATH = "../data/InVivoData_Gemcitabine/dt_Gem_VT_20260209_v5.xlsx"

DOSE_REFERENCE_MG_KG: float = _DK["DOSE_REFERENCE_MG_KG"]


def load_harvest_data(
        excel_path: str,
        harvest_name: str,
        verbose: bool = True,
) -> tuple[dict[int, float], list[tuple[float, float, str, float]], str]:
    df = pd.read_excel(excel_path)
    matches = df[df["harvest"] == harvest_name]
    if matches.empty:
        raise ValueError(f"No rows found with harvest == '{harvest_name}'")

    row = matches.iloc[0]
    cols = list(df.columns)

    start_idx = cols.index("Day_0")
    day_cols = cols[start_idx:]
    values = row[day_cols].dropna()

    dose_mg_kg = float(row["Dose"])
    first_tx = float(row["Day of 1st treatment"])

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
        print(f"  Ploidy CBS name : {ploidy_name}")

    return burdens_cm3, schedule, ploidy_name


SAMPLE_NAME = "SUM159-2N-120-RL_harvest"

_OBSERVED_TUMOR_BURDENS_CM3, OBSERVED_DRUGS_ADMINISTERED, PLOIDY_SAMPLE_NAME = (
    load_harvest_data(EXCEL_PATH, SAMPLE_NAME, verbose=(__name__ == "__main__"))
)

FIRST_TX_DAY: float = next(
    (start for start, _end, drug, _dose in OBSERVED_DRUGS_ADMINISTERED
     if drug != "none"),
    0.0,
)

DRUGS = [
    "gemcitabine",
    "bay1895344",
    "alisertib",
    "ispinesib",
    "none",
]

# Loaded from config/drug_kinetics.json
CYCLE_LENGTHS: dict[str, float] = _DK["CYCLE_LENGTHS"]

# OBSERVED_TUMOR_BURDENS (in cell counts) is built further down, after
# INITIAL_PLOIDY is loaded — the cm^3 -> cells conversion uses a per-sample
# scale derived from the ploidy file, not the hardcoded _CELLS_PER_CM3.

PLOIDY_CSV_PATH: str = "../data/InVivoData_Gemcitabine/all_ploidy.csv"

OXYGEN_COUNT_TSV_PATH: str = (
    "../data/InVivoData_Gemcitabine/"
    "treatment_day_ploidy_counts_fit_invivo_o2_supply_demand_MAP_seed28.tsv"
)


def load_treatment_day_ploidy_from_oxygen(
        count_tsv_path: str,
        harvest_name: str,
        oxygen_haploid_n: int = OXYGEN_HAPLOID_N,
        label_filter: str | None = "live",
        verbose: bool = True,
) -> dict[int, float]:
    df = pd.read_csv(count_tsv_path, sep="\t")
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
        total = sum(result.values())
        tx_day = sub["treatment_day"].iloc[0]
        print(f"  Oxygen model ploidy : {harvest_name}")
        print(f"  Treatment day       : {tx_day}")
        print(f"  Label filter        : {label_filter}")
        print(f"  Ploidy states       : {len(result)}")
        print(f"  Total cells (init)  : {total:.3e}")

    return result


def load_ploidy_distribution(
        csv_path: str = PLOIDY_CSV_PATH,
        sample_name: str = PLOIDY_SAMPLE_NAME,
        verbose: bool = True,
) -> np.ndarray:
    try:
        df = pd.read_csv(csv_path, sep="\t")
    except FileNotFoundError:
        if verbose:
            print(f"  WARNING: ploidy CSV not found at '{csv_path}'. "
                  "OBSERVED_END_PLOIDY_DISTRIBUTION will be empty — "
                  "biopsy likelihood term disabled.")
        return np.array([], dtype=float)

    mask = df["file"].str.contains(sample_name, na=False)
    values = df.loc[mask, "total_chromosomes"].to_numpy(dtype=float)

    if values.size == 0:
        if verbose:
            print(f"  WARNING: no rows matched sample '{sample_name}' in "
                  f"'{csv_path}'. Biopsy likelihood term disabled.")
    else:
        if verbose:
            print(f"  Loaded {values.size} ploidy values for '{sample_name}' "
                  f"(mean={values.mean():.4f}, std={values.std():.4f})")

    return values


OBSERVED_END_PLOIDY_DISTRIBUTION: np.ndarray = load_ploidy_distribution(
    verbose=(__name__ == "__main__")
)

START_BEAM_FROM_OBSERVED_END = False

PLOT_OBSERVED_DATA: bool = False

# Load treatment-day ploidy composition from the oxygen pre-treatment model
# (live cells only — dead states shouldn't seed the ODE as if proliferating).
# The cell counts in this file are trusted as the correct absolute magnitudes
# at FIRST_TX_DAY; they are NOT rescaled to match the observed burden. Instead,
# they anchor the per-sample _CELLS_PER_CM3 used for the cm^3 -> cells
# conversion of observed tumor volumes (see derivation block below).
INITIAL_PLOIDY: dict[int, float] = load_treatment_day_ploidy_from_oxygen(
    OXYGEN_COUNT_TSV_PATH,
    SAMPLE_NAME,
    verbose=(__name__ == "__main__"),
)

# Derive the per-sample cells-per-cm^3 factor from INITIAL_PLOIDY.
#
# The ploidy TSV gives the correct *absolute* cell counts at FIRST_TX_DAY; the
# Excel sheet gives the tumor *volume* at the same day. The ratio of these two
# is the real cells/cm^3 for this sample — and it's what we should use to
# convert every observed tumor volume (cm^3) into a cell count. This replaces
# the previously-used hardcoded _CELLS_PER_CM3 = 1e7.
#
# INITIAL_PLOIDY itself is left alone: its cell counts are already correct.
if INITIAL_PLOIDY and _OBSERVED_TUMOR_BURDENS_CM3:
    _anchor_day = int(FIRST_TX_DAY)
    if _anchor_day not in _OBSERVED_TUMOR_BURDENS_CM3:
        raise ValueError(
            f"No observed tumor burden at FIRST_TX_DAY={_anchor_day} "
            f"(from 'Day of 1st treatment' column). "
            f"Observed days: {sorted(_OBSERVED_TUMOR_BURDENS_CM3.keys())}."
        )
    _initial_cells = sum(INITIAL_PLOIDY.values())
    _initial_vol_cm3 = _OBSERVED_TUMOR_BURDENS_CM3[_anchor_day]
    if _initial_cells > 0 and _initial_vol_cm3 > 0:
        _CELLS_PER_CM3 = _initial_cells / _initial_vol_cm3
        if __name__ == "__main__":
            print(f"  Derived _CELLS_PER_CM3 = {_CELLS_PER_CM3:.3e} cells/cm^3 "
                  f"(= {_initial_cells:.3e} cells / {_initial_vol_cm3:.3g} cm^3 "
                  f"at day {_anchor_day})")
    elif __name__ == "__main__":
        print(f"  WARNING: could not derive _CELLS_PER_CM3 "
              f"(initial cells={_initial_cells:.3e}, "
              f"initial vol={_initial_vol_cm3:.3g} cm^3) — "
              f"falling back to 1e7.")

# Now that _CELLS_PER_CM3 reflects this sample, convert the observed tumor
# volume timeline (cm^3) into cell counts for the likelihood.
OBSERVED_TUMOR_BURDENS = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}

# Loaded from config/drug_kinetics.json (prior_log_mean derived as log(init))
PK_PARAMS_TO_FIT: dict[str, dict[str, dict]] = _DK["PK_PARAMS_TO_FIT"]


def get_cycle_length(drug: str) -> float:
    return CYCLE_LENGTHS.get(drug, DEFAULT_LEN)


GEMCITABINE_INTRACELLULAR_CPEAK_COEFF = 0.00971
# Maximum clinical gemcitabine dose used during beam-search exploration.
# When this is passed as dose_mg_kg, C_peak is derived from the intracellular
# PK/PD model rather than being taken from pk_state.
GEMCITABINE_BEAM_DOSE_MG_KG: float = 120.0


def _pk_overrides_for(
        drug: str,
        pk_state: dict[str, dict] | None,
        dose_mg_kg: float | None = None,
) -> dict:
    """Build PK override kwargs for *drug*.

    For gemcitabine, C_peak handling depends on whether a clinical dose is
    available:
      - *dose_mg_kg* provided → C_peak is computed from the intracellular
        PK/PD linear model (C_peak = 0.00971 × dose_mg_kg) and overwrites
        any value in pk_state.
      - *dose_mg_kg* is None → C_peak is left as the tuned value from
        pk_state (i.e. treated as a free parameter).
    """
    overrides: dict = dict(pk_state.get(drug, {})) if pk_state else {}
    if drug == "gemcitabine" and dose_mg_kg is not None:
        # Dose is known: derive C_peak from the intracellular PK/PD linear
        # model and override whatever pk_state may hold.
        overrides["C_peak"] = GEMCITABINE_INTRACELLULAR_CPEAK_COEFF * dose_mg_kg
    # If dose_mg_kg is None, C_peak is left as the tuned value from pk_state.
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
             beta: float = BETA_INIT,
             p_wgd: float = 0.0) -> tuple[dict, np.ndarray]:
    """ Simulates out 1 cycle
    """
    C_fn = get_concentration_curve(drug, **pk_overrides)
    TIMES = build_times_with_doses(
        (0.0, cycle_len), dt,
        drug_name=drug, include_days=True, eps=1e-8,
    )
    _t, Ns, x0, x1, _xtot, _B0, _B1, _BW, _Z = simulate_karyotype_ode_piecewise(
        ploidy_status, drug,
        t_span=(0.0, cycle_len),
        r=r_base, Kcap=k_cap, beta=beta,
        N_min=10, N_max=90,
        p_wgd=p_wgd,
        boundary="drop",
        C_fn=C_fn, f_param_fn=f, t_eval=TIMES, max_step=dt,
        renormalize_M=True,
        n_tr=3, k_tr=1.0, k_kill=1.0,
    )
    T_mat = x0 + x1
    final_counts = T_mat[:, -1]
    new_status = {int(N): float(c) for N, c in zip(Ns, final_counts) if c > 0}
    trajectory = T_mat.T[1:]
    return new_status, trajectory


def simulate_next_state(ploidy_status, drug, r_base, k_cap=6e10,
                        pk_overrides=None, cycle_len=None, beta=BETA_INIT,
                        p_wgd: float = 0.0):
    """ ODE wrapper
    """
    overrides = pk_overrides or {}
    T = cycle_len if cycle_len is not None else get_cycle_length(drug)
    return _run_ode(ploidy_status, drug, r_base, k_cap, overrides, T,
                    dt=ODE_STEP, beta=beta, p_wgd=p_wgd)


def get_observed_end_ploidy(r_base, k_cap, pk_state=None, beta=BETA_INIT,
                            p_wgd: float = 0.0):
    if not OBSERVED_DRUGS_ADMINISTERED:
        return dict(INITIAL_PLOIDY)
    ploidy = dict(INITIAL_PLOIDY)
    schedule = _fill_gaps_with_none(OBSERVED_DRUGS_ADMINISTERED)
    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            ploidy, _ = simulate_next_state(ploidy, drug, r_base, k_cap,
                                            pk_overrides=overrides,
                                            cycle_len=cycle_len, beta=beta,
                                            p_wgd=p_wgd)
        except Exception:
            return None
        total = sum(ploidy.values())
        if total < MIN_SIZE or not np.isfinite(total):
            return None
    return ploidy


def _simulate_burden_timeline(observed_schedule, r_base, k_cap,
                              pk_state=None, beta=BETA_INIT, p_wgd: float = 0.0):
    """ Simulates entire schedule proposal
    """
    ploidy = dict(INITIAL_PLOIDY)
    # Simulation starts at the treatment day, not at implantation day 0.
    timeline = [(FIRST_TX_DAY, float(sum(ploidy.values())))]
    schedule = _fill_gaps_with_none(observed_schedule)
    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            new_ploidy, seg_traj = simulate_next_state(
                ploidy, drug, r_base, k_cap,
                pk_overrides=overrides, cycle_len=cycle_len, beta=beta,
                p_wgd=p_wgd)
        except Exception:
            return None
        n_tp = len(seg_traj)
        if n_tp > 0:
            for t_idx, burden_by_ploidy in enumerate(seg_traj):
                day = start_day + (t_idx + 1) * cycle_len / n_tp
                total = float(np.asarray(burden_by_ploidy).sum())
                timeline.append((day, total))
        ploidy = new_ploidy
    return timeline


def _simulate_burden_and_ploidy(observed_schedule, r_base, k_cap,
                                pk_state=None, beta=BETA_INIT,
                                p_wgd: float = 0.0):
    """Simulate the full schedule ONCE and return both burden timeline and
    final ploidy distribution.  This replaces the pattern of calling
    _simulate_burden_timeline + get_observed_end_ploidy separately (which
    ran the same ODE twice).

    Returns (timeline, final_ploidy) or (None, None) on failure.
    """
    ploidy = dict(INITIAL_PLOIDY)
    timeline = [(FIRST_TX_DAY, float(sum(ploidy.values())))]
    schedule = _fill_gaps_with_none(observed_schedule)
    for start_day, end_day, drug, dose in schedule:
        cycle_len = end_day - start_day
        overrides = _pk_overrides_for(drug, pk_state, dose_mg_kg=dose)
        try:
            new_ploidy, seg_traj = simulate_next_state(
                ploidy, drug, r_base, k_cap,
                pk_overrides=overrides, cycle_len=cycle_len, beta=beta,
                p_wgd=p_wgd)
        except Exception:
            return None, None
        n_tp = len(seg_traj)
        if n_tp > 0:
            for t_idx, burden_by_ploidy in enumerate(seg_traj):
                day = start_day + (t_idx + 1) * cycle_len / n_tp
                total = float(np.asarray(burden_by_ploidy).sum())
                timeline.append((day, total))
        total = sum(new_ploidy.values())
        if total < MIN_SIZE or not np.isfinite(total):
            return None, None
        ploidy = new_ploidy
    return timeline, ploidy


def _simulate_next_state_wrapper(ploidy, drug, path, traj, r_base, k_cap,
                                 pk_overrides, beta, cycle_len=None,
                                 p_wgd: float = 0.0):
    new_status, seg_traj = simulate_next_state(ploidy, drug, r_base, k_cap,
                                               pk_overrides=pk_overrides,
                                               cycle_len=cycle_len,
                                               beta=beta, p_wgd=p_wgd)
    return new_status, seg_traj, path, traj, drug, cycle_len


def _beam_search_step(current_beams, r_base, k_cap, beam_width,
                      pk_state=None, beta=BETA_INIT, start_day=0.0,
                      p_wgd: float = 0.0):
    # Carry extinct beams forward unchanged so they are never silently dropped.
    next_candidates = [
        (burden, ploidy, path, traj, True)
        for burden, ploidy, path, traj, extinct in current_beams
        if extinct
    ]
    for burden, ploidy, path, traj, extinct in current_beams:
        if extinct:
            continue
        elapsed = start_day + sum(d_len for _, _, d_len in path)
        available_drugs = DRUGS if elapsed >= FIRST_TX_DAY else ["none"]
        for drug in available_drugs:
            # Beam search explores future treatments at the reference dose,
            # except gemcitabine which is always evaluated at its maximum
            # clinical dose so C_peak is derived accordingly.
            beam_dose = (
                GEMCITABINE_BEAM_DOSE_MG_KG
                if drug == "gemcitabine"
                else DOSE_REFERENCE_MG_KG
            )
            overrides = _pk_overrides_for(drug, pk_state,
                                          dose_mg_kg=beam_dose)
            if drug == "none" and elapsed < FIRST_TX_DAY:
                cycle_len = FIRST_TX_DAY - elapsed
            else:
                cycle_len = None  # use drug default
            next_ploidy, seg_traj, old_path, old_traj, drug_out, actual_cycle_len = (
                _simulate_next_state_wrapper(
                    ploidy, drug, path, traj, r_base, k_cap, overrides, beta,
                    cycle_len, p_wgd=p_wgd))
            new_burden = sum(next_ploidy.values())
            resolved_len = actual_cycle_len if actual_cycle_len is not None else get_cycle_length(drug_out)
            segment_info = (drug_out, len(seg_traj), resolved_len)
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
            return (0, days_to_extinction, burden)  # extinct: fewest days, then lowest burden
        else:
            return (1, 0.0, burden)  # alive: always below all extinct beams

    next_candidates.sort(key=_sort_key)
    return next_candidates[:beam_width]


def run_single_beam_search(run_idx, r_base, k_cap, beam_width, max_depth,
                           start_ploidy=None, pk_state=None, beta=BETA_INIT,
                           start_day=0.0, p_wgd: float = 0.0):
    if start_ploidy is None:
        start_ploidy = dict(INITIAL_PLOIDY)
    initial_burden = sum(start_ploidy.values())
    beam = [(initial_burden, start_ploidy, [],
             [np.array(list(start_ploidy.values()))], False)]
    for _ in range(max_depth):
        beam = _beam_search_step(beam, r_base, k_cap,
                                 beam_width, pk_state, beta, start_day,
                                 p_wgd=p_wgd)
        if not beam or all(b[4] for b in beam):
            break
    return beam[0] if beam else None


def _beam_search_worker(i, r_i, k_i, beam_width, max_depth,
                        use_observed_end, pk_state, beta_i,
                        p_wgd_i: float = 0.0):
    # Fill in init values for drugs not updated by the Metropolis-Within_Gibbs
    full_pk_state = {drug: dict(params) for drug, params in pk_state.items()}
    for drug, params in PK_PARAMS_TO_FIT.items():
        if drug not in full_pk_state:
            full_pk_state[drug] = {p: cfg["init"] for p, cfg in params.items()}

    if use_observed_end and OBSERVED_DRUGS_ADMINISTERED:
        sp = get_observed_end_ploidy(r_i, k_i, full_pk_state, beta=beta_i,
                                     p_wgd=p_wgd_i)
        if sp is None:
            sp = dict(INITIAL_PLOIDY)
        # Beam starts at the end of the observed schedule — already past FIRST_TX_DAY.
        start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
    else:
        sp = dict(INITIAL_PLOIDY)
        # Simulation begins at treatment day; beam search explores from there.
        start_day = FIRST_TX_DAY
    return run_single_beam_search(i, r_i, k_i, beam_width, max_depth,
                                  start_ploidy=sp, pk_state=full_pk_state,
                                  beta=beta_i, start_day=start_day,
                                  p_wgd=p_wgd_i)


if __name__ == "__main__":
    start_time = time()

    RESULTS_DIR = "results"
    os.makedirs(RESULTS_DIR, exist_ok=True)

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
    )  # Parameter fitting

    r_base_map = _map_result["r_base"]
    k_cap_map = _map_result["k_cap"]
    beta_map = _map_result["beta"]
    p_wgd_map = _map_result["p_wgd"]
    pk_state_map = _map_result["pk"]

    print(f"\nMCMC fitting done in {time() - start_time:.1f}s")
    print(f"  MAP: r={r_base_map:.5f}  K={k_cap_map:.3e}  β={beta_map:.5f}  "
          f"p_wgd={p_wgd_map:.2e}")
    print(f"  Posterior samples: {len(posterior_samples)}")

    if START_BEAM_FROM_OBSERVED_END and OBSERVED_DRUGS_ADMINISTERED:
        map_start_ploidy = get_observed_end_ploidy(r_base_map, k_cap_map,
                                                   pk_state_map, beta=beta_map,
                                                   p_wgd=p_wgd_map)
        if map_start_ploidy is None:
            map_start_ploidy = dict(INITIAL_PLOIDY)
        # Beam starts at the end of the observed schedule
        baseline_start_day = OBSERVED_DRUGS_ADMINISTERED[-1][1]
    else:
        map_start_ploidy = dict(INITIAL_PLOIDY)
        baseline_start_day = FIRST_TX_DAY

    baseline_res = run_single_beam_search(
        "baseline", r_base_map, k_cap_map,
        BASE_BEAM_WIDTH, BASE_MAX_DEPTH,
        start_ploidy=map_start_ploidy,
        pk_state=pk_state_map,
        beta=beta_map,
        start_day=baseline_start_day,
        p_wgd=p_wgd_map)
    baseline_path = [step[0] for step in baseline_res[2]] if baseline_res else []
    print(f"Baseline path: {baseline_path}")

    rng = np.random.default_rng()
    selected_idx = rng.choice(len(posterior_samples), size=N_SAMPLED_RUNS, replace=True)
    selected_samples = [posterior_samples[i] for i in selected_idx]

    sampled_results, run_weights, sampled_params = [], [], []
    use_observed_end = START_BEAM_FROM_OBSERVED_END and bool(OBSERVED_DRUGS_ADMINISTERED)

    n_full_maxout = 0

    # Beam search run
    _env_workers = int(os.environ.get("HARVEST_MAX_WORKERS", 0))
    _beam_workers = min(N_SAMPLED_RUNS, _env_workers or os.cpu_count())
    with ProcessPoolExecutor(max_workers=_beam_workers) as pool:
        future_map = {
            pool.submit(
                _beam_search_worker, i,
                float(selected_samples[i]["r_base"]),
                float(selected_samples[i]["k_cap"]),
                SAMPLED_BEAM_WIDTH, SAMPLED_MAX_DEPTH,
                use_observed_end,
                selected_samples[i]["pk"],
                float(selected_samples[i]["beta"]),
                float(selected_samples[i]["p_wgd"]),
            ): i
            for i in range(N_SAMPLED_RUNS)
        }
        for future in as_completed(future_map):
            i = future_map[future]
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
          f"({100 * n_full_maxout / N_SAMPLED_RUNS:.1f}%)")

    print("\n" + "=" * 78)
    print(f"{'Cycle':<7} | {'Baseline Drug':<16} | {'Days':<16} | {'Disagreement Rate':>11}")
    print("-" * 78)
    _cycle_day = baseline_start_day
    flip_rate_rows = []
    for i in range(len(baseline_path)):
        target_drug = baseline_path[i]
        cycle_len = baseline_res[2][i][2]
        day_start = int(round(_cycle_day))
        day_end = int(round(_cycle_day + cycle_len))
        days_str = f"Day {day_start}–{day_end}"
        _cycle_day += cycle_len
        unweighted_flip = 0
        active_count = 0
        for res, w in zip(sampled_results, run_weights):
            sampled_path = [step[0] for step in res[2]]
            if i < len(sampled_path):
                active_count += 1
                end_mismatch = sampled_path[i] != target_drug
                start_mismatch = (i > 0 and sampled_path[i - 1] != baseline_path[i - 1])
                if end_mismatch or start_mismatch:
                    unweighted_flip += 1
        raw_rate = (unweighted_flip / active_count) if active_count > 0 else 0.0
        print(f"{i + 1:<7} | {target_drug:<16} | {days_str:<16} | {raw_rate * 100:>9.2f}%")
        flip_rate_rows.append({
            "cycle": i + 1,
            "baseline_drug": target_drug,
            "day_start": day_start,
            "day_end": day_end,
            "active_runs": active_count,
            "disagreement_rate": raw_rate,
        })

    # Plotting/logging:

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
        max_cycles = max(len(res[2]) for res in sampled_results)
        drug_col_w = 14  # width of each per-cycle drug column

        # Fixed-width columns appended after the drug sequence
        param_col_w = 10
        outcome_col_w = 16
        burden_col_w = 12

        cycle_headers = " | ".join(
            f"{'Cycle ' + str(c + 1):<{drug_col_w}}" for c in range(max_cycles)
        )
        param_header = (
            f"{'r_base':<{param_col_w}} | "
            f"{'k_cap':<{param_col_w}} | "
            f"{'beta':<{param_col_w}} | "
            f"{'End Burden':<{burden_col_w}} | "
            f"{'Outcome':<{outcome_col_w}}"
        )
        header = f"{'Run':<6} | {cycle_headers} | {param_header}"
        sep = "-" * len(header)

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
            outcome = "EXTINCT" if extinct else "alive"
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
            row["r_base"] = samp["r_base"]
            row["k_cap"] = samp["k_cap"]
            row["beta"] = samp["beta"]
            row["end_burden"] = final_burden
            row["outcome"] = outcome
            path_rows.append(row)

        print("=" * len(header))

        _path_csv = f"{RESULTS_DIR}/sampled_run_paths.csv"
        _cycle_cols = [f"cycle_{c + 1}" for c in range(max_cycles)]
        _fieldnames = ["run"] + _cycle_cols + ["r_base", "k_cap", "beta",
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
    day = FIRST_TX_DAY

    # Record initial state
    baseline_burden_timeline.append((day, float(sum(ploidy_state.values()))))
    baseline_ploidy_snapshots.append((day, dict(ploidy_state)))

    if baseline_res is not None:
        for cycle_idx, (drug, n_seg_steps, cycle_len) in enumerate(baseline_res[2]):
            overrides = _pk_overrides_for(drug, pk_state_map,
                                          dose_mg_kg=DOSE_REFERENCE_MG_KG)
            C_fn = get_concentration_curve(drug, **overrides)
            TIMES = build_times_with_doses(
                (0.0, cycle_len), ODE_STEP,
                drug_name=drug, include_days=True, eps=1e-8,
            )
            _t, Ns_full, x0_full, x1_full, _xtot_full, _B0, _B1, _BW, _Z_full = simulate_karyotype_ode_piecewise(
                ploidy_state, drug,
                t_span=(0.0, cycle_len),
                r=r_base_map, Kcap=k_cap_map, beta=beta_map,
                N_min=10, N_max=90,
                p_wgd=p_wgd_map,
                boundary="drop",
                C_fn=C_fn, f_param_fn=_f, t_eval=TIMES, max_step=ODE_STEP,
                renormalize_M=True,
                n_tr=3, k_tr=1.0, k_kill=1.0,
            )
            T_mat_full = x0_full + x1_full
            n_tp = T_mat_full.shape[1]
            for t_idx in range(1, n_tp):
                t = day + TIMES[t_idx]
                col = T_mat_full[:, t_idx]
                total = float(col.sum())
                baseline_burden_timeline.append((t, total))
                ploidy_dict_t = {int(n): float(c) for n, c in zip(Ns_full, col) if c > 0}
                baseline_ploidy_snapshots.append((t, ploidy_dict_t))

            day += cycle_len
            ploidy_state = {int(n): float(c)
                    for n, c in zip(Ns_full, T_mat_full[:, -1]) if c > 0}

    all_bins: set[float] = set()
    ploidy_bin_data: list[tuple[float, dict[float, float], dict[float, float]]] = []

    for snap_day, ploidy_dict in baseline_ploidy_snapshots:
        total_cells = sum(ploidy_dict.values())
        if total_cells <= 0:
            ploidy_bin_data.append((snap_day, {}, {}))
            continue
        binned_counts: dict[float, float] = {}
        for chr_count, cell_count in ploidy_dict.items():
            ratio = chr_count / HAPLOID_N
            rounded_bin = round(round(ratio * 10) / 10, 1)  # Round to nearest 0.1
            binned_counts[rounded_bin] = binned_counts.get(rounded_bin, 0.0) + cell_count
        frac_binned = {b: v / total_cells for b, v in binned_counts.items()}
        ploidy_bin_data.append((snap_day, binned_counts, frac_binned))
        all_bins.update(binned_counts.keys())

    sorted_bins = sorted(all_bins)
    snap_days_arr = np.array([d[0] for d in ploidy_bin_data])
    cmap = plt.get_cmap("tab20")

    count_matrix = np.zeros((len(sorted_bins), len(ploidy_bin_data)))
    frac_matrix = np.zeros((len(sorted_bins), len(ploidy_bin_data)))
    for col_idx, (_, counts, fracs) in enumerate(ploidy_bin_data):
        for row_idx, b in enumerate(sorted_bins):
            count_matrix[row_idx, col_idx] = counts.get(b, 0.0)
            frac_matrix[row_idx, col_idx] = fracs.get(b, 0.0)

    drug_colors: dict[str, str] = {
        "gemcitabine": "orange",
        "bay1895344": "red",
        "alisertib": "green",
        "ispinesib": "blue",
        "none": "yellow",
    }


    def _add_drug_shading(ax, path_info, start_day: float = 0.0) -> None:
        """Shade the x-axis by drug cycle """
        if not path_info:
            return
        current_time = float(start_day)
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

    burden_days = np.array([t for t, _ in baseline_burden_timeline])
    burden_values = np.array([b for _, b in baseline_burden_timeline])

    fig1, ax1 = plt.subplots(figsize=(10, 5))

    # Drug-cycle shading (drawn first so data lines sit on top)
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
    ax1.set_ylabel("Number of cells (log scale)")
    ax1.set_title(f"Tumor Burden — Baseline Path\n({SAMPLE_NAME})")
    ax1.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1))
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig(f"{RESULTS_DIR}/baseline_tumor_burden.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: {RESULTS_DIR}/baseline_tumor_burden.png")

    fig1b, ax_combined = plt.subplots(figsize=(12, 6))

    # Drug-cycle shading
    _add_drug_shading(ax_combined, _baseline_path_info, start_day=baseline_start_day)

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

    ax_combined.set_yscale("log")
    ax_combined.set_xlabel("Day")
    ax_combined.set_ylabel("Number of cells (log scale)")
    ax_combined.legend(title="Treatment / Ploidy", loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=6, ncol=2)
    ax_combined.grid(True, alpha=0.3)

    fig1b.suptitle(
        f"Tumor Burden & Ploidy Cell Counts — Baseline Path\n"
        f"({SAMPLE_NAME}  |  ploidy = chromosomes / {HAPLOID_N}, rounded to 0.1)",
        fontsize=10)
    fig1b.tight_layout()
    fig1b.savefig(f"{RESULTS_DIR}/baseline_burden_and_ploidy_counts.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: {RESULTS_DIR}/baseline_burden_and_ploidy_counts.png")

    fig2, ax2 = plt.subplots(figsize=(11, 6))

    # Drug-cycle shading
    _add_drug_shading(ax2, _baseline_path_info, start_day=baseline_start_day)

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
    ax2.legend(title="Treatment / Ploidy", loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=6, ncol=2)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(f"{RESULTS_DIR}/baseline_ploidy_distribution.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: {RESULTS_DIR}/baseline_ploidy_distribution.png")

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
    _add_drug_shading(ax3, _baseline_path_info, start_day=baseline_start_day)

    ax3.plot(snap_days_arr, avg_ploidy_arr,
             color="darkorchid", linewidth=2.0, label="Mean ploidy")

    ax3.set_xlabel("Day")
    ax3.set_ylabel(f"Average ploidy (chromosomes / {HAPLOID_N})")
    ax3.set_title(
        f"Average Ploidy Over Treatment — Baseline Path\n({SAMPLE_NAME})"
    )
    ax3.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1.02, 1))
    ax3.grid(True, alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(f"{RESULTS_DIR}/baseline_average_ploidy.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: {RESULTS_DIR}/baseline_average_ploidy.png")

    print("\nGenerating Gibbs fitting animation...")

    # Subsample iterations to ~200 frames
    N_ANIM_FRAMES = min(200, len(all_iter_params))
    frame_indices = np.linspace(0, len(all_iter_params) - 1,
                                N_ANIM_FRAMES, dtype=int)

    # Pre-compute burden timelines for each frame
    anim_timelines: list[tuple[np.ndarray, np.ndarray, dict]] = []
    for fi, idx in enumerate(frame_indices):
        p = all_iter_params[idx]
        tl = _simulate_burden_timeline(
            OBSERVED_DRUGS_ADMINISTERED,
            p["r_base"], p["k_cap"],
            pk_state=p["pk"], beta=p["beta"], p_wgd=p["p_wgd"],
        )
        if tl is not None:
            days = np.array([t for t, _ in tl])
            vals = np.array([b for _, b in tl])
            anim_timelines.append((days, vals, p))
        else:
            anim_timelines.append((np.array([]), np.array([]), p))
        if (fi + 1) % 50 == 0 or fi == len(frame_indices) - 1:
            print(f"  Simulated frame {fi + 1}/{N_ANIM_FRAMES}")

    # Build animation
    fig_anim, ax_anim = plt.subplots(figsize=(12, 6))
    fig_anim.subplots_adjust(top=0.85, right=0.82)

    # Fixed observed data points
    obs_days = sorted(OBSERVED_TUMOR_BURDENS.keys())
    obs_vals = [OBSERVED_TUMOR_BURDENS[d] for d in obs_days]

    ax_anim.set_yscale("log")
    ax_anim.set_xlabel("Day")
    ax_anim.set_ylabel("Number of cells (log scale)")
    ax_anim.grid(True, alpha=0.3)

    # Set fixed axis limits from the data
    all_burden_vals = [v for _, vals, _ in anim_timelines for v in vals if v > 0]
    if all_burden_vals and obs_vals:
        y_lo = min(min(all_burden_vals), min(obs_vals)) * 0.1
        y_hi = max(max(all_burden_vals), max(obs_vals)) * 10
    else:
        y_lo, y_hi = 1e3, 1e12
    all_day_vals = [d for days, _, _ in anim_timelines for d in days]
    if all_day_vals:
        x_hi = max(max(all_day_vals), max(obs_days)) * 1.05
    else:
        x_hi = max(obs_days) * 1.05
    ax_anim.set_xlim(-1, x_hi)
    ax_anim.set_ylim(y_lo, y_hi)

    # Static scatter (observed)
    ax_anim.scatter(obs_days, obs_vals, color="firebrick", zorder=5,
                    s=50, label="Observed")

    # Current fit line
    fit_line, = ax_anim.plot([], [], color="steelblue", linewidth=2.0,
                             label="Predicted (current iter)")

    title_text = ax_anim.set_title("", fontsize=11, pad=12)
    ax_anim.legend(loc="upper left", bbox_to_anchor=(1.02, 1))

    # Store ghost artists for layering
    ghost_artists = []


    def _anim_init():
        fit_line.set_data([], [])
        title_text.set_text("")
        return [fit_line, title_text]


    def _anim_update(frame_num):
        days, vals, p = anim_timelines[frame_num]
        iter_idx = p["iter"]
        phase = "burn-in" if p["burnin"] else "sampling"

        # Add a ghost of the previous frame with iteration label
        if frame_num > 0:
            prev_days, prev_vals, prev_p = anim_timelines[frame_num - 1]
            if len(prev_days) > 0:
                gl, = ax_anim.plot(prev_days, prev_vals,
                                   color="steelblue", alpha=0.08,
                                   linewidth=0.6, zorder=1)
                ghost_artists.append(gl)
                # Label every 10th ghost with its iteration number
                if frame_num % 10 == 0:
                    lbl = ax_anim.annotate(
                        str(prev_p["iter"] + 1),
                        xy=(prev_days[-1], prev_vals[-1]),
                        fontsize=5, color="steelblue", alpha=0.35,
                        ha="left", va="center",
                        xytext=(3, 0), textcoords="offset points",
                    )
                    ghost_artists.append(lbl)

        # Current line
        if len(days) > 0:
            fit_line.set_data(days, vals)
        else:
            fit_line.set_data([], [])

        title_text.set_text(
            f"Gibbs Fitting — Iteration {iter_idx + 1}/{len(all_iter_params)} "
            f"[{phase}]\n"
            f"r={p['r_base']:.4f}  K={p['k_cap']:.2e}  "
            f"β={p['beta']:.4g}  p_wgd={p['p_wgd']:.2e}  logP={p['logP']:.1f}"
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
