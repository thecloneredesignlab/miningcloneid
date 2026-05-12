
import math
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import least_squares
from scipy.stats import linregress

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]  

EXP_BASE = PROJECT_ROOT / "data" / "InVitroData_Gemcitabine"

PATHS = {
    "counts_raw":      EXP_BASE / "processed" / "counts_by_well_time.parquet",
    "counts_agg":      EXP_BASE / "processed" / "counts_by_well_time_wellAggregated.parquet",
    "platemap":       EXP_BASE / "Gemcitabine_PlateMap_20240111.xlsx",
    "pkpd_constants": EXP_BASE / "drugKinetics" / "GemcitabineExposure_PKPD.xlsx"
}

# Unit / scale constants used by the in vitro gemcitabine model.
#
# The local drugKinetics README documents the heuristic relation
# `C_peak = 0.00971 * drug_intake` for an "intracellular concentration" proxy.
# This script has historically used that same factor as the dose-to-exposure
# signal conversion for the ODE driving term. The exact compartment/assay
# interpretation of that signal is not fully resolved from repo metadata, so the
# code treats it conservatively as a model exposure-signal conversion factor and
# preserves the existing numerical behavior.
MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM = 0.00971

# PK dFdCTP measurements are reported in ng/mL in the workbook. Dividing by the
# molecular weight converts ng/mL to nmol/mL, which is numerically equal to uM.
DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL = 503.1

for name, path in PATHS.items():
    if not path.exists():
        print(f"Warning: {name} not found at {path}")


def slugify_label(value: str) -> str:
    """Converts labels to filesystem-safe filename fragments."""
    return re.sub(r'[^A-Za-z0-9]+', '_', value.strip()).strip('_').lower()


#################### eta, K_decay fitting utilities ####################

def parse_pk_concentration_value(value, censored_strategy: str = "nan") -> Tuple[float, bool]:
    """
    Parses PK concentration values while preserving censoring semantics.

    Returns `(numeric_value, was_censored)`. Missing spreadsheet values
    (`np.nan`, blank cells) return `(np.nan, False)`. Censored strings such as
    `BDL`, `BQL`, `BLQ`, `N/F`, and `<0.5` return `was_censored=True`.
    """
    if censored_strategy not in {"nan", "zero", "half_lod"}:
        raise ValueError(
            "censored_strategy must be one of {'nan', 'zero', 'half_lod'}, "
            f"got {censored_strategy}"
        )

    if pd.isna(value):
        return np.nan, False

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), False

    value_str = str(value).strip()
    if value_str == "":
        return np.nan, False

    try:
        return float(value_str), False
    except ValueError:
        pass

    normalized = value_str.upper()
    censored_tokens = {
        "BDL",
        "BQL",
        "BLQ",
        "N/F",
        "ND",
        "N/D",
        "<LLOQ",
        "<LOQ",
        "BELOW DETECTION",
        "BELOW QUANTIFICATION",
    }
    missing_tokens = {"NA", "N/A", "NONE", "NULL"}

    if normalized in missing_tokens:
        return np.nan, False

    threshold_match = re.search(r'<\s*([0-9]+(?:\.[0-9]+)?)', value_str)
    if threshold_match:
        threshold = float(threshold_match.group(1))
        if censored_strategy == "half_lod":
            return threshold / 2.0, True
        if censored_strategy == "zero":
            return 0.0, True
        return np.nan, True

    if normalized in censored_tokens:
        if censored_strategy == "zero":
            return 0.0, True
        return np.nan, True

    return np.nan, False

def import_and_clean_pkpd(censored_strategy: str = "nan"):
    """
    Imports PK sheets and parses analyte columns with explicit censored-value handling.

    `censored_strategy="nan"` is the conservative default: BDL/BQL/BLQ-style
    entries become `np.nan` and are ignored by downstream means/fits. Set
    `"zero"` to reproduce the prior zero-imputation behavior explicitly, or
    `"half_lod"` to impute half the reported threshold for `<x` values.
    """
    path = PATHS["pkpd_constants"]
    if not path.exists():
        sys.exit(f"TERMINATING: PKPD file not found at: {path}")

    raw_sheets = pd.read_excel(path, sheet_name=None)
    clean_data_dict = {}
    
    for sheet_name, df in raw_sheets.items():
        df.columns = [re.sub(r'\s+', ' ', str(c).strip()) for c in df.columns]
        
        if 'Timepoint' not in df.columns and 'Sample' in df.columns:
            time_regex, rep_regex = r'-(\d+\.?\d*)h', r'h-(\d+)'
            df['Timepoint'] = df['Sample'].apply(lambda x: float(re.search(time_regex, str(x)).group(1)) if re.search(time_regex, str(x)) else 0.0)
            df['replicate'] = df['Sample'].apply(lambda x: int(re.search(rep_regex, str(x)).group(1)) if re.search(rep_regex, str(x)) else 1)
        
        analyte_cols = [c for c in df.columns if '(ng/mL)' in c]
        for col in analyte_cols:
            parsed_values = df[col].apply(
                lambda value: parse_pk_concentration_value(
                    value,
                    censored_strategy=censored_strategy,
                )
            )
            df[col] = parsed_values.apply(lambda item: item[0]).astype(float)
            df[f"{col}__was_censored"] = parsed_values.apply(lambda item: item[1]).astype(bool)
            
        clean_data_dict[sheet_name] = df
    return clean_data_dict

def get_r_eta_parameters(sheets, dfdctp_molecular_weight_ng_per_nmol=DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL):
    """
    Calculates uptake and decay parameters for each cohort (2N, 4N).

    eta:
        Approximate uptake/conversion rate from extracellular gemcitabine dose
        to intracellular dFdCTP accumulation, estimated from the 0h-to-1h increase.

    k_decay:
        Exponential intracellular dFdCTP decay rate, estimated from the post-peak
        decline and converted from 1/hour to 1/day.
    """
    results = {}
    for cohort in ['4N', '2N']:
        if cohort in sheets:
            k_val, t_half, r2 = extract_half_life(sheets[cohort])
            eta_val = extract_eta(
                sheets[cohort],
                dose_muM=1.0,
                dfdctp_molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
            )
            
            results[cohort] = {
                'eta': eta_val,
                'k_decay': k_val * 24 if k_val else 0.0, # Convert 1/hr to 1/day
                'hl_hours': t_half,
                'r2_confidence': r2
            }
    return results

def extract_half_life(df: pd.DataFrame, analyte: str = 'dFdCTP (ng/mL)') -> Tuple[Optional[float], Optional[float], Optional[float]]:
    mean_data = df.groupby('Timepoint')[analyte].mean().reset_index()
    mean_data = mean_data.dropna(subset=[analyte])
    mean_data = mean_data[mean_data[analyte] > 0]
    if len(mean_data) < 2:
        warnings.warn(
            f"Insufficient non-missing positive {analyte} data to estimate half-life.",
            RuntimeWarning,
        )
        return None, None, None

    peak_idx = mean_data[analyte].idxmax()
    decay_phase = mean_data[mean_data['Timepoint'] >= mean_data.loc[peak_idx, 'Timepoint']]
    if len(decay_phase) < 2:
        warnings.warn(
            f"Insufficient post-peak {analyte} data to estimate half-life.",
            RuntimeWarning,
        )
        return None, None, None

    slope, _, r_val, _, _ = linregress(decay_phase['Timepoint'], np.log(decay_phase[analyte]))
    k = -slope
    if k <= 0:
        warnings.warn(f"Non-positive half-life slope estimated for {analyte}.", RuntimeWarning)
        return None, None, None 

    return round(k, 4), round(np.log(2)/k, 2), round(r_val**2, 4)

def converted_extracellular_gem_signal(
    dose_muM: float,
    conversion_factor: float = MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM,
) -> float:
    """
    Converts nominal extracellular gemcitabine dose (uM) to the model's
    exposure-signal scale.
    """
    if dose_muM <= 0:
        raise ValueError(f"dose_muM must be > 0, got {dose_muM}")
    if conversion_factor <= 0:
        raise ValueError(f"conversion_factor must be > 0, got {conversion_factor}")
    return conversion_factor * dose_muM

def dfdctp_ng_per_ml_to_uM(
    dfdctp_ng_per_ml: float,
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> float:
    """Converts dFdCTP concentration from ng/mL to uM using its molecular weight."""
    if molecular_weight_ng_per_nmol <= 0:
        raise ValueError(
            "molecular_weight_ng_per_nmol must be > 0, "
            f"got {molecular_weight_ng_per_nmol}"
        )
    return dfdctp_ng_per_ml / molecular_weight_ng_per_nmol

def estimate_eta_per_day(
    delta_dfdctp_signal_per_hour: float,
    dose_muM: float,
    conversion_factor: float = MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM,
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> float:
    """
    Estimates eta from the observed 0h->1h dFdCTP buildup.

    Formula preserved from the original implementation:
      eta [day^-1] =
        (delta_dFdCTP_ng_per_mL_per_hour / molecular_weight_ng_per_nmol)
        / converted_extracellular_gem_signal(dose_muM)
        * 24

    Under the model convention below, the converted exposure signal and the ODE
    intracellular state are treated on the same uM-equivalent signal scale, so
    eta has units of day^-1.
    """
    exposure_signal = converted_extracellular_gem_signal(
        dose_muM,
        conversion_factor=conversion_factor,
    )
    intracellular_signal_delta_uM_per_hour = dfdctp_ng_per_ml_to_uM(
        delta_dfdctp_signal_per_hour,
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    if exposure_signal == 0:
        raise ValueError("converted extracellular gemcitabine signal must be nonzero")
    return (intracellular_signal_delta_uM_per_hour / exposure_signal) * 24.0

def extract_eta(
    df: pd.DataFrame,
    dose_muM: float = 1.0,
    analyte: str = 'dFdCTP (ng/mL)',
    dfdctp_molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[float]:
    """
    Estimates eta from the 0h->1h increase in PK-measured dFdCTP.

    The workbook reports dFdCTP in ng/mL. This function converts that change to
    a uM-equivalent intracellular signal using the dFdCTP molecular weight, then
    divides by the model's converted extracellular gemcitabine exposure signal.
    The returned eta therefore preserves the previous numeric behavior and is
    interpreted here as day^-1 under the model's shared signal-scale convention.
    """
    mean_data = df.groupby('Timepoint')[analyte].mean()
    if 0.0 not in mean_data.index or 1.0 not in mean_data.index:
        warnings.warn(
            f"Missing required 0h/1h {analyte} means; eta cannot be estimated.",
            RuntimeWarning,
        )
        return None

    delta_c_hourly = mean_data[1.0] - mean_data[0.0]
    if pd.isna(delta_c_hourly):
        warnings.warn(
            f"Observed 0h->1h delta for {analyte} is NaN; eta cannot be estimated.",
            RuntimeWarning,
        )
        return None
    eta_daily = estimate_eta_per_day(
        delta_dfdctp_signal_per_hour=delta_c_hourly,
        dose_muM=dose_muM,
        conversion_factor=MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM,
        molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
    )
    return round(eta_daily, 4)

#################### r, K fitting utilities ####################

def parse_gem_label(val):
    """Extracts standardized Gemcitabine concentration from raw text."""
    if pd.isna(val):
        return np.nan
    val_str = str(val)
    match = re.search(r'(?i)\b([0-9]+\.?[0-9]*)\s*(nM|uM|µM|pM)\b', val_str)
    if not match:
        return np.nan
    
    amount, unit = match.groups()
    unit = unit.lower().replace('µm', 'uM').replace('um', 'uM').replace('nm', 'nM')
    return f"{amount} {unit}"

def load_and_clean_counts(max_days: float = 5.0) -> pd.DataFrame:
    """Loads the aggregated parquet counts data, cleans columns, and filters time."""
    counts_path = PATHS["counts_agg"]
    if not counts_path.exists():
        raise FileNotFoundError(f"Expected counts file at {counts_path}")
        
    df = pd.read_parquet(counts_path)
    
    df = df.rename(columns={"Classifier.Phenotype": "phenotype", "N": "count"})
    
    if 'plate_row' not in df.columns or 'plate_col' not in df.columns:
        extracted = df['well'].str.extract(r'^([A-H])([0-9]{1,2})(?:_([0-9]+))?$')
        df['plate_row'] = extracted[0]
        df['plate_col'] = pd.to_numeric(extracted[1], errors='coerce').astype('Int64')
        
    if 'time_hours' not in df.columns:
        raise ValueError("Expected 'time_hours' in counts table.")
    
    df['time_days'] = df['time_hours'] / 24.0
    df = df[df['time_days'] <= max_days].copy()
    
    df = df[df['phenotype'].isin(['Alive', 'Dead', 'Transitional'])]
    
    df = df.dropna(subset=['plate_row', 'plate_col'])
    df = df[(df['plate_col'] >= 1) & (df['plate_col'] <= 12)]
    
    return df

def load_platemap() -> pd.DataFrame:
    """Loads and parses the Excel platemap to extract ploidy and drug conditions."""
    pm_path = PATHS["platemap"]
    if not pm_path.exists():
         raise FileNotFoundError(f"Expected platemap file at {pm_path}")
        
    pm = pd.read_excel(pm_path, sheet_name=0)
    cols_lower = [str(c).lower() for c in pm.columns]
    
    plate_long = pd.DataFrame()
    
    if 'well' in cols_lower or 'wellid' in cols_lower or ('row' in cols_lower and any(c in cols_lower for c in ['col', 'column'])):
        pass
        
    else:
        row_col = pm.columns[0] 
        pm = pm.rename(columns={row_col: "Row"})
        pm['Row'] = pm['Row'].astype(str).str.upper()
        pm = pm[pm['Row'].isin(list("ABCDEFGH"))]
        
        val_vars = [c for c in pm.columns if str(c).isdigit()]
        plate_long = pd.melt(pm, id_vars="Row", value_vars=val_vars, 
                             var_name="Col", value_name="condition_raw")
                             
        plate_long['plate_row'] = plate_long['Row']
        plate_long['plate_col'] = pd.to_numeric(plate_long['Col'], errors='coerce')
        
        plate_long['ploidy'] = plate_long['condition_raw'].str.extract(r'(?i)\b([24]N)\b')[0].str.upper()
        
        plate_long['gem'] = plate_long['condition_raw'].apply(parse_gem_label)
        
    return plate_long[['plate_row', 'plate_col', 'ploidy', 'gem']].dropna(subset=['plate_row', 'plate_col'])

def assemble_modeling_dataset() -> pd.DataFrame:
    """Combines counts and platemap into a single dataframe."""
    counts_df = load_and_clean_counts()
    platemap_df = load_platemap()
    
    merged = pd.merge(counts_df, platemap_df, on=['plate_row', 'plate_col'], how='left')
    return merged

def get_fitting_data(df: pd.DataFrame, gem_dose: str = "0 nM", ploidy: str = "2N", phenotype: str = "Alive"):
    """
    Extracts time and biological replicates for specific conditions to use in scipy.optimize.
    Returns: time_days (1D array), counts (2D array: rows=time, cols=replicate wells)
    """
    subset = df[(df['gem'] == gem_dose) & 
                (df['ploidy'] == ploidy) & 
                (df['phenotype'] == phenotype)]
                
    if subset.empty:
        raise ValueError(f"No data found for {gem_dose}, {ploidy}, {phenotype}")
        
    pivot = subset.pivot_table(index='time_days', columns=['plate_row', 'plate_col'], values='count')
    
    t_data = pivot.index.values
    y_data = pivot.values # Shape: (n_timepoints, n_replicates)
    
    return t_data, y_data

def logistic_growth(t, N0, r, K):
    """Analytical solution to the logistic growth equation."""
    return (K * N0) / (N0 + (K - N0) * np.exp(-r * t))

# Canonical internal unit convention for the joint in vitro ODE:
# - time `t` is in days
# - `dose_muM` is nominal extracellular gemcitabine dose in uM
# - `C_ext_signal` is the converted exposure signal derived from `dose_muM`
# - `C_dna_signal` is the intracellular dFdCTP state on the same signal scale
#   used to estimate eta from the PK data
# - `eta`, `k_decay`, `k_tr`, and `k_clear` are in day^-1
# - `k_kill` is in day^-1 per unit intracellular dFdCTP signal, so
#   `kappa = k_kill * delayed_C_dna_signal` is a death hazard in day^-1

def decay_rate_from_half_life_days(half_life_days: float) -> float:
    """Converts an extracellular half-life in days to a decay rate in days^-1."""
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days

def extracellular_gemcitabine_concentration(t, dose_muM, k_ext_decay_per_day):
    """Returns the model's exponentially decaying extracellular gemcitabine exposure signal."""
    initial_exposure_signal = converted_extracellular_gem_signal(dose_muM)
    return initial_exposure_signal * np.exp(-k_ext_decay_per_day * t)


def residuals_fixed_N0(params, t, y_true, N0_fixed):
    """Objective function that only fits r and K, keeping N0 fixed."""
    r, K = params
    y_pred = logistic_growth(t, N0_fixed, r, K)
    return y_true - y_pred

def logistic_growth_replicate_n0(t, n0_by_rep, r, K):
    """Evaluates logistic growth with shared r/K and replicate-specific initial counts."""
    t_arr = np.asarray(t, dtype=float)
    n0_arr = np.asarray(n0_by_rep, dtype=float)
    return logistic_growth(t_arr[:, None], n0_arr[None, :], r, K)

def residuals_shared_rk_replicate_n0(params, t_data, y_alive, n0_by_rep, t_max=None):
    """Residuals for shared r/K baseline growth with replicate-specific fixed N0."""
    r, K = params
    y_pred = logistic_growth_replicate_n0(t_data, n0_by_rep, r, K)
    valid_mask = ~np.isnan(y_alive)
    if t_max is not None:
        valid_mask = valid_mask & (np.asarray(t_data)[:, None] <= t_max)
    return (y_alive - y_pred)[valid_mask]

def fit_baseline_shared_rk_replicate_n0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None):
    """Fits shared cohort r/K while fixing each replicate to its observed baseline N0."""
    y_alive_matrix = y_alive[:, None] if y_alive.ndim == 1 else y_alive
    n0_by_replicate = y_alive_matrix[0, :].astype(float)
    replicate_mask = ~np.isnan(n0_by_replicate)
    if not np.any(replicate_mask):
        raise ValueError(f"No valid baseline alive counts found for {ploidy_label}")

    y_alive_matrix = y_alive_matrix[:, replicate_mask]
    n0_by_replicate = n0_by_replicate[replicate_mask]

    residual_mask = ~np.isnan(y_alive_matrix)
    if t_max is not None:
        residual_mask = residual_mask & (np.asarray(t_data)[:, None] <= t_max)
    n_observations = int(np.count_nonzero(residual_mask))
    if n_observations < 2:
        raise ValueError(f"Not enough alive observations to fit baseline growth for {ploidy_label}")

    K_guess = np.nanmax(y_alive_matrix[residual_mask])
    r_guess = 0.5

    res = least_squares(
        residuals_shared_rk_replicate_n0,
        [r_guess, K_guess],
        args=(t_data, y_alive_matrix, n0_by_replicate, t_max),
        bounds=([0, 0], [np.inf, np.inf]),
    )
    r_opt, K_opt = res.x
    n0_mean = float(np.nanmean(n0_by_replicate))
    n0_min = float(np.nanmin(n0_by_replicate))
    n0_max = float(np.nanmax(n0_by_replicate))

    print(f"\n{'='*60}")
    print(f"{ploidy_label} baseline fit")
    print("replicate-specific N0 enabled")
    if t_max is not None:
        print(f"Fitted only on t <= {t_max} Days")
    print(f"Replicates/wells used: {len(n0_by_replicate)}")
    print(f"N0 mean: {n0_mean:.1f}")
    print(f"N0 range: {n0_min:.1f} to {n0_max:.1f}")
    print(f"Observations used: {n_observations}")
    print(f"r = {r_opt:.4f} days^-1")
    print(f"K = {K_opt:.1f}")
    print("=" * 60 + "\n")

    fig, ax1 = plt.subplots(figsize=(8, 5))
    color1 = 'tab:purple'
    ax1.set_xlabel('Time (Days)')
    ax1.set_ylabel('Alive Cell Count', color=color1)

    for rep_idx in range(y_alive_matrix.shape[1]):
        ax1.plot(
            t_data,
            y_alive_matrix[:, rep_idx],
            'o',
            color=color1,
            alpha=0.18,
            label='Alive (Replicates)' if rep_idx == 0 else "",
        )

    ax1.plot(t_data, np.nanmean(y_alive_matrix, axis=1), 's', color='indigo', label='Alive (Mean)')

    t_smooth = np.linspace(min(t_data), max(t_data), 100)
    fitted_curves = logistic_growth_replicate_n0(t_smooth, n0_by_replicate, r_opt, K_opt)
    for rep_idx in range(fitted_curves.shape[1]):
        ax1.plot(
            t_smooth,
            fitted_curves[:, rep_idx],
            '-',
            color='fuchsia',
            alpha=0.25,
            linewidth=1.2,
            label='Replicate-specific Fits' if rep_idx == 0 else "",
        )
    ax1.plot(
        t_smooth,
        logistic_growth(t_smooth, n0_mean, r_opt, K_opt),
        '-',
        color='darkmagenta',
        linewidth=2.2,
        label=f'Mean-N0 Summary Fit (r={r_opt:.3f})',
    )

    if t_max is not None:
        ax1.axvline(x=t_max, color='gray', linestyle='--', alpha=0.7, label=f'Fit Cutoff (t={t_max})')

    ax1.tick_params(axis='y', labelcolor=color1)

    plt.title(f"Phase 1 Control Fit: {ploidy_label} (Replicate-specific N0)")
    plt.legend(loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_dir is not None:
        output_path = output_dir / f"baseline_replicate_n0_{slugify_label(ploidy_label)}.png"
        fig.savefig(output_path, dpi=200, bbox_inches='tight')

    return {
        "N0_by_replicate": n0_by_replicate,
        "N0_mean": n0_mean,
        "r": float(r_opt),
        "K": float(K_opt),
        "ploidy": ploidy_label,
        "n_replicates": int(len(n0_by_replicate)),
        "n_observations": n_observations,
    }

def fit_baseline_locked_N0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None):
    """Backward-compatible wrapper for the replicate-specific-N0 baseline fit."""
    return fit_baseline_shared_rk_replicate_n0(
        t_data,
        y_alive,
        y_dead,
        ploidy_label,
        t_max=t_max,
        output_dir=output_dir,
    )



#################### Cohort fitting utilities ####################

def single_clone_signal_ode_joint(y, t, r, K, dose_muM, eta, k_decay, k_ext_decay_per_day, k_tr, k_kill, k_clear, n_tr):
    """
    ODE system tracking:
    - A: Alive cells (y[0])
    - C_dna_signal: intracellular dFdCTP signal (y[1])
    - Z_1 ... Z_n: Transit compartments (y[2 : 2+n_tr])
    - D_obs: Observable Dead cells (y[-1])
    """
    A = y[0]
    C_dna_signal = y[1]
    D_obs = y[-1]

    # Extracellular gemcitabine starts at the converted concentration and then
    # decays exponentially with a user-specified half-life/rate. Intracellular
    # dFdCTP still follows the existing uptake/decay equation below.
    C_ext_signal = extracellular_gemcitabine_concentration(t, dose_muM, k_ext_decay_per_day)

    # Active drug kinetics: eta * C_ext_signal and k_decay * C_dna_signal both
    # have units of signal/day, so dC_dna_signal/dt does as well.
    dC_dna_signal = eta * C_ext_signal - k_decay * C_dna_signal
    
    # Transit compartments
    dZ = np.zeros(n_tr)
    if n_tr > 0:
        delayed_c_dna_signals = y[2:2+n_tr]
        dZ[0] = k_tr * (C_dna_signal - delayed_c_dna_signals[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (delayed_c_dna_signals[i-1] - delayed_c_dna_signals[i])
            
        kappa = k_kill * delayed_c_dna_signals[-1]
    else:
        kappa = k_kill * C_dna_signal
        
    # Alive cells: Logistic growth minus drug-induced death
    dA = r * A * (1 - A / K) - kappa * A
    
    # Observable dead cells: Accumulate from death, deplete via lysis/clearance
    dD_obs = kappa * A - k_clear * D_obs
    
    return [dA, dC_dna_signal] + dZ.tolist() + [dD_obs]


def simulate_joint_ext(t, params_3, N0, D0, r, K, dose_muM, eta_fixed, k_decay_fixed, k_ext_decay_per_day, n_tr):
    """
    Wrapper to simulate the ODE system. 
    params_3 = [k_tr, k_kill, k_clear]
    """
    k_tr, k_kill, k_clear = params_3 
    
    # Initial conditions: [A0, C_dna_0, Z_1_0, ..., Z_n_0, D_obs_0]
    y0 = [N0, 0.0] + [0.0] * n_tr + [D0] 
    
    sol = odeint(single_clone_signal_ode_joint, y0, t, 
                 args=(r, K, dose_muM, eta_fixed, k_decay_fixed, k_ext_decay_per_day, k_tr, k_kill, k_clear, n_tr))
    
    return sol[:, 0], sol[:, -1]


def residuals_global_joint(params_3, dose_data_list, r_fixed, K_fixed, eta_fixed, k_decay_fixed, k_ext_decay_per_day, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0):
    all_residuals = []
    
    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']
        
        y_alive_pred, y_dead_pred = simulate_joint_ext(
            t_data, params_3, N0, D0, r_fixed, K_fixed, dose_muM, eta_fixed, k_decay_fixed, k_ext_decay_per_day, n_tr_test
        )
        
        # Option 1: Fit mean trajectories first for stability
        if fit_means_only:
            y_alive_obs = np.nanmean(y_alive_data, axis=1) if y_alive_data.ndim > 1 else y_alive_data
            y_dead_obs = np.nanmean(y_dead_data, axis=1) if y_dead_data.ndim > 1 else y_dead_data
            
            scale_alive = np.nanmax(y_alive_obs) if np.nanmax(y_alive_obs) > 0 else 1.0
            scale_dead = np.nanmax(y_dead_obs) if np.nanmax(y_dead_obs) > 0 else 1.0
            
            res_alive = (y_alive_obs - y_alive_pred) / scale_alive
            res_dead = (y_dead_obs - y_dead_pred) / scale_dead
        
        # Option 2: Fit all replicate points
        else:
            scale_alive = np.nanmax(y_alive_data) if np.nanmax(y_alive_data) > 0 else 1.0
            scale_dead = np.nanmax(y_dead_data) if np.nanmax(y_dead_data) > 0 else 1.0
            
            res_alive = (y_alive_data.flatten() - np.repeat(y_alive_pred, y_alive_data.shape[1])) / scale_alive if y_alive_data.ndim > 1 else (y_alive_data - y_alive_pred) / scale_alive
            res_dead = (y_dead_data.flatten() - np.repeat(y_dead_pred, y_dead_data.shape[1])) / scale_dead if y_dead_data.ndim > 1 else (y_dead_data - y_dead_pred) / scale_dead
        
        # Gentler optional dose weighting
        weight = 1.0
        if dose_muM >= 0.2:
            weight = high_dose_weight
            
        all_residuals.extend(res_alive[~np.isnan(res_alive)] * weight)
        all_residuals.extend(res_dead[~np.isnan(res_dead)] * weight)
        
    return np.array(all_residuals)


def plot_global_fit_subplots_joint(
    dose_data_list,
    ploidy_label,
    r,
    K,
    n_tr,
    eta,
    k_decay_fixed,
    k_ext_decay_per_day,
    k_tr,
    k_kill,
    k_clear,
    output_dir: Optional[Path] = None
):
    n_doses = len(dose_data_list)
    cols = 3
    rows = math.ceil(n_doses / cols)
    
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), sharex=True)
    if n_doses == 1: axes = [axes]
    elif n_doses > 1: axes = axes.flatten()
        
    ode_params_3 = [k_tr, k_kill, k_clear]
    
    for i, data in enumerate(dose_data_list):
        ax = axes[i]
        t_data = data['t']
        y_alive_raw = data['y_alive']
        y_dead_raw = data['y_dead']
        N0 = data['N0']
        D0 = data['D0']
        dose_muM = data['dose_muM']
        dose_label = data['dose_label']
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            A_sim, D_sim = simulate_joint_ext(
                t_data, ode_params_3, N0, D0, r, K, dose_muM, eta, k_decay_fixed, k_ext_decay_per_day, n_tr
            )
            
        if y_alive_raw.ndim > 1:
            for j in range(y_alive_raw.shape[1]):
                ax.scatter(t_data, y_alive_raw[:, j], color='limegreen', alpha=0.3, s=20)
        else:
            ax.scatter(t_data, y_alive_raw, color='limegreen', alpha=0.5, s=20)
        ax.plot(t_data, A_sim, color='darkgreen', linewidth=2, label='Alive Model' if i==0 else "")
        
        if y_dead_raw.ndim > 1:
            for j in range(y_dead_raw.shape[1]):
                ax.scatter(t_data, y_dead_raw[:, j], color='lightcoral', alpha=0.3, s=20)
        else:
            ax.scatter(t_data, y_dead_raw, color='lightcoral', alpha=0.5, s=20)
        ax.plot(t_data, D_sim, color='darkred', linewidth=2, linestyle='--', label='Dead Model' if i==0 else "")
        
        ax.set_title(f"Dose: {dose_label}", fontsize=12, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.6)
        
        if i % cols == 0: ax.set_ylabel("Cell Count", fontsize=11)
        if i >= (rows - 1) * cols: ax.set_xlabel("Time (Days)", fontsize=11)
        if i == 0: ax.legend(loc='upper left', fontsize=10)
            
    for j in range(n_doses, len(axes)): fig.delaxes(axes[j])
        
    param_text = (f"Optima: n_tr = {n_tr} | $\\eta$ (fixed) = {eta:.3f} | $k_{{decay}}$ (fixed) = {k_decay_fixed:.3f}\n"
                  f"$k_{{tr}}$ = {k_tr:.3f} | $k_{{kill}}$ = {k_kill:.3f} | $k_{{clear}}$ = {k_clear:.3f}")
    fig.suptitle(f"{ploidy_label} Cohort - Joint Live/Dead Kinetic Fit\n{param_text}", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    if output_dir is not None:
        output_path = output_dir / f"cohort_joint_fit_{slugify_label(ploidy_label)}.png"
        fig.savefig(output_path, dpi=200, bbox_inches='tight')

def get_fitting_data_one_row(df, gem_dose: str, ploidy: str, phenotype: str, plate_row: str):
    """
    Extracts time and counts for one plate row replicate.

    Example:
        plate_row='A' for 2N replicate A
        plate_row='E' for 4N replicate E
    """
    subset = df[(df['gem'] == gem_dose) & 
                (df['ploidy'] == ploidy) & 
                (df['phenotype'] == phenotype) &
                (df['plate_row'] == plate_row)]
                
    if subset.empty:
        raise ValueError(f"No data found for {gem_dose}, {ploidy}, {phenotype}, row {plate_row}")
        
    pivot = subset.pivot_table(index='time_days', values='count')
    
    t_data = pivot.index.values
    y_data = pivot['count'].values
    
    return t_data, y_data

def fit_joint_one_replicate(dose_data_list, r_opt, K_opt, eta_fixed, k_decay_fixed, k_ext_decay_per_day, ploidy, rep_idx):
    """
    Fits one global joint live/dead model across all doses for one replicate.
    """
    best_cost, best_params, best_n_tr = np.inf, None, 1
    n_range = range(2, 10)
    
    guess_grid = [
        [0.2, 10.0, 0.2],
        [0.5, 25.0, 0.5],
        [1.0, 50.0, 1.0],
        [2.0, 100.0, 1.0],
        [5.0, 50.0, 0.5],
        [10.0, 20.0, 0.2],
    ]
    
    bounds = ([1e-4, 1e-4, 0.0], [50.0, 500.0, 5.0])
    
    for n_test in n_range:
        for guess in guess_grid:
            res = least_squares(residuals_global_joint, guess, bounds=bounds,
                                args=(dose_data_list, r_opt, K_opt, eta_fixed, k_decay_fixed,
                                      k_ext_decay_per_day, n_test, False, 1.0),
                                loss="soft_l1", f_scale=0.2, max_nfev=3000)
            
            if res.cost < best_cost:
                best_cost, best_params, best_n_tr = res.cost, res.x, n_test

    if best_params is None:
        return None

    k_tr_opt, k_kill_opt, k_clear_opt = best_params
    
    print(f"\n--- {ploidy} Replicate {rep_idx} ---")
    print(f"Optimal n_tr: {best_n_tr}")
    print(f"Optimal k_tr: {k_tr_opt:.4f}")
    print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
    print(f"Optimal k_clear: {k_clear_opt:.4f}")
    print(f"Best cost: {best_cost:.4f}")

    return {
        'ploidy': ploidy,
        'rep_idx': rep_idx,
        'n_tr': best_n_tr,
        'k_tr': k_tr_opt,
        'k_kill': k_kill_opt,
        'k_clear': k_clear_opt,
        'cost': best_cost
    }

def build_row_replicate_dose_data_list(df, gem_doses, ploidy, plate_row, t_max):
    """
    Builds a dose_data_list for one row replicate across all gemcitabine doses.
    """
    dose_data_list = []

    for gem_dose in gem_doses:
        try:
            t_part_A, y_part_A = get_fitting_data_one_row(
                df, gem_dose=gem_dose, ploidy=ploidy, phenotype="Alive", plate_row=plate_row
            )
            t_part_D, y_part_D = get_fitting_data_one_row(
                df, gem_dose=gem_dose, ploidy=ploidy, phenotype="Dead", plate_row=plate_row
            )
            
            time_mask = t_part_A <= t_max
            t_fit = t_part_A[time_mask]
            
            y_alive_fit = y_part_A[time_mask]
            y_dead_fit = y_part_D[time_mask]
            
            dose_muM = float(gem_dose.split()[0]) / 1000.0
            N0 = y_alive_fit[0]
            D0 = y_dead_fit[0]
            
            dose_data_list.append({
                'dose_label': gem_dose,
                'dose_muM': dose_muM,
                't': t_fit,
                'y_alive': y_alive_fit,
                'y_dead': y_dead_fit,
                'N0': N0,
                'D0': D0,
                'plate_row': plate_row
            })
            
        except ValueError as e:
            print(f"  Skipping {gem_dose}, row {plate_row}: {e}")

    return dose_data_list

def plot_replicate_parameter_summary(replicate_fit_df, output_dir: Optional[Path] = None):
    """
    Plots fitted kinetic parameters by ploidy, with one point per row replicate.
    """
    params_to_plot = ['k_tr', 'k_kill', 'k_clear']
    param_labels = {
        'k_tr': r'$k_{tr}$',
        'k_kill': r'$k_{kill}$',
        'k_clear': r'$k_{clear}$'
    }
    
    x_map = {'2N': 0, '4N': 1}
    
    for param in params_to_plot:
        fig, ax = plt.subplots(figsize=(6, 5))
        
        for ploidy in ['2N', '4N']:
            subset = replicate_fit_df[replicate_fit_df['ploidy'] == ploidy].copy()
            
            if subset.empty:
                continue
            
            x_center = x_map[ploidy]
            jitter = np.linspace(-0.08, 0.08, len(subset))
            x_vals = x_center + jitter
            y_vals = subset[param].values
            
            ax.scatter(x_vals, y_vals, s=80, alpha=0.8, label=ploidy)
            
            # Optional: label each point by row replicate
            if 'plate_row' in subset.columns:
                for x, y, row_label in zip(x_vals, y_vals, subset['plate_row'].values):
                    ax.text(x + 0.02, y, str(row_label), fontsize=9, va='center')
            
            # Plot mean line for each ploidy
            y_mean = np.nanmean(y_vals)
            ax.hlines(y_mean, x_center - 0.18, x_center + 0.18, colors='black', linewidth=2)
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['2N', '4N'], fontsize=12)
        ax.set_ylabel(param_labels[param], fontsize=13)
        ax.set_title(f"Per-Replicate Estimates of {param_labels[param]}", fontsize=14, fontweight='bold')
        ax.grid(True, axis='y', linestyle='--', alpha=0.4)
        plt.tight_layout()

        if output_dir is not None:
            output_path = output_dir / f"replicate_parameter_summary_{param}.png"
            fig.savefig(output_path, dpi=200, bbox_inches='tight')

if __name__ == "__main__":
    output_dir = SCRIPT_PATH.parent / "invitro_fitting_outputs" / datetime.now().strftime("%Y%m%dT%H%M%S")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving figures to: {output_dir}")

    gem_ext_half_life_days = 1.0
    k_ext_decay_per_day = decay_rate_from_half_life_days(gem_ext_half_life_days)
    pk_censored_strategy = "nan"

    dfdctp_molecular_weight_ng_per_nmol = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL
    pk_sheets = import_and_clean_pkpd(censored_strategy=pk_censored_strategy)
    fitted_map = get_r_eta_parameters(
        pk_sheets,
        dfdctp_molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
    )

    try:
        eta_2N, k_decay_2N = fitted_map['2N']['eta'], fitted_map['2N']['k_decay']
        eta_4N, k_decay_4N = fitted_map['4N']['eta'], fitted_map['4N']['k_decay']
        print(f"\n{'='*60}")
        print(f"Successfully extracted PK:")
        print(f"  censored PK strategy: {pk_censored_strategy}")
        print(
            "  2N -> "
            f"eta: {eta_2N} day^-1 "
            "(from dFdCTP ng/mL -> uM-equivalent signal), "
            f"k_decay: {k_decay_2N:.4f} day^-1"
        )
        print(
            "  4N -> "
            f"eta: {eta_4N} day^-1 "
            "(from dFdCTP ng/mL -> uM-equivalent signal), "
            f"k_decay: {k_decay_4N:.4f} day^-1"
        )
        print(f"{'='*60}")
    except KeyError as e:
        sys.exit(f"TERMINATING: Required cohort data {e} missing from results.")
        
    df = assemble_modeling_dataset()

    # Pull 2N Data
    t_2N, y_alive_2N = get_fitting_data(df, gem_dose="0 nM", ploidy="2N", phenotype="Alive")
    _, y_dead_2N = get_fitting_data(df, gem_dose="0 nM", ploidy="2N", phenotype="Dead")

    # Pull 4N Data
    t_4N, y_alive_4N = get_fitting_data(df, gem_dose="0 nM", ploidy="4N", phenotype="Alive")
    _, y_dead_4N = get_fitting_data(df, gem_dose="0 nM", ploidy="4N", phenotype="Dead")
    
    params_2N_baseline = fit_baseline_shared_rk_replicate_n0(
        t_2N, y_alive_2N, y_dead_2N, "2N Cohort", t_max=3.8, output_dir=output_dir
    )
    params_4N_baseline = fit_baseline_shared_rk_replicate_n0(
        t_4N, y_alive_4N, y_dead_4N, "4N Cohort", t_max=3.8, output_dir=output_dir
    )
    
    ploidy_options = [p for p in df['ploidy'].unique() if pd.notna(p)]
    gem_doses = sorted([d for d in df['gem'].unique() if pd.notna(d) and d != "0 nM"], key=lambda x: float(x.split()[0]))
    
    t_max = 5 

    for ploidy in ploidy_options:
        print(f"\n{'='*60}")
        print(f"3-PARAM FIT | COHORT: {ploidy}")
        print(f"{'='*60}")
        
        if ploidy == "2N":
            r_opt, K_opt = params_2N_baseline["r"], params_2N_baseline["K"]
            eta_fixed = eta_2N
            k_decay_fixed = k_decay_2N
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
            eta_fixed = eta_4N
            k_decay_fixed = k_decay_4N
        else:
            continue
            
        dose_data_list = []

        print(f"Gathering and truncating joint data...")
        for gem_dose in gem_doses:
            try:
                t_part_A, y_part_A = get_fitting_data(df, gem_dose=gem_dose, ploidy=ploidy, phenotype="Alive")
                t_part_D, y_part_D = get_fitting_data(df, gem_dose=gem_dose, ploidy=ploidy, phenotype="Dead")
                
                time_mask = t_part_A <= t_max
                t_fit = t_part_A[time_mask]
                
                y_alive_fit = y_part_A[time_mask, :] if y_part_A.ndim > 1 else y_part_A[time_mask]
                y_dead_fit = y_part_D[time_mask, :] if y_part_D.ndim > 1 else y_part_D[time_mask]
                
                dose_muM = float(gem_dose.split()[0]) / 1000.0
                # Cohort-level dose fits compare mean trajectories, so they use
                # the mean initial alive count across wells for initialization.
                N0 = np.nanmean(y_alive_fit[0, :]) if y_alive_fit.ndim > 1 else y_alive_fit[0]
                D0 = np.nanmean(y_dead_fit[0, :]) if y_dead_fit.ndim > 1 else y_dead_fit[0]
                
                dose_data_list.append({
                    'dose_label': gem_dose,
                    'dose_muM': dose_muM,
                    't': t_fit,
                    'y_alive': y_alive_fit,
                    'y_dead': y_dead_fit,
                    'N0': N0,
                    'D0': D0
                })
            except ValueError as e:
                print(f"  Skipping {gem_dose}: {e}")

        if len(dose_data_list) > 0:
            best_cost, best_params, best_n_tr = np.inf, None, 1
            n_range = range(2, 10)
            
            guess_grid = [
                [0.2, 10.0, 0.2],
                [0.5, 25.0, 0.5],
                [1.0, 50.0, 1.0],
                [2.0, 100.0, 1.0],
                [5.0, 50.0, 0.5],
                [10.0, 20.0, 0.2],
            ]
            
            # Tighter k_clear bound than before. Old upper bound of 100 was too loose.
            bounds = ([1e-4, 1e-4, 0.0], [50.0, 500.0, 5.0])
            
            for n_test in n_range:
                for guess in guess_grid:
                    res = least_squares(residuals_global_joint, guess, bounds=bounds,
                                        args=(dose_data_list, r_opt, K_opt, eta_fixed, k_decay_fixed,
                                            k_ext_decay_per_day, n_test, True, 1.0),
                                        loss="soft_l1", f_scale=0.2, max_nfev=3000)
                    
                    if res.cost < best_cost:
                        best_cost, best_params, best_n_tr = res.cost, res.x, n_test

            k_tr_opt, k_kill_opt, k_clear_opt = best_params
            
            print(f"Optimal n_tr: {best_n_tr}")
            print(f"Optimal k_tr: {k_tr_opt:.4f}")
            print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
            print(f"Optimal k_clear: {k_clear_opt:.4f}")
            print(f"Best cost: {best_cost:.4f}")

            plot_global_fit_subplots_joint(
                dose_data_list, ploidy, r_opt, K_opt, best_n_tr, 
                eta_fixed, k_decay_fixed, k_ext_decay_per_day, k_tr_opt, k_kill_opt, k_clear_opt, output_dir=output_dir
            )
        else:
            print(f"No valid data found for {ploidy}. Skipping.")
    replicate_rows = {
        "2N": ["A", "B", "C", "D"],
        "4N": ["E", "F", "G", "H"]
    }
    replicate_fit_results = []

    for ploidy in ploidy_options:
        print(f"\n{'='*60}")
        print(f"ROW-REPLICATE GLOBAL JOINT FIT | COHORT: {ploidy}")
        print(f"{'='*60}")
        
        if ploidy == "2N":
            r_opt, K_opt = params_2N_baseline["r"], params_2N_baseline["K"]
            eta_fixed = eta_2N
            k_decay_fixed = k_decay_2N
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
            eta_fixed = eta_4N
            k_decay_fixed = k_decay_4N
        else:
            continue

        for plate_row in replicate_rows[ploidy]:
            print(f"\nGathering data for {ploidy}, row replicate {plate_row}...")
            
            dose_data_list = build_row_replicate_dose_data_list(
                df=df,
                gem_doses=gem_doses,
                ploidy=ploidy,
                plate_row=plate_row,
                t_max=t_max
            )

            if len(dose_data_list) < 3:
                print(f"  Skipping {ploidy} row {plate_row}: not enough dose curves")
                continue

            result = fit_joint_one_replicate(
                dose_data_list=dose_data_list,
                r_opt=r_opt,
                K_opt=K_opt,
                eta_fixed=eta_fixed,
                k_decay_fixed=k_decay_fixed,
                k_ext_decay_per_day=k_ext_decay_per_day,
                ploidy=ploidy,
                rep_idx=plate_row
            )

            if result is not None:
                result['plate_row'] = plate_row
                replicate_fit_results.append(result)

                # plot_global_fit_subplots_joint(
                #     dose_data_list, 
                #     f"{ploidy} Row {plate_row}", 
                #     r_opt, 
                #     K_opt, 
                #     result['n_tr'], 
                #     eta_fixed, 
                #     k_decay_fixed, 
                #     result['k_tr'], 
                #     result['k_kill'], 
                #     result['k_clear']
                # )

    replicate_fit_df = pd.DataFrame(replicate_fit_results)

    print("\nRow-replicate fit summary:")
    print(replicate_fit_df)
    replicate_fit_df.to_csv(output_dir / "row_replicate_fit_summary.tsv", sep="\t", index=False)

    plot_replicate_parameter_summary(replicate_fit_df, output_dir=output_dir)
    plt.show()    
        
