
import math
import re
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union
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

# Legacy scale constant retained only for explicit fallback/diagnostic paths.
#
# The local drugKinetics README documents the heuristic relation
# `C_peak = 0.00971 * drug_intake` for an "intracellular concentration" proxy.
# That heuristic is no longer used in the default empirical extracellular PK
# path, but is preserved for comparison and legacy fallback behavior.
MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM = 0.00971

# The PK workbook reports "Gemcitabine (ng/mL)". The repo does not explicitly
# state whether the assay is reported as free-base or hydrochloride equivalent.
# We use Gemcitabine free-base MW because the analyte is labeled Gemcitabine
# rather than Gemcitabine HCl; this should be revisited if assay metadata say
# otherwise.
GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL = 263.2

# PK dFdCTP measurements are reported in ng/mL in the workbook. Dividing by the
# molecular weight converts ng/mL to nmol/mL, which is numerically equal to uM.
DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL = 503.1

# Current project guidance is that the main 2N/4N PK sheets correspond to
# 1.0 uM Gemcitabine reference dosing, while the explicitly named
# low-initial-Gemcitabine sheets correspond to 0.1 uM.
PKPD_REFERENCE_DOSE_BY_SHEET_UM = {
    "2N": 1.0,
    "4N": 1.0,
    "2N_lowInitialGemcitabine": 0.1,
    "4N_lowInitialGemcitabine": 0.1,
}

for name, path in PATHS.items():
    if not path.exists():
        print(f"Warning: {name} not found at {path}")


def slugify_label(value: str) -> str:
    """Converts labels to filesystem-safe filename fragments."""
    return re.sub(r'[^A-Za-z0-9]+', '_', value.strip()).strip('_').lower()


def get_pk_reference_dose_uM(sheet_name: str) -> float:
    """Returns the documented administered Gemcitabine dose for a PK workbook sheet."""
    if sheet_name not in PKPD_REFERENCE_DOSE_BY_SHEET_UM:
        raise KeyError(
            f"No documented PK reference dose configured for sheet '{sheet_name}'. "
            "Update PKPD_REFERENCE_DOSE_BY_SHEET_UM to match the workbook metadata."
        )
    return PKPD_REFERENCE_DOSE_BY_SHEET_UM[sheet_name]


@dataclass
class ExtracellularExposureCurve:
    """
    Dose-scaled extracellular Gemcitabine concentration curve in uM.

    For empirical mode the in-range segment is interpolated from measured
    Gemcitabine PK concentrations after conversion from ng/mL to uM, while times
    after the last measurement use the fitted exponential tail in the same uM
    concentration units.
    """
    mode: str
    reference_dose_muM: float
    analyte_column: Optional[str]
    observed_time_days: Optional[np.ndarray]
    observed_concentration_uM_values: Optional[np.ndarray]
    censored_time_days: Optional[np.ndarray]
    censored_upper_bound_uM_values: Optional[np.ndarray]
    reference_time_days: Optional[np.ndarray]
    reference_concentration_uM_values: Optional[np.ndarray]
    tail_c0_refdose_uM: float
    k_ext_decay_per_day: float
    half_life_days: float
    fit_r2: Optional[float]
    fit_rmse: Optional[float] = None
    n_uncensored_observations: int = 0
    n_censored_observations: int = 0
    n_numeric_censored_observations: int = 0
    n_censoring_violations: int = 0
    source_ploidy: Optional[str] = None

    def __call__(self, t_days, dose_muM):
        if dose_muM < 0:
            raise ValueError(f"dose_muM must be >= 0, got {dose_muM}")

        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)

        if self.mode == "half_life_fallback":
            ref_values = self.tail_c0_refdose_uM * np.exp(-self.k_ext_decay_per_day * t_eval)
        elif self.mode == "fitted_exponential":
            ref_values = self.tail_c0_refdose_uM * np.exp(-self.k_ext_decay_per_day * t_eval)
        elif self.mode == "constrained_bolus_exponential":
            ref_values = self.tail_c0_refdose_uM * np.exp(-self.k_ext_decay_per_day * t_eval)
        elif self.mode == "empirical + exponential tail":
            assert self.reference_time_days is not None
            assert self.reference_concentration_uM_values is not None
            ref_values = np.interp(
                t_eval,
                self.reference_time_days,
                self.reference_concentration_uM_values,
                left=self.reference_concentration_uM_values[0],
                right=np.nan,
            )
            tail_mask = t_eval > self.reference_time_days[-1]
            if np.any(tail_mask):
                ref_values[tail_mask] = self.tail_c0_refdose_uM * np.exp(-self.k_ext_decay_per_day * t_eval[tail_mask])
        elif self.mode == "legacy_half_life_fallback":
            ref_values = self.tail_c0_refdose_uM * np.exp(-self.k_ext_decay_per_day * t_eval)
        else:
            raise ValueError(f"Unsupported exposure curve mode: {self.mode}")

        dose_scale = dose_muM / self.reference_dose_muM if self.reference_dose_muM > 0 else 0.0
        scaled_values = ref_values * dose_scale
        if scalar_input:
            return float(scaled_values[0])
        return scaled_values


#################### eta, K_decay fitting utilities ####################

def parse_pk_concentration_value(
    value,
    censored_strategy: str = "nan",
) -> Tuple[float, bool, float]:
    """
    Parses PK concentration values while preserving censoring semantics.

    Returns `(numeric_value, was_censored, censor_upper_bound)`. Missing spreadsheet values
    (`np.nan`, blank cells) return `(np.nan, False)`. Censored strings such as
    `BDL`, `BQL`, `BLQ`, `N/F`, and `<0.5` return `was_censored=True`.

    For strings like `<0.5`, `censor_upper_bound` preserves the numeric upper
    bound in the original assay units. Token-only censored values retain
    `censor_upper_bound=np.nan` unless assay metadata provide a numeric LLOQ.
    """
    if censored_strategy not in {"nan", "zero", "half_lod"}:
        raise ValueError(
            "censored_strategy must be one of {'nan', 'zero', 'half_lod'}, "
            f"got {censored_strategy}"
        )

    if pd.isna(value):
        return np.nan, False, np.nan

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), False, np.nan

    value_str = str(value).strip()
    if value_str == "":
        return np.nan, False, np.nan

    try:
        return float(value_str), False, np.nan
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
        return np.nan, False, np.nan

    threshold_match = re.search(r'<\s*([0-9]+(?:\.[0-9]+)?)', value_str)
    if threshold_match:
        threshold = float(threshold_match.group(1))
        if censored_strategy == "half_lod":
            return threshold / 2.0, True, threshold
        if censored_strategy == "zero":
            return 0.0, True, threshold
        return np.nan, True, threshold

    if normalized in censored_tokens:
        if censored_strategy == "zero":
            return 0.0, True, np.nan
        return np.nan, True, np.nan

    return np.nan, False, np.nan

def import_and_clean_pkpd(censored_strategy: str = "nan"):
    """
    Imports PK sheets and parses analyte columns with explicit censored-value handling.

    `censored_strategy="nan"` is the conservative default: BDL/BQL/BLQ-style
    entries become `np.nan` and are ignored by downstream means/fits. Set
    `"zero"` to reproduce the prior zero-imputation behavior explicitly, or
    `"half_lod"` to impute half the reported threshold for `<x` values. The
    original censoring metadata and any numeric censor upper bounds are retained
    in companion columns for downstream PK fitting.
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
            df[f"{col}__censor_upper_bound"] = parsed_values.apply(lambda item: item[2]).astype(float)
            
        clean_data_dict[sheet_name] = df
    return clean_data_dict

def identify_pk_analyte_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Finds key PK analyte columns by normalized name."""
    analyte_map: Dict[str, str] = {}
    for col in df.columns:
        normalized = str(col).strip().lower()
        if '(ng/ml)' not in normalized:
            continue
        if 'gemcitabine' in normalized and 'gemcitabine' not in analyte_map:
            analyte_map['gemcitabine'] = col
        elif 'dfdu' in normalized and 'dfdu' not in analyte_map:
            analyte_map['dfdu'] = col
        elif 'dfdctp' in normalized and 'dfdctp' not in analyte_map:
            analyte_map['dfdctp'] = col
        elif 'dfdcdp' in normalized and 'dfdcdp' not in analyte_map:
            analyte_map['dfdcdp'] = col
        elif 'dfdcmp' in normalized and 'dfdcmp' not in analyte_map:
            analyte_map['dfdcmp'] = col
    return analyte_map

def fit_exponential_pk_decay_model(time_days, concentrations) -> Optional[Dict[str, float]]:
    """
    Fits C(t) = C0 * exp(-k * t) on positive, non-missing concentrations.

    The fit uses the decay phase beginning at the observed peak concentration.
    """
    time_arr = np.asarray(time_days, dtype=float)
    conc_arr = np.asarray(concentrations, dtype=float)
    valid_mask = np.isfinite(time_arr) & np.isfinite(conc_arr) & (conc_arr > 0)
    if np.count_nonzero(valid_mask) < 2:
        return None

    time_valid = time_arr[valid_mask]
    conc_valid = conc_arr[valid_mask]
    peak_idx = int(np.argmax(conc_valid))
    tail_time = time_valid[peak_idx:]
    tail_conc = conc_valid[peak_idx:]
    if len(tail_time) < 2:
        return None

    slope, intercept, r_val, _, _ = linregress(tail_time, np.log(tail_conc))
    k_ext_decay_per_day = -slope
    if k_ext_decay_per_day <= 0:
        return None

    c0 = float(np.exp(intercept))
    return {
        "C0": c0,
        "k_ext_decay_per_day": float(k_ext_decay_per_day),
        "half_life_days": float(np.log(2.0) / k_ext_decay_per_day),
        "r2": float(r_val ** 2),
    }

def gemcitabine_ng_per_ml_to_uM(
    gemcitabine_ng_per_ml: Union[float, np.ndarray],
    molecular_weight_ng_per_nmol: float = GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Union[float, np.ndarray]:
    """Converts extracellular Gemcitabine concentration from ng/mL to uM."""
    if molecular_weight_ng_per_nmol <= 0:
        raise ValueError(
            "molecular_weight_ng_per_nmol must be > 0, "
            f"got {molecular_weight_ng_per_nmol}"
        )
    values = np.asarray(gemcitabine_ng_per_ml, dtype=float) / molecular_weight_ng_per_nmol
    if np.ndim(gemcitabine_ng_per_ml) == 0:
        return float(values)
    return values

def fit_constrained_extracellular_gem_decay(
    time_days,
    observed_gem_uM,
    censored_mask,
    censor_upper_bound_uM,
    reference_dose_muM: float,
) -> Optional[Dict[str, float]]:
    """
    Fits C_ext(t) = C0 * exp(-k_ext * t) with C0 fixed to the administered dose.

    Uncensored positive points contribute log-space residuals. Numeric censored
    upper bounds contribute one-sided penalties only when the model exceeds the
    reported upper bound.
    """
    if reference_dose_muM <= 0:
        raise ValueError(f"reference_dose_muM must be > 0, got {reference_dose_muM}")

    time_arr = np.asarray(time_days, dtype=float)
    observed_arr = np.asarray(observed_gem_uM, dtype=float)
    censored_arr = np.asarray(censored_mask, dtype=bool)
    upper_arr = np.asarray(censor_upper_bound_uM, dtype=float)

    valid_time_mask = np.isfinite(time_arr) & (time_arr >= 0)
    time_arr = time_arr[valid_time_mask]
    observed_arr = observed_arr[valid_time_mask]
    censored_arr = censored_arr[valid_time_mask]
    upper_arr = upper_arr[valid_time_mask]
    if len(time_arr) == 0:
        return None

    uncensored_mask = (~censored_arr) & np.isfinite(observed_arr) & (observed_arr > 0)
    numeric_censored_mask = censored_arr & np.isfinite(upper_arr) & (upper_arr > 0)
    if np.count_nonzero(uncensored_mask) == 0 and np.count_nonzero(numeric_censored_mask) == 0:
        return None

    k_guess = decay_rate_from_half_life_days(1.0)
    if np.count_nonzero(uncensored_mask) >= 2:
        guess_fit = fit_exponential_pk_decay_model(time_arr[uncensored_mask], observed_arr[uncensored_mask])
        if guess_fit is not None:
            k_guess = guess_fit["k_ext_decay_per_day"]

    def model_concentration(eval_time_days, k_ext_decay_per_day):
        return reference_dose_muM * np.exp(-k_ext_decay_per_day * np.asarray(eval_time_days, dtype=float))

    def residuals(log_k_ext):
        k_ext_decay_per_day = float(np.exp(log_k_ext[0]))
        residual_list: List[np.ndarray] = []
        if np.count_nonzero(uncensored_mask) > 0:
            model_uM = model_concentration(time_arr[uncensored_mask], k_ext_decay_per_day)
            obs_uM = observed_arr[uncensored_mask]
            residual_list.append(np.log(model_uM) - np.log(obs_uM))
        if np.count_nonzero(numeric_censored_mask) > 0:
            model_uM = model_concentration(time_arr[numeric_censored_mask], k_ext_decay_per_day)
            ub_uM = upper_arr[numeric_censored_mask]
            residual_list.append(np.maximum(0.0, np.log(model_uM) - np.log(ub_uM)))
        if not residual_list:
            return np.array([], dtype=float)
        return np.concatenate(residual_list)

    res = least_squares(
        residuals,
        x0=[np.log(k_guess)],
        bounds=([np.log(1e-8)], [np.log(1e4)]),
        loss="linear",
        max_nfev=5000,
    )
    k_ext_decay_per_day = float(np.exp(res.x[0]))
    predicted_uncensored = model_concentration(time_arr[uncensored_mask], k_ext_decay_per_day) if np.count_nonzero(uncensored_mask) > 0 else np.array([], dtype=float)
    predicted_numeric_censored = model_concentration(time_arr[numeric_censored_mask], k_ext_decay_per_day) if np.count_nonzero(numeric_censored_mask) > 0 else np.array([], dtype=float)
    uncensored_fit_r2 = np.nan
    uncensored_fit_rmse = np.nan
    if np.count_nonzero(uncensored_mask) >= 2:
        obs_uncensored = observed_arr[uncensored_mask]
        ss_res = float(np.sum((predicted_uncensored - obs_uncensored) ** 2))
        ss_tot = float(np.sum((obs_uncensored - np.mean(obs_uncensored)) ** 2))
        uncensored_fit_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        uncensored_fit_rmse = float(np.sqrt(np.mean((predicted_uncensored - obs_uncensored) ** 2)))

    n_censoring_violations = 0
    if np.count_nonzero(numeric_censored_mask) > 0:
        n_censoring_violations = int(np.count_nonzero(predicted_numeric_censored > upper_arr[numeric_censored_mask]))

    return {
        "mode": "constrained_bolus_exponential",
        "C0_uM": float(reference_dose_muM),
        "k_ext_decay_per_day": k_ext_decay_per_day,
        "half_life_days": float(np.log(2.0) / k_ext_decay_per_day),
        "r2": float(uncensored_fit_r2) if np.isfinite(uncensored_fit_r2) else np.nan,
        "rmse": float(uncensored_fit_rmse) if np.isfinite(uncensored_fit_rmse) else np.nan,
        "n_uncensored": int(np.count_nonzero(uncensored_mask)),
        "n_censored": int(np.count_nonzero(censored_arr)),
        "n_numeric_censored": int(np.count_nonzero(numeric_censored_mask)),
        "n_censoring_violations": n_censoring_violations,
    }

def build_extracellular_exposure_curve_from_profile(
    time_days,
    gem_concentration_ng_per_ml,
    gem_was_censored,
    gem_censor_upper_bound_ng_per_ml,
    reference_dose_muM: float,
    fallback_half_life_days: float,
    analyte_column: Optional[str],
    source_ploidy: Optional[str],
    preferred_mode: str = "constrained_bolus_exponential",
) -> ExtracellularExposureCurve:
    """
    Builds a dose-scaled extracellular Gemcitabine concentration curve in uM.

    The default path fits a constrained bolus-decay model in which
    `C_ext(0) = reference_dose_muM` and measured/censored Gemcitabine PK values
    constrain `k_ext_decay_per_day`.
    """
    if preferred_mode not in {"constrained_bolus_exponential", "half_life_fallback"}:
        raise ValueError(f"Unsupported preferred_mode: {preferred_mode}")

    fallback_k = decay_rate_from_half_life_days(fallback_half_life_days)

    time_arr = np.asarray(time_days, dtype=float)
    conc_arr = np.asarray(gem_concentration_ng_per_ml, dtype=float)
    censored_mask = np.asarray(gem_was_censored, dtype=bool)
    upper_bound_arr = np.asarray(gem_censor_upper_bound_ng_per_ml, dtype=float)
    valid_time_mask = np.isfinite(time_arr) & (time_arr >= 0)

    if np.count_nonzero(valid_time_mask) == 0 or preferred_mode == "half_life_fallback":
        return ExtracellularExposureCurve(
            mode="legacy_half_life_fallback",
            reference_dose_muM=reference_dose_muM,
            analyte_column=analyte_column,
            observed_time_days=None,
            observed_concentration_uM_values=None,
            censored_time_days=None,
            censored_upper_bound_uM_values=None,
            reference_time_days=None,
            reference_concentration_uM_values=None,
            tail_c0_refdose_uM=converted_extracellular_gem_signal(reference_dose_muM),
            k_ext_decay_per_day=fallback_k,
            half_life_days=fallback_half_life_days,
            fit_r2=None,
            fit_rmse=None,
            n_uncensored_observations=0,
            n_censored_observations=0,
            n_numeric_censored_observations=0,
            n_censoring_violations=0,
            source_ploidy=source_ploidy,
        )

    time_valid = time_arr[valid_time_mask]
    conc_valid = conc_arr[valid_time_mask]
    censored_valid = censored_mask[valid_time_mask]
    upper_valid = upper_bound_arr[valid_time_mask]
    order = np.argsort(time_valid)
    time_valid = time_valid[order]
    conc_valid = conc_valid[order]
    censored_valid = censored_valid[order]
    upper_valid = upper_valid[order]

    observed_conc_uM = gemcitabine_ng_per_ml_to_uM(conc_valid)
    upper_uM = gemcitabine_ng_per_ml_to_uM(upper_valid)
    fit_result = fit_constrained_extracellular_gem_decay(
        time_valid,
        observed_conc_uM,
        censored_valid,
        upper_uM,
        reference_dose_muM=reference_dose_muM,
    )
    if fit_result is not None:
        uncensored_mask = (~censored_valid) & np.isfinite(observed_conc_uM) & (observed_conc_uM > 0)
        numeric_censored_mask = censored_valid & np.isfinite(upper_uM) & (upper_uM > 0)
        return ExtracellularExposureCurve(
            mode=fit_result["mode"],
            reference_dose_muM=reference_dose_muM,
            analyte_column=analyte_column,
            observed_time_days=time_valid[uncensored_mask],
            observed_concentration_uM_values=observed_conc_uM[uncensored_mask],
            censored_time_days=time_valid[numeric_censored_mask],
            censored_upper_bound_uM_values=upper_uM[numeric_censored_mask],
            reference_time_days=None,
            reference_concentration_uM_values=None,
            tail_c0_refdose_uM=fit_result["C0_uM"],
            k_ext_decay_per_day=fit_result["k_ext_decay_per_day"],
            half_life_days=fit_result["half_life_days"],
            fit_r2=fit_result["r2"],
            fit_rmse=fit_result["rmse"],
            n_uncensored_observations=fit_result["n_uncensored"],
            n_censored_observations=fit_result["n_censored"],
            n_numeric_censored_observations=fit_result["n_numeric_censored"],
            n_censoring_violations=fit_result["n_censoring_violations"],
            source_ploidy=source_ploidy,
        )

    return ExtracellularExposureCurve(
        mode="legacy_half_life_fallback",
        reference_dose_muM=reference_dose_muM,
        analyte_column=analyte_column,
        observed_time_days=None,
        observed_concentration_uM_values=None,
        censored_time_days=None,
        censored_upper_bound_uM_values=None,
        reference_time_days=None,
        reference_concentration_uM_values=None,
        tail_c0_refdose_uM=converted_extracellular_gem_signal(reference_dose_muM),
        k_ext_decay_per_day=fallback_k,
        half_life_days=fallback_half_life_days,
        fit_r2=None,
        fit_rmse=None,
        n_uncensored_observations=0,
        n_censored_observations=0,
        n_numeric_censored_observations=0,
        n_censoring_violations=0,
        source_ploidy=source_ploidy,
    )

def build_extracellular_exposure_curve_from_pk_sheet(
    df: pd.DataFrame,
    ploidy_label: str,
    reference_dose_muM: float = 1.0,
    fallback_half_life_days: float = 1.0,
    preferred_mode: str = "constrained_bolus_exponential",
) -> ExtracellularExposureCurve:
    analyte_columns = identify_pk_analyte_columns(df)
    gem_col = analyte_columns.get("gemcitabine")
    if gem_col is None or "Timepoint" not in df.columns:
        return build_extracellular_exposure_curve_from_profile(
            time_days=np.array([], dtype=float),
            gem_concentration_ng_per_ml=np.array([], dtype=float),
            gem_was_censored=np.array([], dtype=bool),
            gem_censor_upper_bound_ng_per_ml=np.array([], dtype=float),
            reference_dose_muM=reference_dose_muM,
            fallback_half_life_days=fallback_half_life_days,
            analyte_column=gem_col,
            source_ploidy=ploidy_label,
            preferred_mode="half_life_fallback",
        )
    censor_col = f"{gem_col}__was_censored"
    upper_col = f"{gem_col}__censor_upper_bound"
    return build_extracellular_exposure_curve_from_profile(
        time_days=df["Timepoint"].to_numpy(dtype=float) / 24.0,
        gem_concentration_ng_per_ml=df[gem_col].to_numpy(dtype=float),
        gem_was_censored=df[censor_col].to_numpy(dtype=bool) if censor_col in df.columns else np.zeros(len(df), dtype=bool),
        gem_censor_upper_bound_ng_per_ml=df[upper_col].to_numpy(dtype=float) if upper_col in df.columns else np.full(len(df), np.nan, dtype=float),
        reference_dose_muM=reference_dose_muM,
        fallback_half_life_days=fallback_half_life_days,
        analyte_column=gem_col,
        source_ploidy=ploidy_label,
        preferred_mode=preferred_mode,
    )

def print_pk_workbook_summary(pk_sheets: Dict[str, pd.DataFrame]) -> None:
    print("PK workbook:")
    for sheet_name, df in pk_sheets.items():
        analyte_columns = identify_pk_analyte_columns(df)
        analyte_names = []
        for key in ["gemcitabine", "dfdu", "dfdctp", "dfdcdp", "dfdcmp"]:
            if key in analyte_columns:
                analyte_names.append(analyte_columns[key].replace(" (ng/mL)", ""))
        analyte_summary = ", ".join(analyte_names) if analyte_names else "no recognized analytes"
        print(f"  {sheet_name} sheet: found {analyte_summary}")

def print_exposure_curve_summary(ploidy_label: str, curve: ExtracellularExposureCurve) -> None:
    c0 = curve(0.0, curve.reference_dose_muM)
    c1h = curve(1.0 / 24.0, curve.reference_dose_muM)
    c2h = curve(2.0 / 24.0, curve.reference_dose_muM)
    c24h = curve(1.0, curve.reference_dose_muM)
    c48h = curve(2.0, curve.reference_dose_muM)
    c5 = curve(5.0, curve.reference_dose_muM)
    first_positive_time = np.nan
    if curve.observed_time_days is not None and len(curve.observed_time_days) > 0:
        first_positive_time = float(np.min(curve.observed_time_days))
    print(f"{ploidy_label} extracellular exposure:")
    print(f"  mode: {curve.mode}")
    print(f"  analyte: {curve.analyte_column or 'fallback'}")
    print(f"  reference dose: {curve.reference_dose_muM:.3f} uM")
    if curve.fit_r2 is not None:
        print(f"  R^2: {curve.fit_r2:.4f}")
    if curve.fit_rmse is not None and np.isfinite(curve.fit_rmse):
        print(f"  RMSE: {curve.fit_rmse:.4f} uM")
    print(f"  uncensored observations: {curve.n_uncensored_observations}")
    print(f"  censored observations: {curve.n_censored_observations}")
    print(f"  numeric censored upper bounds: {curve.n_numeric_censored_observations}")
    print(
        "  first positive uncensored time: "
        f"{first_positive_time:.4f} days" if np.isfinite(first_positive_time) else
        "  first positive uncensored time: none"
    )
    print(f"  censoring violations: {curve.n_censoring_violations}")
    print(f"  half-life: {curve.half_life_days:.4f} days")
    print(f"  C_ext(0d): {c0:.4f} uM")
    print(f"  C_ext(1h): {c1h:.4f} uM")
    print(f"  C_ext(2h): {c2h:.4f} uM")
    print(f"  C_ext(24h): {c24h:.4f} uM")
    print(f"  C_ext(48h): {c48h:.4f} uM")
    print(f"  C_ext(5d): {c5:.4f} uM")

def print_extracellular_pk_scale_diagnostic(
    cohort_label: str,
    pk_df: pd.DataFrame,
    curve: ExtracellularExposureCurve,
) -> None:
    """Prints the empirical-vs-legacy scale comparison requested for PK review."""
    if curve.analyte_column is None or "Timepoint" not in pk_df.columns:
        print(f"{cohort_label} extracellular PK scale diagnostic: no measured Gemcitabine data")
        return

    gem_mean = pk_df.groupby("Timepoint")[curve.analyte_column].mean().dropna()
    gem_mean = gem_mean[gem_mean > 0]
    if gem_mean.empty:
        print(f"{cohort_label} extracellular PK scale diagnostic: no positive measured Gemcitabine data")
        return

    first_ng_per_ml = float(gem_mean.iloc[0])
    first_uM = float(gemcitabine_ng_per_ml_to_uM(first_ng_per_ml))
    legacy_value = float(converted_extracellular_gem_signal(curve.reference_dose_muM))
    legacy_to_empirical_ratio = legacy_value / first_uM if first_uM > 0 else np.nan

    print(f"{cohort_label} extracellular PK scale diagnostic:")
    print(f"  reference PK dose: {curve.reference_dose_muM:.3f} uM")
    print(f"  first observed Gemcitabine: {first_ng_per_ml:.4f} ng/mL")
    print(f"  first observed Gemcitabine: {first_uM:.6f} uM")
    print(f"  legacy heuristic only: {legacy_value:.6f} uM-equivalent")
    print(f"  legacy / empirical ratio: {legacy_to_empirical_ratio:.4f}")

def plot_extracellular_exposure_curve(
    ploidy_label: str,
    curve: ExtracellularExposureCurve,
    t_max_days: float = 5.0,
    output_dir: Optional[Path] = None,
):
    """Plots measured PK-derived extracellular Gemcitabine concentration against the model exposure curve."""
    t_grid = np.linspace(0.0, t_max_days, 300)
    curve_values = curve(t_grid, curve.reference_dose_muM)
    fit_quality_parts = [
        f"mode: {curve.mode}",
        f"C0 fixed: {curve.tail_c0_refdose_uM:.3f} uM",
        f"half-life: {curve.half_life_days:.3f} d",
    ]
    if curve.fit_r2 is not None:
        fit_quality_parts.append(f"R^2: {curve.fit_r2:.3f}")
    else:
        fit_quality_parts.append("R^2: n/a")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(
        t_grid,
        curve_values,
        color="darkblue",
        linewidth=2.2,
        label=f"Exposure Curve ({curve.mode})",
    )

    if curve.observed_time_days is not None and curve.observed_concentration_uM_values is not None:
        ax.scatter(
            curve.observed_time_days,
            curve.observed_concentration_uM_values,
            color="orange",
            edgecolor="black",
            linewidth=0.5,
            s=55,
            zorder=3,
            label="Uncensored Gemcitabine PK (uM)",
        )
    if curve.censored_time_days is not None and curve.censored_upper_bound_uM_values is not None and len(curve.censored_time_days) > 0:
        ax.scatter(
            curve.censored_time_days,
            curve.censored_upper_bound_uM_values,
            marker="v",
            facecolors="none",
            edgecolors="crimson",
            linewidth=1.0,
            s=70,
            zorder=3,
            label="Censored upper bounds (uM)",
        )
    ax.axhline(curve.tail_c0_refdose_uM, color="gray", linestyle=":", linewidth=1.0, alpha=0.6)
    ax.text(0.98, 0.02, "C0 tied to administered bolus dose", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Extracellular Gemcitabine Concentration (uM)")
    ax.set_title(
        f"{ploidy_label} Extracellular Gemcitabine PK\n"
        f"{curve.analyte_column or 'Fallback'} | ref dose {curve.reference_dose_muM:.3f} uM"
    )
    ax.text(
        0.02,
        0.98,
        "\n".join(fit_quality_parts),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()

    if output_dir is not None:
        output_path = output_dir / f"extracellular_pk_curve_{slugify_label(ploidy_label)}.png"
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

def extract_mean_dfdctp_signal_profile(
    df: pd.DataFrame,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Tuple[np.ndarray, np.ndarray]:
    """Returns mean dFdCTP measurements converted from ng/mL to intracellular uM."""
    if "Timepoint" not in df.columns or analyte not in df.columns:
        return np.array([], dtype=float), np.array([], dtype=float)

    mean_data = df.groupby("Timepoint")[analyte].mean().dropna()
    if mean_data.empty:
        return np.array([], dtype=float), np.array([], dtype=float)

    time_days = mean_data.index.to_numpy(dtype=float) / 24.0
    concentration_uM_values = dfdctp_ng_per_ml_to_uM(
        mean_data.to_numpy(dtype=float),
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    return time_days, np.asarray(concentration_uM_values, dtype=float)

def simulate_intracellular_dfdctp_signal(
    time_days,
    exposure_curve: ExtracellularExposureCurve,
    dose_muM: float,
    eta_per_day: float,
    k_decay_per_day: float,
    initial_signal: float,
) -> np.ndarray:
    """Simulates intracellular dFdCTP concentration in uM using the PK-derived exposure curve."""
    t_eval = np.asarray(time_days, dtype=float)
    if t_eval.ndim != 1:
        raise ValueError("time_days must be a 1D array")
    if len(t_eval) == 0:
        return np.array([], dtype=float)
    if not np.all(np.isfinite(t_eval)):
        raise ValueError("time_days must be finite")
    if not np.all(np.diff(t_eval) >= 0):
        raise ValueError("time_days must be sorted in nondecreasing order")

    def intracellular_signal_rhs(c_dna_uM, t):
        c_ext_uM = float(exposure_curve(t, dose_muM))
        return eta_per_day * c_ext_uM - k_decay_per_day * c_dna_uM

    simulated = odeint(
        intracellular_signal_rhs,
        [float(initial_signal)],
        t_eval,
    ).ravel()
    return simulated

def plot_intracellular_dfdctp_pk_fit(
    ploidy_label: str,
    pk_df: pd.DataFrame,
    exposure_curve: ExtracellularExposureCurve,
    eta_per_day: float,
    k_decay_per_day: float,
    reference_dose_muM: float,
    output_dir: Optional[Path] = None,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
):
    """
    Plots measured dFdCTP PK against the intracellular model implied by
    dC_dna/dt = eta * C_ext(t) - k_decay * C_dna.
    """
    time_days, observed_signal = extract_mean_dfdctp_signal_profile(
        pk_df,
        analyte=analyte,
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    if len(time_days) == 0:
        return

    initial_signal = float(observed_signal[0])
    simulated_signal = simulate_intracellular_dfdctp_signal(
        time_days,
        exposure_curve=exposure_curve,
        dose_muM=reference_dose_muM,
        eta_per_day=eta_per_day,
        k_decay_per_day=k_decay_per_day,
        initial_signal=initial_signal,
    )

    residuals = observed_signal - simulated_signal
    valid_mask = np.isfinite(observed_signal) & np.isfinite(simulated_signal)
    if np.count_nonzero(valid_mask) >= 2:
        ss_res = float(np.sum((residuals[valid_mask]) ** 2))
        ss_tot = float(np.sum((observed_signal[valid_mask] - np.mean(observed_signal[valid_mask])) ** 2))
        fit_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rmse = float(np.sqrt(np.mean((residuals[valid_mask]) ** 2)))
    else:
        fit_r2 = np.nan
        rmse = np.nan

    t_grid = np.linspace(float(time_days.min()), float(time_days.max()), 300)
    curve_grid = simulate_intracellular_dfdctp_signal(
        t_grid,
        exposure_curve=exposure_curve,
        dose_muM=reference_dose_muM,
        eta_per_day=eta_per_day,
        k_decay_per_day=k_decay_per_day,
        initial_signal=initial_signal,
    )

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(
        t_grid,
        curve_grid,
        color="darkgreen",
        linewidth=2.2,
        label="Intracellular PK Model",
    )
    ax.scatter(
        time_days,
        observed_signal,
        color="goldenrod",
        edgecolor="black",
        linewidth=0.5,
        s=55,
        zorder=3,
        label="Measured dFdCTP PK (uM)",
    )
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Intracellular dFdCTP Concentration (uM)")
    ax.set_title(
        f"{ploidy_label} Intracellular dFdCTP PK\n"
        f"{analyte} | ref dose {reference_dose_muM:.3f} uM"
    )
    ax.text(
        0.02,
        0.98,
        "\n".join([
            f"eta: {eta_per_day:.3f} day^-1",
            f"k_decay: {k_decay_per_day:.3f} day^-1",
            f"R^2: {fit_r2:.3f}" if np.isfinite(fit_r2) else "R^2: n/a",
            f"RMSE: {rmse:.4f}" if np.isfinite(rmse) else "RMSE: n/a",
        ]),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    plt.tight_layout()

    if output_dir is not None:
        output_path = output_dir / f"intracellular_dfdctp_pk_fit_{slugify_label(ploidy_label)}.png"
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

def fit_intracellular_pk_parameters(
    df: pd.DataFrame,
    exposure_curve: ExtracellularExposureCurve,
    reference_dose_muM: float,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[Dict[str, float]]:
    """
    Jointly fits eta and k_decay to the full mean dFdCTP trajectory for one cohort
    with the extracellular Gemcitabine exposure curve held fixed.
    """
    time_days, observed_signal = extract_mean_dfdctp_signal_profile(
        df,
        analyte=analyte,
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    valid_mask = np.isfinite(time_days) & np.isfinite(observed_signal)
    time_days = time_days[valid_mask]
    observed_signal = observed_signal[valid_mask]
    if len(time_days) < 2:
        warnings.warn(
            f"Insufficient {analyte} observations to jointly fit eta and k_decay.",
            RuntimeWarning,
        )
        return None

    initial_signal = float(observed_signal[0])
    eta_guess = extract_eta(
        df,
        dose_muM=reference_dose_muM,
        analyte=analyte,
        dfdctp_molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    legacy_k_hourly, legacy_half_life_hours, legacy_r2 = extract_half_life(df, analyte=analyte)
    k_decay_guess = legacy_k_hourly * 24.0 if legacy_k_hourly is not None else decay_rate_from_half_life_days(1.0)
    eta_guess = float(eta_guess) if eta_guess is not None and eta_guess > 0 else 1.0
    k_decay_guess = float(k_decay_guess) if np.isfinite(k_decay_guess) and k_decay_guess > 0 else decay_rate_from_half_life_days(1.0)

    def residuals_intracellular_pk(params):
        eta_per_day, k_decay_per_day = params
        simulated_signal = simulate_intracellular_dfdctp_signal(
            time_days,
            exposure_curve=exposure_curve,
            dose_muM=reference_dose_muM,
            eta_per_day=eta_per_day,
            k_decay_per_day=k_decay_per_day,
            initial_signal=initial_signal,
        )
        return observed_signal - simulated_signal

    res = least_squares(
        residuals_intracellular_pk,
        x0=[eta_guess, k_decay_guess],
        bounds=([1e-8, 1e-8], [1e4, 1e3]),
        loss="linear",
        max_nfev=5000,
    )
    eta_opt, k_decay_opt = [float(x) for x in res.x]
    simulated_signal = simulate_intracellular_dfdctp_signal(
        time_days,
        exposure_curve=exposure_curve,
        dose_muM=reference_dose_muM,
        eta_per_day=eta_opt,
        k_decay_per_day=k_decay_opt,
        initial_signal=initial_signal,
    )
    residuals = observed_signal - simulated_signal
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((observed_signal - np.mean(observed_signal)) ** 2))
    fit_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    rmse = float(np.sqrt(np.mean(residuals ** 2)))
    half_life_hours = float(np.log(2.0) / k_decay_opt * 24.0)

    return {
        "eta": eta_opt,
        "k_decay": k_decay_opt,
        "hl_hours": half_life_hours,
        "r2_confidence": fit_r2,
        "rmse": rmse,
        "initial_signal": initial_signal,
        "n_points": int(len(time_days)),
        "legacy_eta_0_to_1h": float(eta_guess),
        "legacy_k_decay_day": float(k_decay_guess),
        "legacy_hl_hours": float(legacy_half_life_hours) if legacy_half_life_hours is not None else np.nan,
        "legacy_r2_confidence": float(legacy_r2) if legacy_r2 is not None else np.nan,
    }

def get_r_eta_parameters(
    sheets,
    exposure_curve_by_ploidy: Dict[str, ExtracellularExposureCurve],
    reference_dose_by_cohort_uM: Dict[str, float],
    dfdctp_molecular_weight_ng_per_nmol=DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
):
    """
    Calculates uptake and decay parameters for each cohort (2N, 4N).

    eta and k_decay are jointly fit to the full mean dFdCTP trajectory for each
    cohort, with the cohort-specific extracellular Gemcitabine exposure curve
    held fixed.
    """
    results = {}
    for cohort in ['4N', '2N']:
        if cohort in sheets and cohort in exposure_curve_by_ploidy:
            fit_result = fit_intracellular_pk_parameters(
                sheets[cohort],
                exposure_curve=exposure_curve_by_ploidy[cohort],
                reference_dose_muM=reference_dose_by_cohort_uM[cohort],
                molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
            )
            if fit_result is not None:
                results[cohort] = fit_result
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
    Legacy helper converting nominal extracellular Gemcitabine dose (uM) to the
    historical heuristic exposure scale.
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
    delta_dfdctp_uM_per_hour: float,
    dose_muM: float,
) -> float:
    """
    Estimates eta from the observed 0h->1h dFdCTP buildup in physical uM units.

    Formula:
      eta [day^-1] = (delta_dFdCTP_uM_per_hour / dose_uM) * 24
    """
    if dose_muM <= 0:
        raise ValueError(f"dose_muM must be > 0, got {dose_muM}")
    return (delta_dfdctp_uM_per_hour / dose_muM) * 24.0

def extract_eta(
    df: pd.DataFrame,
    dose_muM: float = 1.0,
    analyte: str = 'dFdCTP (ng/mL)',
    dfdctp_molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[float]:
    """
    Estimates eta from the 0h->1h increase in PK-measured dFdCTP.

    The workbook reports dFdCTP in ng/mL. This function converts that change to
    intracellular dFdCTP concentration in uM using the dFdCTP molecular weight,
    then divides by the administered extracellular Gemcitabine dose in uM. The
    returned eta therefore has units of day^-1 under the physical uM convention.
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
    delta_c_uM_per_hour = dfdctp_ng_per_ml_to_uM(
        delta_c_hourly,
        molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
    )
    eta_daily = estimate_eta_per_day(
        delta_dfdctp_uM_per_hour=delta_c_uM_per_hour,
        dose_muM=dose_muM,
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
# - `C_ext_uM` is extracellular Gemcitabine concentration in uM
# - `C_dna_uM` is intracellular dFdCTP concentration/effective concentration in uM
# - `eta`, `k_decay`, `k_tr`, and `k_clear` are in day^-1
# - `k_kill` is in day^-1 per uM intracellular dFdCTP, so
#   `kappa = k_kill * delayed_C_dna_uM` is a death hazard in day^-1

def decay_rate_from_half_life_days(half_life_days: float) -> float:
    """Converts an extracellular half-life in days to a decay rate in days^-1."""
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days

def extracellular_gemcitabine_concentration(t, dose_muM, k_ext_decay_per_day):
    """Legacy fallback: exponentially decaying extracellular heuristic exposure scale."""
    initial_exposure_uM_equivalent = converted_extracellular_gem_signal(dose_muM)
    return initial_exposure_uM_equivalent * np.exp(-k_ext_decay_per_day * t)


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

def single_clone_signal_ode_joint(
    y,
    t,
    r,
    K,
    dose_muM,
    eta,
    k_decay,
    exposure_curve: Callable[[Union[float, np.ndarray], float], Union[float, np.ndarray]],
    k_tr,
    k_kill,
    k_clear,
    n_tr,
):
    """
    ODE system tracking:
    - A: Alive cells (y[0])
    - C_dna_uM: intracellular dFdCTP concentration (y[1])
    - Z_1 ... Z_n: Transit compartments (y[2 : 2+n_tr])
    - D_obs: Observable Dead cells (y[-1])
    """
    A = y[0]
    C_dna_uM = y[1]
    D_obs = y[-1]

    # Extracellular gemcitabine is supplied by a precomputed PK exposure curve
    # (empirical, fitted exponential, or half-life fallback) evaluated at model
    # time in days. Intracellular dFdCTP still follows the existing uptake/decay
    # equation below.
    C_ext_uM = float(exposure_curve(t, dose_muM))

    # Active drug kinetics: eta * C_ext_uM and k_decay * C_dna_uM both have
    # units of uM/day, so dC_dna_uM/dt does as well.
    dC_dna_uM = eta * C_ext_uM - k_decay * C_dna_uM
    
    # Transit compartments
    dZ = np.zeros(n_tr)
    if n_tr > 0:
        delayed_c_dna_uM = y[2:2+n_tr]
        dZ[0] = k_tr * (C_dna_uM - delayed_c_dna_uM[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (delayed_c_dna_uM[i-1] - delayed_c_dna_uM[i])
            
        kappa = k_kill * delayed_c_dna_uM[-1]
    else:
        kappa = k_kill * C_dna_uM
        
    # Alive cells: Logistic growth minus drug-induced death
    dA = r * A * (1 - A / K) - kappa * A
    
    # Observable dead cells: Accumulate from death, deplete via lysis/clearance
    dD_obs = kappa * A - k_clear * D_obs
    
    return [dA, dC_dna_uM] + dZ.tolist() + [dD_obs]


def simulate_joint_ext(t, params_3, N0, D0, r, K, dose_muM, eta_fixed, k_decay_fixed, exposure_curve, n_tr):
    """
    Wrapper to simulate the ODE system. 
    params_3 = [k_tr, k_kill, k_clear]
    """
    k_tr, k_kill, k_clear = params_3 
    
    # Initial conditions: [A0, C_dna_0, Z_1_0, ..., Z_n_0, D_obs_0]
    y0 = [N0, 0.0] + [0.0] * n_tr + [D0] 
    
    sol = odeint(single_clone_signal_ode_joint, y0, t, 
                 args=(r, K, dose_muM, eta_fixed, k_decay_fixed, exposure_curve, k_tr, k_kill, k_clear, n_tr))
    
    return sol[:, 0], sol[:, -1]


def residuals_global_joint(params_3, dose_data_list, r_fixed, K_fixed, eta_fixed, k_decay_fixed, exposure_curve, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0):
    """
    Residuals for the joint live/dead fit in raw observed cell-count units.

    `high_dose_weight` remains an optional caller-specified dose weight. It is
    not a normalization term and does not rescale residuals by within-dose
    maxima.
    """
    all_residuals = []
    
    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']
        
        y_alive_pred, y_dead_pred = simulate_joint_ext(
            t_data, params_3, N0, D0, r_fixed, K_fixed, dose_muM, eta_fixed, k_decay_fixed, exposure_curve, n_tr_test
        )
        
        # Option 1: Fit mean trajectories first for stability
        if fit_means_only:
            y_alive_obs = np.nanmean(y_alive_data, axis=1) if y_alive_data.ndim > 1 else y_alive_data
            y_dead_obs = np.nanmean(y_dead_data, axis=1) if y_dead_data.ndim > 1 else y_dead_data
            res_alive = y_alive_obs - y_alive_pred
            res_dead = y_dead_obs - y_dead_pred
        
        # Option 2: Fit all replicate points
        else:
            res_alive = y_alive_data.flatten() - np.repeat(y_alive_pred, y_alive_data.shape[1]) if y_alive_data.ndim > 1 else (y_alive_data - y_alive_pred)
            res_dead = y_dead_data.flatten() - np.repeat(y_dead_pred, y_dead_data.shape[1]) if y_dead_data.ndim > 1 else (y_dead_data - y_dead_pred)
        
        # Optional caller-specified dose weighting.
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
    exposure_curve,
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
                t_data, ode_params_3, N0, D0, r, K, dose_muM, eta, k_decay_fixed, exposure_curve, n_tr
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

def fit_joint_one_replicate(dose_data_list, r_opt, K_opt, eta_fixed, k_decay_fixed, exposure_curve, ploidy, rep_idx):
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
                                      exposure_curve, n_test, False, 1.0),
                                loss="linear", max_nfev=3000)
            
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
    print(f"Best raw-count least-squares cost: {best_cost:.4f}")

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

    skip_row_replicate_fitting = True
    gem_ext_half_life_days = 1.0
    pk_censored_strategy = "nan"
    exposure_curve_mode = "constrained_bolus_exponential"

    dfdctp_molecular_weight_ng_per_nmol = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL
    pk_sheets = import_and_clean_pkpd(censored_strategy=pk_censored_strategy)
    print_pk_workbook_summary(pk_sheets)

    exposure_curve_by_ploidy: Dict[str, ExtracellularExposureCurve] = {}
    reference_pk_dose_by_ploidy_uM: Dict[str, float] = {}
    for ploidy_key in ["2N", "4N"]:
        reference_pk_dose_uM = get_pk_reference_dose_uM(ploidy_key)
        reference_pk_dose_by_ploidy_uM[ploidy_key] = reference_pk_dose_uM
        curve = build_extracellular_exposure_curve_from_pk_sheet(
            pk_sheets.get(ploidy_key, pd.DataFrame()),
            ploidy_label=ploidy_key,
            reference_dose_muM=reference_pk_dose_uM,
            fallback_half_life_days=gem_ext_half_life_days,
            preferred_mode=exposure_curve_mode,
        )
        exposure_curve_by_ploidy[ploidy_key] = curve

    fitted_map = get_r_eta_parameters(
        pk_sheets,
        exposure_curve_by_ploidy=exposure_curve_by_ploidy,
        reference_dose_by_cohort_uM=reference_pk_dose_by_ploidy_uM,
        dfdctp_molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
    )

    try:
        eta_2N, k_decay_2N = fitted_map['2N']['eta'], fitted_map['2N']['k_decay']
        eta_4N, k_decay_4N = fitted_map['4N']['eta'], fitted_map['4N']['k_decay']
        print(f"\n{'='*60}")
        print(f"Successfully extracted PK:")
        print(f"  censored PK strategy: {pk_censored_strategy}")
        print(
            "  documented PK reference doses: "
            f"2N={reference_pk_dose_by_ploidy_uM['2N']:.3f} uM, "
            f"4N={reference_pk_dose_by_ploidy_uM['4N']:.3f} uM"
        )
        print(
            "  2N -> "
            f"eta: {eta_2N:.4f} day^-1, "
            f"k_decay: {k_decay_2N:.4f} day^-1, "
            f"dFdCTP fit R^2: {fitted_map['2N']['r2_confidence']:.4f}, "
            f"RMSE: {fitted_map['2N']['rmse']:.4f}"
        )
        print(
            "  4N -> "
            f"eta: {eta_4N:.4f} day^-1, "
            f"k_decay: {k_decay_4N:.4f} day^-1, "
            f"dFdCTP fit R^2: {fitted_map['4N']['r2_confidence']:.4f}, "
            f"RMSE: {fitted_map['4N']['rmse']:.4f}"
        )
        print("  one-step dFdCTP initializers kept for comparison:")
        print(
            "    2N -> "
            f"eta_0to1h: {fitted_map['2N']['legacy_eta_0_to_1h']:.4f} day^-1, "
            f"k_decay_terminal: {fitted_map['2N']['legacy_k_decay_day']:.4f} day^-1"
        )
        print(
            "    4N -> "
            f"eta_0to1h: {fitted_map['4N']['legacy_eta_0_to_1h']:.4f} day^-1, "
            f"k_decay_terminal: {fitted_map['4N']['legacy_k_decay_day']:.4f} day^-1"
        )
        print(f"{'='*60}")
        print_exposure_curve_summary("2N", exposure_curve_by_ploidy["2N"])
        print_exposure_curve_summary("4N", exposure_curve_by_ploidy["4N"])
        print_extracellular_pk_scale_diagnostic("2N", pk_sheets["2N"], exposure_curve_by_ploidy["2N"])
        print_extracellular_pk_scale_diagnostic("4N", pk_sheets["4N"], exposure_curve_by_ploidy["4N"])
        plot_extracellular_exposure_curve("2N", exposure_curve_by_ploidy["2N"], t_max_days=5.0, output_dir=output_dir)
        plot_extracellular_exposure_curve("4N", exposure_curve_by_ploidy["4N"], t_max_days=5.0, output_dir=output_dir)
        plot_intracellular_dfdctp_pk_fit(
            "2N",
            pk_sheets["2N"],
            exposure_curve_by_ploidy["2N"],
            eta_2N,
            k_decay_2N,
            reference_pk_dose_by_ploidy_uM["2N"],
            output_dir=output_dir,
            molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
        )
        plot_intracellular_dfdctp_pk_fit(
            "4N",
            pk_sheets["4N"],
            exposure_curve_by_ploidy["4N"],
            eta_4N,
            k_decay_4N,
            reference_pk_dose_by_ploidy_uM["4N"],
            output_dir=output_dir,
            molecular_weight_ng_per_nmol=dfdctp_molecular_weight_ng_per_nmol,
        )
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
            exposure_curve = exposure_curve_by_ploidy["2N"]
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
            eta_fixed = eta_4N
            k_decay_fixed = k_decay_4N
            exposure_curve = exposure_curve_by_ploidy["4N"]
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
                                            exposure_curve, n_test, True, 1.0),
                                        loss="linear", max_nfev=3000)
                    
                    if res.cost < best_cost:
                        best_cost, best_params, best_n_tr = res.cost, res.x, n_test

            k_tr_opt, k_kill_opt, k_clear_opt = best_params
            
            print(f"Optimal n_tr: {best_n_tr}")
            print(f"Optimal k_tr: {k_tr_opt:.4f}")
            print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
            print(f"Optimal k_clear: {k_clear_opt:.4f}")
            print(f"Best raw-count least-squares cost: {best_cost:.4f}")

            plot_global_fit_subplots_joint(
                dose_data_list, ploidy, r_opt, K_opt, best_n_tr, 
                eta_fixed, k_decay_fixed, exposure_curve, k_tr_opt, k_kill_opt, k_clear_opt, output_dir=output_dir
            )
        else:
            print(f"No valid data found for {ploidy}. Skipping.")
    if skip_row_replicate_fitting:
        print("\nSkipping row-replicate fitting (skip_row_replicate_fitting=True).")
    else:
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
                exposure_curve = exposure_curve_by_ploidy["2N"]
            elif ploidy == "4N":
                r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
                eta_fixed = eta_4N
                k_decay_fixed = k_decay_4N
                exposure_curve = exposure_curve_by_ploidy["4N"]
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
                    exposure_curve=exposure_curve,
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
        
