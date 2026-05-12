"""
Retained for diagnostics/backward comparison only; not used by the preferred
live/dead model.

The preferred live/dead model in `invitro_fitting.py` is driven directly by
intracellular dFdCTP signal. Parent Gemcitabine PK is not used to generate the
live/dead drug signal. This module preserves the older parent-Gemcitabine
exposure and eta/k_decay utilities for legacy plots, comparisons, and tests.
"""

import re
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import least_squares
from scipy.stats import linregress

# Legacy scale constant retained only for explicit fallback/diagnostic paths.
MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM = 0.00971

# The PK workbook reports "Gemcitabine (ng/mL)". The repo does not explicitly
# state whether the assay is reported as free-base or hydrochloride equivalent.
GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL = 263.2

# PK dFdCTP measurements are reported in ng/mL in the workbook. Dividing by the
# molecular weight converts ng/mL to nmol/mL, which is numerically equal to uM.
DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL = 503.1


def decay_rate_from_half_life_days(half_life_days: float) -> float:
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days


def slugify_label(value: str) -> str:
    return re.sub(r'[^A-Za-z0-9]+', '_', value.strip()).strip('_').lower()


def identify_pk_analyte_columns(df: pd.DataFrame) -> Dict[str, str]:
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


@dataclass
class ExtracellularExposureCurve:
    """Legacy parent-Gemcitabine concentration curve in uM."""
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

        if self.mode in {"half_life_fallback", "fitted_exponential", "constrained_bolus_exponential", "legacy_half_life_fallback"}:
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
        else:
            raise ValueError(f"Unsupported exposure curve mode: {self.mode}")

        dose_scale = dose_muM / self.reference_dose_muM if self.reference_dose_muM > 0 else 0.0
        scaled_values = ref_values * dose_scale
        if scalar_input:
            return float(scaled_values[0])
        return scaled_values


def gemcitabine_ng_per_ml_to_uM(
    gemcitabine_ng_per_ml: Union[float, np.ndarray],
    molecular_weight_ng_per_nmol: float = GEMCITABINE_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Union[float, np.ndarray]:
    if molecular_weight_ng_per_nmol <= 0:
        raise ValueError(
            "molecular_weight_ng_per_nmol must be > 0, "
            f"got {molecular_weight_ng_per_nmol}"
        )
    values = np.asarray(gemcitabine_ng_per_ml, dtype=float) / molecular_weight_ng_per_nmol
    if np.ndim(gemcitabine_ng_per_ml) == 0:
        return float(values)
    return values


def converted_extracellular_gem_signal(
    dose_muM: float,
    conversion_factor: float = MODEL_GEM_EXPOSURE_SIGNAL_PER_DOSE_UM,
) -> float:
    if dose_muM <= 0:
        raise ValueError(f"dose_muM must be > 0, got {dose_muM}")
    if conversion_factor <= 0:
        raise ValueError(f"conversion_factor must be > 0, got {conversion_factor}")
    return conversion_factor * dose_muM


def extracellular_gemcitabine_concentration(t, dose_muM, k_ext_decay_per_day):
    initial_exposure_uM_equivalent = converted_extracellular_gem_signal(dose_muM)
    return initial_exposure_uM_equivalent * np.exp(-k_ext_decay_per_day * t)


def fit_constrained_extracellular_gem_decay(
    time_days,
    observed_gem_uM,
    censored_mask,
    censor_upper_bound_uM,
    reference_dose_muM: float,
) -> Optional[Dict[str, float]]:
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
    output_dir=None,
    close_fig: bool = True,
):
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
    if close_fig:
        plt.close(fig)


def dfdctp_ng_per_ml_to_uM(
    dfdctp_ng_per_ml: float,
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> float:
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
    if dose_muM <= 0:
        raise ValueError(f"dose_muM must be > 0, got {dose_muM}")
    return (delta_dfdctp_uM_per_hour / dose_muM) * 24.0


def extract_eta(
    df: pd.DataFrame,
    dose_muM: float = 1.0,
    analyte: str = 'dFdCTP (ng/mL)',
    dfdctp_molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[float]:
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


def simulate_intracellular_dfdctp_signal(
    time_days,
    exposure_curve: ExtracellularExposureCurve,
    dose_muM: float,
    eta_per_day: float,
    k_decay_per_day: float,
    initial_signal: float,
) -> np.ndarray:
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
    output_dir=None,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
    close_fig: bool = True,
):
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
    ax.plot(t_grid, curve_grid, color="darkgreen", linewidth=2.2, label="Intracellular PK Model")
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
    if close_fig:
        plt.close(fig)


def fit_intracellular_pk_parameters(
    df: pd.DataFrame,
    exposure_curve: ExtracellularExposureCurve,
    reference_dose_muM: float,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[Dict[str, float]]:
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


def extract_mean_dfdctp_signal_profile(
    df: pd.DataFrame,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Tuple[np.ndarray, np.ndarray]:
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
