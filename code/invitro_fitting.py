
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
from scipy.optimize import least_squares, minimize
from scipy.special import gammaln
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

# Current modeling convention:
# The preferred live/dead model is driven directly by intracellular dFdCTP.
# Parent Gemcitabine PK is not used to generate the live/dead drug signal.

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
class DfdctpDoseProfile:
    """
    One calibrated dose-specific intracellular dFdCTP profile in uM.

    Each PK sheet keeps its own timing, peak value, and terminal tail. The
    higher-level surface interpolates between these dose-specific profiles in
    log-dose space.
    """
    sheet_name: str
    dose_uM: float
    analyte_column: str
    time_days: np.ndarray
    raw_signal_uM_values: np.ndarray
    induced_signal_uM_values: np.ndarray
    peak_signal_uM: float
    peak_time_days: float
    tail_half_life_days: float
    tail_decay_per_day: float
    tail_fit_r2: Optional[float] = None
    warning_message: Optional[str] = None

    def evaluate(self, t_days):
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)
        if len(self.time_days) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        else:
            values = np.interp(
                t_eval,
                self.time_days,
                self.induced_signal_uM_values,
                left=float(self.induced_signal_uM_values[0]),
                right=np.nan,
            )
            tail_mask = t_eval > self.time_days[-1]
            if np.any(tail_mask):
                last_value = float(self.induced_signal_uM_values[-1])
                values[tail_mask] = last_value * np.exp(
                    -self.tail_decay_per_day * (t_eval[tail_mask] - float(self.time_days[-1]))
                )
            values = np.clip(values, 0.0, np.inf)
        if scalar_input:
            return float(values[0])
        return values


@dataclass
class DfdctpSignalSurface:
    """
    Dose-dependent intracellular dFdCTP surface C_dFdCTP(t, dose; ploidy).

    For each calibrated dose, the corresponding PK sheet keeps its own temporal
    profile. Evaluation proceeds by:
    1. evaluating each calibrated dose profile at time t, including its own tail
    2. interpolating the resulting dFdCTP concentrations across administered dose
       in log-dose space

    Outside the calibrated dose range:
    - below the minimum calibrated dose, the minimum-dose profile timing/shape is
      preserved and its amplitude is scaled proportionally with dose
    - above the maximum calibrated dose, the current default is explicit linear
      extrapolation in log-dose space, with warnings reported in summaries and
      plots
    """
    source_ploidy: str
    analyte_column: str
    calibration_profiles_by_dose: Dict[float, DfdctpDoseProfile]
    calibration_sheet_names: Tuple[str, ...]
    extrapolation_policy: str = (
        "within-range log-dose interpolation; below-range proportional scaling "
        "from minimum calibrated profile; above-range log-dose linear extrapolation"
    )
    warning_messages: Tuple[str, ...] = ()

    @property
    def calibration_doses_uM(self) -> np.ndarray:
        return np.array(sorted(self.calibration_profiles_by_dose.keys()), dtype=float)

    @property
    def min_calibration_dose_uM(self) -> float:
        doses = self.calibration_doses_uM
        return float(np.min(doses)) if len(doses) > 0 else np.nan

    @property
    def max_calibration_dose_uM(self) -> float:
        doses = self.calibration_doses_uM
        return float(np.max(doses)) if len(doses) > 0 else np.nan

    def get_matching_calibrated_dose(
        self,
        dose_uM: float,
        rtol: float = 1e-9,
        atol: float = 1e-12,
    ) -> Optional[float]:
        for calibrated_dose in self.calibration_doses_uM:
            if np.isclose(float(dose_uM), float(calibrated_dose), rtol=rtol, atol=atol):
                return float(calibrated_dose)
        return None

    def profile_for_dose(self, dose_uM: float) -> Optional[DfdctpDoseProfile]:
        matched_dose = self.get_matching_calibrated_dose(dose_uM)
        if matched_dose is None:
            return None
        return self.calibration_profiles_by_dose[matched_dose]

    def calibration_status_for_doses(self, dose_uM_values: Sequence[float]) -> str:
        dose_arr = np.asarray(dose_uM_values, dtype=float)
        if len(dose_arr) == 0 or len(self.calibration_profiles_by_dose) == 0:
            return "no live/dead doses provided"
        min_dose = self.min_calibration_dose_uM
        max_dose = self.max_calibration_dose_uM
        below = np.any(dose_arr < min_dose)
        above = np.any(dose_arr > max_dose)
        if below and above:
            return (
                "live/dead doses extend below and above the PK calibration range; "
                "below-range doses use proportional scaling of the minimum-dose profile"
            )
        if below:
            return (
                "live/dead doses extend below the PK calibration range; "
                "below-range doses use proportional scaling of the minimum-dose profile"
            )
        if above:
            return "live/dead doses extend above the PK calibration range"
        return "live/dead doses lie within the PK calibration range"

    def evaluate_at_calibrated_doses(self, t_days) -> Tuple[np.ndarray, np.ndarray]:
        doses = self.calibration_doses_uM
        values = np.vstack(
            [
                np.atleast_1d(
                    np.asarray(self.calibration_profiles_by_dose[float(dose)].evaluate(t_days), dtype=float)
                )
                for dose in doses
            ]
        )
        return doses, values

    def evaluate_below_range_dose(self, t_days, dose_muM):
        if dose_muM < 0:
            raise ValueError(f"dose_muM must be >= 0, got {dose_muM}")
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)
        if dose_muM == 0 or len(self.calibration_profiles_by_dose) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        else:
            min_dose = self.min_calibration_dose_uM
            min_profile = self.calibration_profiles_by_dose[min_dose]
            values = min_profile.evaluate(t_eval) * (float(dose_muM) / float(min_dose))
        if scalar_input:
            return float(np.asarray(values, dtype=float).reshape(-1)[0])
        return np.asarray(values, dtype=float)

    def __call__(self, t_days, dose_muM):
        if dose_muM < 0:
            raise ValueError(f"dose_muM must be >= 0, got {dose_muM}")
        doses = self.calibration_doses_uM
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)

        if len(doses) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        elif dose_muM == 0:
            values = np.zeros_like(t_eval, dtype=float)
        elif len(doses) == 1:
            only_dose = float(doses[-1])
            if dose_muM < only_dose:
                values = self.evaluate_below_range_dose(t_eval, dose_muM)
            else:
                values = self.calibration_profiles_by_dose[only_dose].evaluate(t_eval)
        else:
            matched_dose = self.get_matching_calibrated_dose(dose_muM)
            if matched_dose is not None:
                values = self.calibration_profiles_by_dose[matched_dose].evaluate(t_eval)
            elif dose_muM < self.min_calibration_dose_uM:
                values = self.evaluate_below_range_dose(t_eval, dose_muM)
            else:
                dose_grid, profile_values = self.evaluate_at_calibrated_doses(t_eval)
                log_doses = np.log(dose_grid)
                target_log_dose = np.log(max(float(dose_muM), 1e-12))
                if target_log_dose > log_doses[-1]:
                    idx0, idx1 = len(log_doses) - 2, len(log_doses) - 1
                else:
                    idx1 = int(np.searchsorted(log_doses, target_log_dose, side="right"))
                    idx0 = idx1 - 1
                x0, x1 = float(log_doses[idx0]), float(log_doses[idx1])
                y0 = profile_values[idx0, :]
                y1 = profile_values[idx1, :]
                if x1 == x0:
                    values = y0.copy()
                else:
                    frac = (target_log_dose - x0) / (x1 - x0)
                    values = y0 + frac * (y1 - y0)
                values = np.clip(values, 0.0, np.inf)

        if scalar_input:
            return float(np.asarray(values, dtype=float).reshape(-1)[0])
        return np.asarray(values, dtype=float)


def assert_preferred_live_dead_driver(driver: Any) -> DfdctpSignalSurface:
    """
    Guards the preferred live/dead modeling path.

    The preferred path must be driven directly by a `DfdctpSignalSurface`, not
    by parent Gemcitabine PK or any legacy extracellular-exposure object.
    """
    if not isinstance(driver, DfdctpSignalSurface):
        raise TypeError(
            "Preferred live/dead model driver must be a DfdctpSignalSurface; "
            f"got {type(driver).__name__}"
        )
    return driver


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

def print_dfdctp_signal_curve_summary(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    live_dead_dose_uM_values: Optional[Sequence[float]] = None,
) -> None:
    sheet_dose_summary = ", ".join(
        f"{sheet}={get_pk_reference_dose_uM(sheet):.4f} uM" for sheet in curve.calibration_sheet_names
    )
    print(f"{ploidy_label} intracellular dFdCTP signal:")
    print(f"  analyte driver: {curve.analyte_column}")
    print(f"  PK sheets used: {', '.join(curve.calibration_sheet_names)}")
    print(f"  PK sheet reference doses: {sheet_dose_summary}")
    print(
        "  calibration dose range: "
        f"{curve.min_calibration_dose_uM:.4f} to {curve.max_calibration_dose_uM:.4f} uM"
    )
    print(f"  extrapolation policy: {curve.extrapolation_policy}")
    for dose in curve.calibration_doses_uM:
        profile = curve.calibration_profiles_by_dose[float(dose)]
        print(
            f"  dose {dose:.4f} uM -> "
            f"peak {profile.peak_signal_uM:.6f} uM at {profile.peak_time_days:.4f} d; "
            f"tail t1/2 {profile.tail_half_life_days:.4f} d"
            + (
                f"; tail R^2 {profile.tail_fit_r2:.4f}"
                if profile.tail_fit_r2 is not None else
                "; tail R^2 n/a"
            )
        )
        if profile.warning_message:
            print(f"    warning: {profile.warning_message}")
    for warning_message in curve.warning_messages:
        print(f"  warning: {warning_message}")
    if live_dead_dose_uM_values:
        dose_status = curve.calibration_status_for_doses(live_dead_dose_uM_values)
        print(f"  live/dead dose coverage: {dose_status}")
        dose_arr = np.asarray(live_dead_dose_uM_values, dtype=float)
        if np.any(dose_arr < curve.min_calibration_dose_uM):
            print(
                "  below-range policy: doses below "
                f"{curve.min_calibration_dose_uM:.4f} uM use proportional scaling "
                f"of the {curve.min_calibration_dose_uM:.4f} uM dFdCTP profile "
                f"({np.min(dose_arr):.6f} to {curve.min_calibration_dose_uM:.4f} uM)"
            )
        if np.any(dose_arr > curve.max_calibration_dose_uM):
            print(
                "  above-range policy: doses above the PK calibration range use "
                "log-dose linear extrapolation "
                f"({curve.max_calibration_dose_uM:.4f} to {np.max(dose_arr):.6f} uM)"
            )

def plot_dfdctp_signal_curve(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
):
    """
    Plots measured and modeled baseline-subtracted intracellular dFdCTP signal.
    """
    if len(curve.calibration_profiles_by_dose) == 0:
        return
    max_time = max(float(profile.time_days.max()) for profile in curve.calibration_profiles_by_dose.values() if len(profile.time_days) > 0)
    t_grid = np.linspace(0.0, max(5.0, max_time), 400)

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for idx, dose in enumerate(curve.calibration_doses_uM):
        profile = curve.calibration_profiles_by_dose[float(dose)]
        modeled_signal = curve(t_grid, float(dose))
        ax.plot(
            t_grid,
            modeled_signal,
            linewidth=2.0,
            label=f"Modeled {dose:.3f} uM profile",
        )
        if len(profile.time_days) > 0:
            ax.scatter(
                profile.time_days,
                profile.raw_signal_uM_values,
                s=35,
                alpha=0.55,
                label="Measured dFdCTP PK (raw uM)" if idx == 0 else None,
                color="goldenrod",
            )
            ax.scatter(
                profile.time_days,
                profile.induced_signal_uM_values,
                s=55,
                alpha=0.9,
                label=f"Baseline-subtracted {profile.sheet_name}",
                color="darkorange" if dose == np.max(curve.calibration_doses_uM) else "sandybrown",
                edgecolor="black",
                linewidth=0.4,
            )
    if curve.min_calibration_dose_uM > 0:
        preview_dose = curve.min_calibration_dose_uM / 10.0
        preview_signal = curve(t_grid, preview_dose)
        ax.plot(
            t_grid,
            preview_signal,
            linestyle="--",
            linewidth=1.5,
            color="slateblue",
            alpha=0.9,
            label=f"Below-range policy preview ({preview_dose:.3f} uM)",
        )
    ax.set_xlabel("Time (Days)")
    ax.set_ylabel("Intracellular dFdCTP Signal (uM)")
    ax.set_title(f"{ploidy_label} dFdCTP Signal Driver")
    ax.text(
        0.02,
        0.98,
        "\n".join([
            f"calibration range: {curve.min_calibration_dose_uM:.3f}-{curve.max_calibration_dose_uM:.3f} uM",
            f"policy: {curve.extrapolation_policy}",
        ]),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    if output_dir is not None:
        output_path = output_dir / f"dfdctp_signal_curve_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches="tight")
    elif close_fig:
        plt.close(fig)

def plot_dfdctp_amplitude_scaling(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
):
    """
    Plots calibrated PK-sheet peak dFdCTP signal by administered dose.
    """
    calibration_doses = curve.calibration_doses_uM
    calibration_peaks = np.array(
        [curve.calibration_profiles_by_dose[float(dose)].peak_signal_uM for dose in calibration_doses],
        dtype=float,
    )
    if len(calibration_doses) == 0:
        return

    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.scatter(
        calibration_doses,
        calibration_peaks,
        color="crimson",
        edgecolor="black",
        linewidth=0.4,
        s=60,
        zorder=3,
        label="PK-sheet peak dFdCTP",
    )
    ax.set_xlabel("Administered Gemcitabine Dose (uM)")
    ax.set_ylabel("Peak dFdCTP Signal (uM)")
    ax.set_title(f"{ploidy_label} dFdCTP Calibration Peaks")
    ax.text(
        0.02,
        0.98,
        curve.extrapolation_policy,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()
    if output_dir is not None:
        output_path = output_dir / f"dfdctp_amplitude_scaling_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches="tight")
    elif close_fig:
        plt.close(fig)

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

def baseline_subtract_treatment_induced_signal(signal_uM) -> np.ndarray:
    """
    Baseline-subtracts a dFdCTP signal profile so the curve represents the
    treatment-induced active metabolite above the 0h level.
    """
    signal_arr = np.asarray(signal_uM, dtype=float)
    if signal_arr.ndim != 1:
        raise ValueError("signal_uM must be a 1D array")
    if len(signal_arr) == 0:
        return signal_arr.copy()
    baseline_value = float(signal_arr[0]) if np.isfinite(signal_arr[0]) else 0.0
    return np.clip(signal_arr - baseline_value, 0.0, np.inf)

def build_dfdctp_profile_from_sheet(
    df: pd.DataFrame,
    sheet_name: str,
    reference_dose_uM: float,
    fallback_half_life_days: float = 1.0,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[DfdctpDoseProfile]:
    """
    Builds a mean baseline-subtracted dFdCTP profile from one PK sheet.
    """
    time_days, raw_signal_uM = extract_mean_dfdctp_signal_profile(
        df,
        analyte=analyte,
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    if len(time_days) == 0:
        return None

    induced_signal_uM = baseline_subtract_treatment_induced_signal(raw_signal_uM)
    if not np.any(np.isfinite(induced_signal_uM)):
        return None

    peak_idx = int(np.nanargmax(induced_signal_uM))
    peak_signal_uM = float(induced_signal_uM[peak_idx])
    time_of_peak_days = float(time_days[peak_idx])
    tail_fit = fit_exponential_pk_decay_model(time_days, induced_signal_uM)
    if tail_fit is not None:
        tail_decay_per_day = float(tail_fit["k_ext_decay_per_day"])
        tail_half_life_days = float(tail_fit["half_life_days"])
        fit_r2 = float(tail_fit["r2"])
        warning_message = None
    else:
        tail_half_life_days = float(fallback_half_life_days)
        tail_decay_per_day = float(decay_rate_from_half_life_days(tail_half_life_days))
        fit_r2 = None
        warning_message = (
            f"Using fallback dFdCTP tail half-life of {tail_half_life_days:.3f} days "
            f"for {sheet_name} because the terminal fit was not identifiable."
        )

    return DfdctpDoseProfile(
        sheet_name=sheet_name,
        dose_uM=float(reference_dose_uM),
        analyte_column=analyte,
        time_days=np.asarray(time_days, dtype=float),
        raw_signal_uM_values=np.asarray(raw_signal_uM, dtype=float),
        induced_signal_uM_values=np.asarray(induced_signal_uM, dtype=float),
        peak_signal_uM=peak_signal_uM,
        peak_time_days=time_of_peak_days,
        tail_half_life_days=tail_half_life_days,
        tail_decay_per_day=tail_decay_per_day,
        tail_fit_r2=fit_r2,
        warning_message=warning_message,
    )

def build_dfdctp_signal_curve_for_ploidy(
    pk_sheets: Dict[str, pd.DataFrame],
    ploidy_label: str,
    fallback_half_life_days: float = 1.0,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> DfdctpSignalSurface:
    """
    Builds a ploidy-specific dose-dependent dFdCTP surface from all usable PK
    sheets for that ploidy.
    """
    sheet_candidates = [ploidy_label, f"{ploidy_label}_lowInitialGemcitabine"]
    profiles_by_dose: Dict[float, DfdctpDoseProfile] = {}
    sheet_names_used: List[str] = []
    warning_messages: List[str] = []

    for sheet_name in sheet_candidates:
        sheet_df = pk_sheets.get(sheet_name)
        if sheet_df is None or analyte not in sheet_df.columns:
            continue
        reference_dose = get_pk_reference_dose_uM(sheet_name)
        profile = build_dfdctp_profile_from_sheet(
            sheet_df,
            sheet_name=sheet_name,
            reference_dose_uM=reference_dose,
            fallback_half_life_days=fallback_half_life_days,
            analyte=analyte,
            molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
        )
        if profile is None:
            continue
        profiles_by_dose[reference_dose] = profile
        sheet_names_used.append(sheet_name)
        if profile.warning_message:
            warning_messages.append(profile.warning_message)

    if len(profiles_by_dose) == 0:
        raise ValueError(f"No usable {analyte} PK profile found for ploidy {ploidy_label}")

    return DfdctpSignalSurface(
        source_ploidy=ploidy_label,
        analyte_column=analyte,
        calibration_profiles_by_dose=dict(sorted(profiles_by_dose.items())),
        calibration_sheet_names=tuple(sheet_names_used),
        extrapolation_policy="log-dose linear extrapolation",
        warning_messages=tuple(warning_messages),
    )


def build_preferred_dfdctp_signal_surfaces(
    pk_sheets: Dict[str, pd.DataFrame],
    ploidy_keys: Sequence[str] = ("2N", "4N"),
    fallback_half_life_days: float = 1.0,
) -> Dict[str, DfdctpSignalSurface]:
    """
    Builds the preferred live/dead model drivers from dFdCTP PK only.

    Parent Gemcitabine PK is intentionally not used here. If `Gemcitabine
    (ng/mL)` is absent but dFdCTP is present, the preferred driver still builds.
    """
    return {
        ploidy_key: build_dfdctp_signal_curve_for_ploidy(
            pk_sheets,
            ploidy_label=ploidy_key,
            fallback_half_life_days=fallback_half_life_days,
        )
        for ploidy_key in ploidy_keys
    }

def single_clone_dfdctp_ode_joint(
    y,
    t,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve: Callable[[Union[float, np.ndarray], float], Union[float, np.ndarray]],
    k_tr,
    k_kill,
    k_clear,
    n_tr,
):
    """
    ODE system driven directly by intracellular dFdCTP signal.

    State order:
    - A: Alive cells
    - Z_1 ... Z_n: transit compartments smoothing/delaying dFdCTP signal
    - D_obs: observable dead cells
    """
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    A = y[0]
    D_obs = y[-1]
    c_dfdctp_signal_uM = float(np.asarray(dfdctp_signal_curve(t, dose_muM), dtype=float).reshape(-1)[0])

    dZ = np.zeros(n_tr)
    if n_tr > 0:
        delayed_signal = y[1:1+n_tr]
        dZ[0] = k_tr * (c_dfdctp_signal_uM - delayed_signal[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (delayed_signal[i - 1] - delayed_signal[i])
        kappa = k_kill * delayed_signal[-1]
    else:
        kappa = k_kill * c_dfdctp_signal_uM

    dA = r * A * (1 - A / K) - kappa * A
    dD_obs = kappa * A - k_clear * D_obs
    return [dA] + dZ.tolist() + [dD_obs]

def simulate_joint_dfdctp(t, params_3, N0, D0, r, K, dose_muM, dfdctp_signal_curve, n_tr):
    """
    Simulates the live/dead ODE with dFdCTP supplied directly as the drug-driver
    curve instead of integrating a separate intracellular drug state.
    """
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    k_tr, k_kill, k_clear = params_3
    y0 = [N0] + [0.0] * n_tr + [D0]
    sol = odeint(
        single_clone_dfdctp_ode_joint,
        y0,
        t,
        args=(r, K, dose_muM, dfdctp_signal_curve, k_tr, k_kill, k_clear, n_tr),
    )
    return sol[:, 0], sol[:, -1]
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

def save_and_maybe_close(fig, output_path: Optional[Path] = None, close: bool = True, **savefig_kwargs) -> None:
    """Saves a Matplotlib figure and closes it by default for batch-script safety."""
    if output_path is not None:
        fig.savefig(output_path, **savefig_kwargs)
    if close and (output_path is not None or close):
        plt.close(fig)

def _build_phenotype_pivot(
    df: pd.DataFrame,
    gem_dose: str,
    ploidy: str,
    phenotype: str,
    plate_row: Optional[str] = None,
) -> pd.DataFrame:
    """
    Aggregates counts by time/replicate and returns a pivoted matrix.

    Duplicate entries for the same time/replicate/phenotype are summed because
    the count column is an observed cell count.
    """
    subset = df[
        (df['gem'] == gem_dose) &
        (df['ploidy'] == ploidy) &
        (df['phenotype'] == phenotype)
    ].copy()
    if plate_row is not None:
        subset = subset[subset['plate_row'] == plate_row].copy()
    if subset.empty:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No data found for {gem_dose}, {ploidy}, {phenotype}{scope}")

    grouped = (
        subset
        .groupby(['time_days', 'plate_row', 'plate_col'], as_index=False)['count']
        .sum()
    )
    return grouped.pivot_table(
        index='time_days',
        columns=['plate_row', 'plate_col'],
        values='count',
        aggfunc='sum',
    ).sort_index()

def get_aligned_live_dead_data(
    df: pd.DataFrame,
    gem_dose: str,
    ploidy: str,
    plate_row: Optional[str] = None,
    t_max: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Inner-joins Alive and Dead trajectories by shared time and replicate identity.

    Alignment is performed on `time_days` and the replicate keys
    `plate_row`/`plate_col`. Time points or replicate columns that do not appear
    in both phenotypes are dropped explicitly and reported.
    """
    alive_pivot = _build_phenotype_pivot(df, gem_dose, ploidy, "Alive", plate_row=plate_row)
    dead_pivot = _build_phenotype_pivot(df, gem_dose, ploidy, "Dead", plate_row=plate_row)

    common_times = alive_pivot.index.intersection(dead_pivot.index).sort_values()
    common_cols = alive_pivot.columns.intersection(dead_pivot.columns)
    if len(common_times) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead time points for {gem_dose}, {ploidy}{scope}")
    if len(common_cols) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead replicate columns for {gem_dose}, {ploidy}{scope}")

    alive_time_set = set(alive_pivot.index.tolist())
    dead_time_set = set(dead_pivot.index.tolist())
    dropped_timepoints = len(alive_time_set.union(dead_time_set)) - len(common_times)
    dropped_replicates = len(set(alive_pivot.columns.tolist()).union(set(dead_pivot.columns.tolist()))) - len(common_cols)

    if t_max is not None:
        common_times = common_times[common_times <= t_max]
    if len(common_times) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead time points remain after t_max for {gem_dose}, {ploidy}{scope}")

    alive_aligned = alive_pivot.loc[common_times, common_cols]
    dead_aligned = dead_pivot.loc[common_times, common_cols]

    return {
        "t": common_times.to_numpy(dtype=float),
        "y_alive": alive_aligned.to_numpy(dtype=float),
        "y_dead": dead_aligned.to_numpy(dtype=float),
        "replicate_columns": list(common_cols),
        "dropped_timepoints": int(dropped_timepoints),
        "dropped_replicates": int(dropped_replicates),
        "plate_row": plate_row,
    }

def get_fitting_data(df: pd.DataFrame, gem_dose: str = "0 nM", ploidy: str = "2N", phenotype: str = "Alive"):
    """
    Backward-compatible wrapper returning one phenotype matrix after explicit
    Alive/Dead alignment.
    """
    aligned = get_aligned_live_dead_data(df, gem_dose=gem_dose, ploidy=ploidy)
    if phenotype == "Alive":
        return aligned["t"], aligned["y_alive"]
    if phenotype == "Dead":
        return aligned["t"], aligned["y_dead"]
    raise ValueError(f"Unsupported phenotype {phenotype}; expected Alive or Dead")

def logistic_growth(t, N0, r, K):
    """Analytical solution to the logistic growth equation."""
    return (K * N0) / (N0 + (K - N0) * np.exp(-r * t))

# Canonical internal unit convention for the joint in vitro ODE:
# - time `t` is in days
# - `dose_muM` is nominal extracellular gemcitabine dose in uM
# - `C_dfdctp_uM` is intracellular dFdCTP concentration/effective concentration in uM
# - `k_tr` and `k_clear` are in day^-1
# - `k_kill` is in day^-1 per uM intracellular dFdCTP, so
#   `kappa = k_kill * delayed_C_dfdctp_uM` is a death hazard in day^-1

def decay_rate_from_half_life_days(half_life_days: float) -> float:
    """Converts an extracellular half-life in days to a decay rate in days^-1."""
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days


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

def fit_baseline_shared_rk_replicate_n0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None, close_fig: bool = True):
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
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')
    elif close_fig:
        plt.close(fig)

    return {
        "N0_by_replicate": n0_by_replicate,
        "N0_mean": n0_mean,
        "r": float(r_opt),
        "K": float(K_opt),
        "ploidy": ploidy_label,
        "n_replicates": int(len(n0_by_replicate)),
        "n_observations": n_observations,
    }

def fit_baseline_locked_N0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None, close_fig: bool = True):
    """Backward-compatible wrapper for the replicate-specific-N0 baseline fit."""
    return fit_baseline_shared_rk_replicate_n0(
        t_data,
        y_alive,
        y_dead,
        ploidy_label,
        t_max=t_max,
        output_dir=output_dir,
        close_fig=close_fig,
    )



#################### Cohort fitting utilities ####################

def single_clone_signal_ode_joint(
    y,
    t,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve: Callable[[Union[float, np.ndarray], float], Union[float, np.ndarray]],
    k_tr,
    k_kill,
    k_clear,
    n_tr,
):
    """
    Backward-compatible name for the dFdCTP-driven live/dead ODE.
    """
    return single_clone_dfdctp_ode_joint(
        y=y,
        t=t,
        r=r,
        K=K,
        dose_muM=dose_muM,
        dfdctp_signal_curve=dfdctp_signal_curve,
        k_tr=k_tr,
        k_kill=k_kill,
        k_clear=k_clear,
        n_tr=n_tr,
    )


def simulate_joint_ext(t, params_3, N0, D0, r, K, dose_muM, dfdctp_signal_curve, n_tr):
    """
    Backward-compatible wrapper for the dFdCTP-driven live/dead simulation.
    """
    return simulate_joint_dfdctp(
        t=t,
        params_3=params_3,
        N0=N0,
        D0=D0,
        r=r,
        K=K,
        dose_muM=dose_muM,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr=n_tr,
    )


def poisson_nll(y, mu):
    """Poisson negative log-likelihood for nonnegative count observations."""
    y_arr = np.asarray(y, dtype=float)
    mu_arr = np.clip(np.asarray(mu, dtype=float), 1e-9, np.inf)
    if np.any(y_arr < 0):
        raise ValueError("Poisson observations must be nonnegative")
    return float(np.sum(mu_arr - y_arr * np.log(mu_arr) + gammaln(y_arr + 1.0)))


def negative_binomial_nll(y, mu, theta):
    """
    Negative binomial negative log-likelihood with mean/dispersion parameterization.

    Var(y) = mu + mu^2 / theta
    """
    if theta <= 0:
        raise ValueError(f"theta must be > 0, got {theta}")
    y_arr = np.asarray(y, dtype=float)
    mu_arr = np.clip(np.asarray(mu, dtype=float), 1e-9, np.inf)
    if np.any(y_arr < 0):
        raise ValueError("Negative binomial observations must be nonnegative")

    theta_val = float(theta)
    ll = (
        gammaln(y_arr + theta_val)
        - gammaln(theta_val)
        - gammaln(y_arr + 1.0)
        + theta_val * (np.log(theta_val) - np.log(theta_val + mu_arr))
        + y_arr * (np.log(mu_arr) - np.log(theta_val + mu_arr))
    )
    return float(-np.sum(ll))


def collect_alive_dead_observations(
    dose_data_list,
    params_3,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    fit_means_only=True,
    high_dose_weight=1.0,
):
    """
    Collects observed/predicted alive and dead counts for objective evaluation.

    `high_dose_weight` is an optional caller-specified dose weight. It is not a
    normalization term.
    """
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    records = []

    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']

        y_alive_pred, y_dead_pred = simulate_joint_ext(
            t_data, params_3, N0, D0, r_fixed, K_fixed, dose_muM, dfdctp_signal_curve, n_tr_test
        )

        if fit_means_only:
            y_alive_obs = np.nanmean(y_alive_data, axis=1) if np.ndim(y_alive_data) > 1 else np.asarray(y_alive_data, dtype=float)
            y_dead_obs = np.nanmean(y_dead_data, axis=1) if np.ndim(y_dead_data) > 1 else np.asarray(y_dead_data, dtype=float)
            y_alive_pred_expanded = np.asarray(y_alive_pred, dtype=float)
            y_dead_pred_expanded = np.asarray(y_dead_pred, dtype=float)
        else:
            y_alive_obs = y_alive_data.flatten() if np.ndim(y_alive_data) > 1 else np.asarray(y_alive_data, dtype=float)
            y_dead_obs = y_dead_data.flatten() if np.ndim(y_dead_data) > 1 else np.asarray(y_dead_data, dtype=float)
            if np.ndim(y_alive_data) > 1:
                y_alive_pred_expanded = np.repeat(np.asarray(y_alive_pred, dtype=float), y_alive_data.shape[1])
                y_dead_pred_expanded = np.repeat(np.asarray(y_dead_pred, dtype=float), y_dead_data.shape[1])
            else:
                y_alive_pred_expanded = np.asarray(y_alive_pred, dtype=float)
                y_dead_pred_expanded = np.asarray(y_dead_pred, dtype=float)

        weight = high_dose_weight if dose_muM >= 0.2 else 1.0

        alive_mask = ~np.isnan(y_alive_obs)
        dead_mask = ~np.isnan(y_dead_obs)
        records.append(
            {
                "dose_muM": dose_muM,
                "weight": float(weight),
                "alive_obs": np.asarray(y_alive_obs[alive_mask], dtype=float),
                "alive_pred": np.asarray(y_alive_pred_expanded[alive_mask], dtype=float),
                "dead_obs": np.asarray(y_dead_obs[dead_mask], dtype=float),
                "dead_pred": np.asarray(y_dead_pred_expanded[dead_mask], dtype=float),
            }
        )
    return records


def residuals_global_joint(params_3, dose_data_list, r_fixed, K_fixed, dfdctp_signal_curve, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0):
    """
    Residuals for the joint live/dead fit in raw observed cell-count units.

    `high_dose_weight` remains an optional caller-specified dose weight. It is
    not a normalization term and does not rescale residuals by within-dose
    maxima.
    """
    all_residuals = []

    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        params_3=params_3,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr_test,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        all_residuals.extend((record["alive_obs"] - record["alive_pred"]) * weight)
        all_residuals.extend((record["dead_obs"] - record["dead_pred"]) * weight)

    return np.array(all_residuals)


def live_dead_objective_nll(
    params_3,
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    objective="negative_binomial",
    fit_means_only=True,
    high_dose_weight=1.0,
    theta_alive=None,
    theta_dead=None,
):
    """
    Likelihood objective for alive/dead counts using predicted ODE means.
    """
    if objective not in {"negative_binomial", "poisson"}:
        raise ValueError(f"Unsupported likelihood objective {objective}")

    if objective == "negative_binomial":
        if theta_alive is None or theta_dead is None:
            raise ValueError("theta_alive and theta_dead are required for negative_binomial objective")
        if theta_alive <= 0 or theta_dead <= 0:
            raise ValueError("theta_alive and theta_dead must be > 0")

    total_nll = 0.0
    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        params_3=params_3,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr_test,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        if objective == "poisson":
            total_nll += weight * poisson_nll(record["alive_obs"], record["alive_pred"])
            total_nll += weight * poisson_nll(record["dead_obs"], record["dead_pred"])
        else:
            total_nll += weight * negative_binomial_nll(record["alive_obs"], record["alive_pred"], theta_alive)
            total_nll += weight * negative_binomial_nll(record["dead_obs"], record["dead_pred"], theta_dead)
    return float(total_nll)


def summarize_live_dead_fit(
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr,
    params_3,
    objective,
    fit_means_only=True,
    high_dose_weight=1.0,
    theta_alive=None,
    theta_dead=None,
):
    records = collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        params_3=params_3,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    )
    residual_vector = np.concatenate(
        [
            (record["alive_obs"] - record["alive_pred"]) * record["weight"]
            for record in records
        ] + [
            (record["dead_obs"] - record["dead_pred"]) * record["weight"]
            for record in records
        ]
    ) if records else np.array([], dtype=float)

    n_obs = int(sum(len(record["alive_obs"]) + len(record["dead_obs"]) for record in records))
    rmse = float(np.sqrt(np.mean(residual_vector ** 2))) if residual_vector.size > 0 else np.nan
    summary = {
        "objective": objective,
        "n_observations": n_obs,
        "rmse": rmse,
        "theta_alive": theta_alive,
        "theta_dead": theta_dead,
    }
    if objective == "least_squares":
        summary["objective_value"] = float(0.5 * np.sum(residual_vector ** 2))
        summary["nll"] = None
        summary["aic"] = None
        summary["bic"] = None
        summary["n_parameters"] = 3
    else:
        nll = live_dead_objective_nll(
            params_3=params_3,
            dose_data_list=dose_data_list,
            r_fixed=r_fixed,
            K_fixed=K_fixed,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=n_tr,
            objective=objective,
            fit_means_only=fit_means_only,
            high_dose_weight=high_dose_weight,
            theta_alive=theta_alive,
            theta_dead=theta_dead,
        )
        n_parameters = 5 if objective == "negative_binomial" else 3
        summary["objective_value"] = float(nll)
        summary["nll"] = float(nll)
        summary["n_parameters"] = n_parameters
        summary["aic"] = float(2 * n_parameters + 2 * nll)
        summary["bic"] = float(n_parameters * np.log(max(n_obs, 1)) + 2 * nll)
    return summary


def fit_live_dead_model(
    dose_data_list,
    r_opt,
    K_opt,
    dfdctp_signal_curve,
    n_tr,
    objective="negative_binomial",
    fit_means_only=True,
    high_dose_weight=1.0,
    max_nfev=3000,
):
    """
    Fits live/dead model parameters under the requested objective.

    Available objectives:
    - `negative_binomial` (default)
    - `poisson`
    - `least_squares` (legacy raw-residual objective)
    """
    if objective == "least_squares":
        guess_grid = [
            [0.2, 10.0, 0.2],
            [0.5, 25.0, 0.5],
            [1.0, 50.0, 1.0],
            [2.0, 100.0, 1.0],
            [5.0, 50.0, 0.5],
            [10.0, 20.0, 0.2],
        ]
        bounds = ([1e-4, 1e-4, 0.0], [50.0, 500.0, 5.0])
        best_result = None
        best_summary = None
        best_cost = np.inf
        for guess in guess_grid:
            res = least_squares(
                residuals_global_joint,
                guess,
                bounds=bounds,
                args=(dose_data_list, r_opt, K_opt, dfdctp_signal_curve, n_tr, fit_means_only, high_dose_weight),
                loss="linear",
                max_nfev=max_nfev,
            )
            if res.cost < best_cost:
                best_cost = float(res.cost)
                best_result = res
                best_summary = summarize_live_dead_fit(
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_tr,
                    params_3=res.x,
                    objective=objective,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                )
        if best_result is None or best_summary is None:
            return None
        return {
            "success": bool(best_result.success),
            "params_3": np.asarray(best_result.x, dtype=float),
            "objective": objective,
            "summary": best_summary,
            "optimizer_result": best_result,
        }

    if objective not in {"negative_binomial", "poisson"}:
        raise ValueError(f"Unsupported objective {objective}")

    guess_grid = [
        [0.2, 10.0, 0.2],
        [0.5, 25.0, 0.5],
        [1.0, 50.0, 1.0],
        [2.0, 100.0, 1.0],
        [5.0, 50.0, 0.5],
        [10.0, 20.0, 0.2],
    ]
    theta_guess_grid = [(20.0, 20.0), (50.0, 50.0), (100.0, 30.0)]

    param_lower = np.array([1e-4, 1e-4, 1e-8], dtype=float)
    param_upper = np.array([50.0, 500.0, 5.0], dtype=float)
    theta_lower = np.array([1e-3, 1e-3], dtype=float)
    theta_upper = np.array([1e6, 1e6], dtype=float)

    best_result = None
    best_summary = None
    best_nll = np.inf

    for guess in guess_grid:
        theta_guesses = theta_guess_grid if objective == "negative_binomial" else [(None, None)]
        for theta_guess_alive, theta_guess_dead in theta_guesses:
            if objective == "negative_binomial":
                x0 = np.log(np.array([guess[0], guess[1], max(guess[2], 1e-8), theta_guess_alive, theta_guess_dead], dtype=float))
                lower = np.log(np.concatenate([param_lower, theta_lower]))
                upper = np.log(np.concatenate([param_upper, theta_upper]))
            else:
                x0 = np.log(np.array([guess[0], guess[1], max(guess[2], 1e-8)], dtype=float))
                lower = np.log(param_lower)
                upper = np.log(param_upper)

            def objective_fn(log_params):
                natural_params = np.exp(np.asarray(log_params, dtype=float))
                params_3 = natural_params[:3]
                if objective == "negative_binomial":
                    theta_alive = float(natural_params[3])
                    theta_dead = float(natural_params[4])
                else:
                    theta_alive = None
                    theta_dead = None
                return live_dead_objective_nll(
                    params_3=params_3,
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr_test=n_tr,
                    objective=objective,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                    theta_alive=theta_alive,
                    theta_dead=theta_dead,
                )

            res = minimize(
                objective_fn,
                x0=x0,
                method="L-BFGS-B",
                bounds=list(zip(lower, upper)),
                options={"maxiter": max_nfev},
            )
            if res.fun < best_nll:
                natural_params = np.exp(np.asarray(res.x, dtype=float))
                params_3 = natural_params[:3]
                theta_alive = float(natural_params[3]) if objective == "negative_binomial" else None
                theta_dead = float(natural_params[4]) if objective == "negative_binomial" else None
                best_nll = float(res.fun)
                best_result = res
                best_summary = summarize_live_dead_fit(
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_tr,
                    params_3=params_3,
                    objective=objective,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                    theta_alive=theta_alive,
                    theta_dead=theta_dead,
                )

    if best_result is None or best_summary is None:
        return None

    return {
        "success": bool(best_result.success),
        "params_3": np.asarray(np.exp(np.asarray(best_result.x, dtype=float))[:3], dtype=float),
        "objective": objective,
        "summary": best_summary,
        "optimizer_result": best_result,
    }


def plot_global_fit_subplots_joint(
    dose_data_list,
    ploidy_label,
    r,
    K,
    n_tr,
    dfdctp_signal_curve,
    k_tr,
    k_kill,
    k_clear,
    fit_summary: Optional[Dict[str, Any]] = None,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
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
                t_data, ode_params_3, N0, D0, r, K, dose_muM, dfdctp_signal_curve, n_tr
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
        
    param_lines = [
        f"Optima: n_tr = {n_tr}",
        f"$k_{{tr}}$ = {k_tr:.3f} | $k_{{kill}}$ = {k_kill:.3f} | $k_{{clear}}$ = {k_clear:.3f}",
    ]
    if fit_summary is not None:
        param_lines.append(f"Fit objective: {fit_summary['objective']}")
        if fit_summary.get("theta_alive") is not None and fit_summary.get("theta_dead") is not None:
            param_lines.append(
                f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f} | "
                f"$\\theta_{{dead}}$ = {fit_summary['theta_dead']:.3f}"
            )
        if fit_summary.get("nll") is not None:
            param_lines.append(
                f"NLL = {fit_summary['nll']:.3f} | AIC = {fit_summary['aic']:.3f} | "
                f"BIC = {fit_summary['bic']:.3f}"
            )
        else:
            param_lines.append(f"Legacy raw-count least-squares cost = {fit_summary['objective_value']:.3f}")
        if np.isfinite(fit_summary.get("rmse", np.nan)):
            param_lines.append(f"RMSE = {fit_summary['rmse']:.3f}")
    param_text = "\n".join(param_lines)
    fig.suptitle(f"{ploidy_label} Cohort - Joint Live/Dead Kinetic Fit\n{param_text}", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    if output_dir is not None:
        output_path = output_dir / f"cohort_joint_fit_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')
    elif close_fig:
        plt.close(fig)

def get_fitting_data_one_row(df, gem_dose: str, ploidy: str, phenotype: str, plate_row: str):
    """
    Backward-compatible wrapper returning one phenotype vector after explicit
    Alive/Dead alignment for a single row replicate.

    Example:
        plate_row='A' for 2N replicate A
        plate_row='E' for 4N replicate E
    """
    aligned = get_aligned_live_dead_data(df, gem_dose=gem_dose, ploidy=ploidy, plate_row=plate_row)
    matrix = aligned["y_alive"] if phenotype == "Alive" else aligned["y_dead"] if phenotype == "Dead" else None
    if matrix is None:
        raise ValueError(f"Unsupported phenotype {phenotype}; expected Alive or Dead")
    if matrix.shape[1] != 1:
        raise ValueError(
            f"Expected exactly one aligned replicate column for {gem_dose}, {ploidy}, row {plate_row}; "
            f"got {matrix.shape[1]}"
        )
    return aligned["t"], matrix[:, 0]

def fit_joint_one_replicate(
    dose_data_list,
    r_opt,
    K_opt,
    dfdctp_signal_curve,
    ploidy,
    rep_idx,
    objective="negative_binomial",
):
    """
    Fits one global joint live/dead model across all doses for one replicate.
    """
    best_objective_value, best_fit, best_n_tr = np.inf, None, 1
    n_range = range(2, 10)

    for n_test in n_range:
        fit_result = fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=r_opt,
            K_opt=K_opt,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=n_test,
            objective=objective,
            fit_means_only=False,
            high_dose_weight=1.0,
            max_nfev=3000,
        )
        if fit_result is None:
            continue
        objective_value = float(fit_result["summary"]["objective_value"])
        if objective_value < best_objective_value:
            best_objective_value, best_fit, best_n_tr = objective_value, fit_result, n_test

    if best_fit is None:
        return None

    k_tr_opt, k_kill_opt, k_clear_opt = best_fit["params_3"]
    fit_summary = best_fit["summary"]
    
    print(f"\n--- {ploidy} Replicate {rep_idx} ---")
    print(f"Optimal n_tr: {best_n_tr}")
    print(f"Optimal k_tr: {k_tr_opt:.4f}")
    print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
    print(f"Optimal k_clear: {k_clear_opt:.4f}")
    print(f"Fit objective: {objective}")
    if fit_summary["nll"] is not None:
        print(f"Best negative log-likelihood: {fit_summary['nll']:.4f}")
        print(f"AIC: {fit_summary['aic']:.4f}")
        print(f"BIC: {fit_summary['bic']:.4f}")
        print(f"Alive dispersion theta: {fit_summary['theta_alive']:.4f}")
        print(f"Dead dispersion theta: {fit_summary['theta_dead']:.4f}")
    else:
        print(f"Best raw-count least-squares cost: {fit_summary['objective_value']:.4f}")
    print(f"RMSE: {fit_summary['rmse']:.4f}")

    return {
        'ploidy': ploidy,
        'rep_idx': rep_idx,
        'n_tr': best_n_tr,
        'k_tr': k_tr_opt,
        'k_kill': k_kill_opt,
        'k_clear': k_clear_opt,
        'cost': best_objective_value,
        'objective': objective,
        'theta_alive': fit_summary['theta_alive'],
        'theta_dead': fit_summary['theta_dead'],
        'nll': fit_summary['nll'],
        'aic': fit_summary['aic'],
        'bic': fit_summary['bic'],
        'rmse': fit_summary['rmse'],
    }

def build_row_replicate_dose_data_list(df, gem_doses, ploidy, plate_row, t_max):
    """
    Builds a dose_data_list for one row replicate across all gemcitabine doses.
    """
    dose_data_list = []

    for gem_dose in gem_doses:
        try:
            aligned = get_aligned_live_dead_data(
                df,
                gem_dose=gem_dose,
                ploidy=ploidy,
                plate_row=plate_row,
                t_max=t_max,
            )
            if aligned["dropped_timepoints"] > 0:
                print(
                    f"  {gem_dose}, row {plate_row}: dropped {aligned['dropped_timepoints']} "
                    "non-overlapping Alive/Dead time points"
                )
            if aligned["y_alive"].shape[1] != 1:
                raise ValueError(
                    f"Expected one aligned replicate column for {gem_dose}, {ploidy}, row {plate_row}; "
                    f"got {aligned['y_alive'].shape[1]}"
                )
            t_fit = aligned["t"]
            y_alive_fit = aligned["y_alive"][:, 0]
            y_dead_fit = aligned["y_dead"][:, 0]
            
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

def plot_replicate_parameter_summary(replicate_fit_df, output_dir: Optional[Path] = None, close_fig: bool = True):
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
            save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')
        elif close_fig:
            plt.close(fig)

if __name__ == "__main__":
    output_dir = SCRIPT_PATH.parent / "invitro_fitting_outputs" / datetime.now().strftime("%Y%m%dT%H%M%S")
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving figures to: {output_dir}")

    skip_row_replicate_fitting = True
    live_dead_objective = "negative_binomial"
    dfdctp_tail_fallback_half_life_days = 1.0
    pk_censored_strategy = "nan"
    pk_sheets = import_and_clean_pkpd(censored_strategy=pk_censored_strategy)
    print_pk_workbook_summary(pk_sheets)

    dfdctp_signal_curve_by_ploidy = build_preferred_dfdctp_signal_surfaces(
        pk_sheets,
        ploidy_keys=("2N", "4N"),
        fallback_half_life_days=dfdctp_tail_fallback_half_life_days,
    )
    for ploidy_key, driver in dfdctp_signal_curve_by_ploidy.items():
        assert_preferred_live_dead_driver(driver)
    print(f"\n{'='*60}")
    print("Built PK-driven live/dead signal curves")
    print(f"  analyte driver: dFdCTP (ng/mL) -> uM")
    print(f"  censored PK strategy: {pk_censored_strategy}")
    print("  live/dead model driver: intracellular dFdCTP surface; parent Gemcitabine PK not used.")
    print(f"  live/dead fit objective: {live_dead_objective}")
    print(f"{'='*60}")
        
    df = assemble_modeling_dataset()

    control_2N = get_aligned_live_dead_data(df, gem_dose="0 nM", ploidy="2N")
    control_4N = get_aligned_live_dead_data(df, gem_dose="0 nM", ploidy="4N")

    if control_2N["dropped_timepoints"] > 0:
        print(f"2N control: dropped {control_2N['dropped_timepoints']} non-overlapping Alive/Dead time points")
    if control_4N["dropped_timepoints"] > 0:
        print(f"4N control: dropped {control_4N['dropped_timepoints']} non-overlapping Alive/Dead time points")

    t_2N, y_alive_2N, y_dead_2N = control_2N["t"], control_2N["y_alive"], control_2N["y_dead"]
    t_4N, y_alive_4N, y_dead_4N = control_4N["t"], control_4N["y_alive"], control_4N["y_dead"]
    
    params_2N_baseline = fit_baseline_shared_rk_replicate_n0(
        t_2N, y_alive_2N, y_dead_2N, "2N Cohort", t_max=3.8, output_dir=output_dir
    )
    params_4N_baseline = fit_baseline_shared_rk_replicate_n0(
        t_4N, y_alive_4N, y_dead_4N, "4N Cohort", t_max=3.8, output_dir=output_dir
    )
    
    ploidy_options = [p for p in df['ploidy'].unique() if pd.notna(p)]
    gem_doses = sorted([d for d in df['gem'].unique() if pd.notna(d) and d != "0 nM"], key=lambda x: float(x.split()[0]))
    live_dead_dose_uM_values = [float(dose.split()[0]) / 1000.0 for dose in gem_doses]

    for ploidy_key in ["2N", "4N"]:
        print_dfdctp_signal_curve_summary(
            ploidy_key,
            dfdctp_signal_curve_by_ploidy[ploidy_key],
            live_dead_dose_uM_values=live_dead_dose_uM_values,
        )
        plot_dfdctp_signal_curve(ploidy_key, dfdctp_signal_curve_by_ploidy[ploidy_key], output_dir=output_dir)
        plot_dfdctp_amplitude_scaling(ploidy_key, dfdctp_signal_curve_by_ploidy[ploidy_key], output_dir=output_dir)
    
    t_max = 5 

    for ploidy in ploidy_options:
        print(f"\n{'='*60}")
        print(f"3-PARAM FIT | COHORT: {ploidy}")
        print(f"{'='*60}")
        
        if ploidy == "2N":
            r_opt, K_opt = params_2N_baseline["r"], params_2N_baseline["K"]
            dfdctp_signal_curve = dfdctp_signal_curve_by_ploidy["2N"]
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
            dfdctp_signal_curve = dfdctp_signal_curve_by_ploidy["4N"]
        else:
            continue
            
        dose_data_list = []

        print(f"Gathering and truncating joint data...")
        for gem_dose in gem_doses:
            try:
                aligned = get_aligned_live_dead_data(df, gem_dose=gem_dose, ploidy=ploidy, t_max=t_max)
                if aligned["dropped_timepoints"] > 0:
                    print(f"  {gem_dose}: dropped {aligned['dropped_timepoints']} non-overlapping Alive/Dead time points")
                t_fit = aligned["t"]
                y_alive_fit = aligned["y_alive"]
                y_dead_fit = aligned["y_dead"]
                
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
            best_objective_value, best_fit, best_n_tr = np.inf, None, 1
            n_range = range(2, 10)

            for n_test in n_range:
                fit_result = fit_live_dead_model(
                    dose_data_list=dose_data_list,
                    r_opt=r_opt,
                    K_opt=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_test,
                    objective=live_dead_objective,
                    fit_means_only=True,
                    high_dose_weight=1.0,
                    max_nfev=3000,
                )
                if fit_result is None:
                    continue
                objective_value = float(fit_result["summary"]["objective_value"])
                if objective_value < best_objective_value:
                    best_objective_value, best_fit, best_n_tr = objective_value, fit_result, n_test

            if best_fit is None:
                print(f"No successful {live_dead_objective} fit found for {ploidy}.")
                continue

            k_tr_opt, k_kill_opt, k_clear_opt = best_fit["params_3"]
            fit_summary = best_fit["summary"]
            
            print(f"Optimal n_tr: {best_n_tr}")
            print(f"Optimal k_tr: {k_tr_opt:.4f}")
            print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
            print(f"Optimal k_clear: {k_clear_opt:.4f}")
            print(f"Fit objective: {live_dead_objective}")
            print(f"Observations used: {fit_summary['n_observations']}")
            print(f"Fitted parameter count: {fit_summary['n_parameters']}")
            if fit_summary["nll"] is not None:
                print(f"Best negative log-likelihood: {fit_summary['nll']:.4f}")
                print(f"AIC: {fit_summary['aic']:.4f}")
                print(f"BIC: {fit_summary['bic']:.4f}")
                print(f"Alive dispersion theta: {fit_summary['theta_alive']:.4f}")
                print(f"Dead dispersion theta: {fit_summary['theta_dead']:.4f}")
            else:
                print(f"Best raw-count least-squares cost: {fit_summary['objective_value']:.4f}")
            print(f"RMSE: {fit_summary['rmse']:.4f}")

            plot_global_fit_subplots_joint(
                dose_data_list, ploidy, r_opt, K_opt, best_n_tr, 
                dfdctp_signal_curve, k_tr_opt, k_kill_opt, k_clear_opt,
                fit_summary=fit_summary,
                output_dir=output_dir,
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
                dfdctp_signal_curve = dfdctp_signal_curve_by_ploidy["2N"]
            elif ploidy == "4N":
                r_opt, K_opt = params_4N_baseline["r"], params_4N_baseline["K"]
                dfdctp_signal_curve = dfdctp_signal_curve_by_ploidy["4N"]
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
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    ploidy=ploidy,
                    rep_idx=plate_row,
                    objective=live_dead_objective,
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
        
