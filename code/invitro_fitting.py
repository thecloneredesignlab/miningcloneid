
import math
import re
import sys
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import least_squares, minimize
from scipy.special import gammaln
from scipy.stats import linregress

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[1]  

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

ALLOWED_MODEL_VARIANTS = {
    "delayed_death_only",
    "immediate_cytostasis_delayed_death",
    "delayed_cytostasis_delayed_death",
}


def validate_model_variant(model_variant: str) -> str:
    if model_variant not in ALLOWED_MODEL_VARIANTS:
        raise ValueError(
            "model_variant must be one of "
            f"{sorted(ALLOWED_MODEL_VARIANTS)}, got {model_variant!r}"
        )
    return model_variant


BASE_PARAMETER_NAMES = ("r", "K", "k_tr", "k_kill", "k_clear")


def get_parameter_names(model_variant: str) -> Tuple[str, ...]:
    model_variant = validate_model_variant(model_variant)
    if model_variant == "delayed_death_only":
        return BASE_PARAMETER_NAMES
    return BASE_PARAMETER_NAMES + ("k_cyto",)


def get_treatment_parameter_names(model_variant: str) -> Tuple[str, ...]:
    parameter_names = get_parameter_names(model_variant)
    return tuple(name for name in parameter_names if name not in {"r", "K"})


def unpack_treatment_params(params: Sequence[float], model_variant: str) -> Dict[str, float]:
    model_variant = validate_model_variant(model_variant)
    expected = 3 if model_variant == "delayed_death_only" else 4
    params_arr = np.asarray(params, dtype=float)
    if len(params_arr) != expected:
        raise ValueError(
            f"Expected {expected} treatment parameters for {model_variant}, got {len(params_arr)}"
        )
    out = {
        "k_tr": float(params_arr[0]),
        "k_kill": float(params_arr[1]),
        "k_clear": float(params_arr[2]),
    }
    if model_variant != "delayed_death_only":
        out["k_cyto"] = float(params_arr[3])
    return out


@dataclass(frozen=True)
class ExperimentPaths:
    project_root: Path
    exp_base: Path
    counts_raw: Path
    counts_agg: Path
    platemap: Path
    pkpd_constants: Path
    output_dir: Path


@dataclass(frozen=True)
class PKConfig:
    dfdctp_molecular_weight_ng_per_nmol: float = 503.1
    censored_strategy: str = "nan"
    fallback_half_life_days: float = 1.0
    reference_dose_by_sheet_uM: Dict[str, float] = field(default_factory=lambda: dict(PKPD_REFERENCE_DOSE_BY_SHEET_UM))


@dataclass(frozen=True)
class JointFitConfig:
    max_days: float = 5.0
    fit_t_max: float = 5.0
    objective: str = "negative_binomial"
    observation_channels: str = "alive_dead"
    model_variant: str = "immediate_cytostasis_delayed_death"
    fit_means_only: bool = False
    high_dose_weight: float = 1.0
    n_tr_values: Tuple[int, ...] = tuple(range(2, 10))
    max_nfev: int = 3000
    large_objective_penalty: float = 1e30
    optimizer_method: str = "L-BFGS-B"
    solver_method: str = "LSODA"
    solver_rtol: float = 1e-6
    solver_atol: float = 1e-8
    prior_sd_log_r: float = 0.75
    prior_sd_log_K: float = 0.75
    prior_sd_log_k_tr: float = 1.00
    prior_sd_log_k_kill: float = 1.00
    prior_sd_log_k_clear: float = 1.00
    prior_sd_log_k_cyto: float = 1.00
    lower_bounds: Dict[str, float] = field(default_factory=lambda: {
        "r": 1e-8,
        "K": 1e-8,
        "k_tr": 1e-4,
        "k_kill": 1e-4,
        "k_clear": 1e-8,
        "k_cyto": 1e-8,
        "theta_alive": 1e-3,
        "theta_dead": 1e-3,
    })
    upper_bounds: Dict[str, float] = field(default_factory=lambda: {
        "r": 100.0,
        "K": 1e9,
        "k_tr": 50.0,
        "k_kill": 500.0,
        "k_clear": 5.0,
        "k_cyto": 1e6,
        "theta_alive": 1e6,
        "theta_dead": 1e6,
    })

    def __post_init__(self):
        validate_model_variant(self.model_variant)


@dataclass
class PKTailFitDiagnostics:
    sheet_name: str
    ploidy: str
    dose_uM: float
    analyte: str
    n_timepoints_total: int
    n_timepoints_positive_finite: int
    peak_idx: int
    peak_time_days: float
    peak_signal_uM: float
    tail_n_points: int
    tail_time_days: Tuple[float, ...]
    tail_signal_uM: Tuple[float, ...]
    tail_fit_used: bool
    fallback_used: bool
    tail_decay_per_day: float
    tail_half_life_days: float
    tail_fit_r2: Optional[float]
    warning_message: Optional[str]


@dataclass
class SimulationResult:
    alive: np.ndarray
    dead: np.ndarray
    success: bool
    message: str = ""


@dataclass
class ReplicateTrajectory:
    ploidy: str
    dose_label: str
    dose_uM: float
    replicate_id: Tuple[str, int]
    t: np.ndarray
    alive: np.ndarray
    dead: np.ndarray
    N0: float
    D0: float


def trim_finite_live_dead_observations(
    t: Sequence[float],
    alive: Sequence[float],
    dead: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Keep only timepoints with finite paired alive/dead observations.

    The joint count-likelihood fit requires paired finite observations for both
    channels at a given timepoint. Missing trailing observations are dropped
    here rather than allowed to propagate into the likelihood as NaNs.
    """
    t_arr = np.asarray(t, dtype=float)
    alive_arr = np.asarray(alive, dtype=float)
    dead_arr = np.asarray(dead, dtype=float)
    if t_arr.shape != alive_arr.shape or t_arr.shape != dead_arr.shape:
        raise ValueError("t, alive, and dead must have matching shapes")
    valid_mask = np.isfinite(t_arr) & np.isfinite(alive_arr) & np.isfinite(dead_arr)
    if not np.any(valid_mask):
        raise ValueError("No finite paired alive/dead observations remain")
    first_valid_idx = int(np.flatnonzero(valid_mask)[0])
    if first_valid_idx != 0:
        raise ValueError("Replicate is missing a finite baseline alive/dead observation at the first timepoint")
    dropped = int(np.count_nonzero(~valid_mask))
    return t_arr[valid_mask], alive_arr[valid_mask], dead_arr[valid_mask], dropped


def default_experiment_paths(project_root: Optional[Path] = None) -> ExperimentPaths:
    root = (project_root or PROJECT_ROOT).resolve()
    exp_base = root / "data" / "InVitroData_Gemcitabine"
    output_dir = SCRIPT_PATH.parent / "invitro_fitting_outputs" / datetime.now().strftime("%Y%m%dT%H%M%S")
    return ExperimentPaths(
        project_root=root,
        exp_base=exp_base,
        counts_raw=exp_base / "processed" / "counts_by_well_time.parquet",
        counts_agg=exp_base / "processed" / "counts_by_well_time_wellAggregated.parquet",
        platemap=exp_base / "Gemcitabine_PlateMap_20240111.xlsx",
        pkpd_constants=exp_base / "drugKinetics" / "GemcitabineExposure_PKPD.xlsx",
        output_dir=output_dir,
    )


def validate_paths(paths: ExperimentPaths, required: Sequence[str]) -> None:
    missing = []
    for attr_name in required:
        if not hasattr(paths, attr_name):
            raise ValueError(f"Unknown ExperimentPaths attribute '{attr_name}'")
        path_value = getattr(paths, attr_name)
        if not Path(path_value).exists():
            missing.append(f"{attr_name}={path_value}")
    if missing:
        raise FileNotFoundError("Missing required experiment paths: " + ", ".join(missing))


_DEFAULT_PATHS = default_experiment_paths()
PATHS = {
    "counts_raw": _DEFAULT_PATHS.counts_raw,
    "counts_agg": _DEFAULT_PATHS.counts_agg,
    "platemap": _DEFAULT_PATHS.platemap,
    "pkpd_constants": _DEFAULT_PATHS.pkpd_constants,
}


def slugify_label(value: str) -> str:
    """Converts labels to filesystem-safe filename fragments."""
    return re.sub(r'[^A-Za-z0-9]+', '_', value.strip()).strip('_').lower()


def normalize_column_name(name: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def parse_well_id(well: Any) -> Tuple[str, int]:
    match = re.fullmatch(r"\s*([A-Ha-h])\s*0*([1-9]|1[0-2])\s*", str(well))
    if match is None:
        raise ValueError(f"Unsupported well identifier '{well}'")
    return match.group(1).upper(), int(match.group(2))


def parse_ploidy_label(value: Any) -> Optional[str]:
    if pd.isna(value):
        return None
    match = re.search(r"(?i)\b([24]N)\b", str(value))
    if match is None:
        return None
    return match.group(1).upper()


def get_pk_reference_dose_uM(sheet_name: str, pk_config: Optional[PKConfig] = None) -> float:
    """Returns the documented administered Gemcitabine dose for a PK workbook sheet."""
    config = pk_config or PKConfig()
    if sheet_name not in config.reference_dose_by_sheet_uM:
        raise KeyError(
            f"No documented PK reference dose configured for sheet '{sheet_name}'. "
            "Update PKConfig.reference_dose_by_sheet_uM to match the workbook metadata."
        )
    return config.reference_dose_by_sheet_uM[sheet_name]


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
    tail_fit_diagnostics: Optional[PKTailFitDiagnostics] = None

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

def import_and_clean_pkpd(
    censored_strategy: Optional[str] = None,
    paths: Optional[ExperimentPaths] = None,
    pk_config: Optional[PKConfig] = None,
):
    """
    Imports PK sheets and parses analyte columns with explicit censored-value handling.

    `censored_strategy="nan"` is the conservative default: BDL/BQL/BLQ-style
    entries become `np.nan` and are ignored by downstream means/fits. Set
    `"zero"` to reproduce the prior zero-imputation behavior explicitly, or
    `"half_lod"` to impute half the reported threshold for `<x` values. The
    original censoring metadata and any numeric censor upper bounds are retained
    in companion columns for downstream PK fitting.
    """
    paths = paths or default_experiment_paths()
    pk_config = pk_config or PKConfig()
    validate_paths(paths, required=("pkpd_constants",))
    strategy = censored_strategy or pk_config.censored_strategy

    raw_sheets = pd.read_excel(paths.pkpd_constants, sheet_name=None)
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
                    censored_strategy=strategy,
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
        "tail_time_days": tuple(float(x) for x in tail_time),
        "tail_conc": tuple(float(x) for x in tail_conc),
        "tail_n_points": int(len(tail_time)),
        "peak_idx": int(peak_idx),
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


def collect_pk_tail_diagnostics(
    curves_by_ploidy: Dict[str, DfdctpSignalSurface]
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for curve in curves_by_ploidy.values():
        for dose, profile in curve.calibration_profiles_by_dose.items():
            diagnostics = profile.tail_fit_diagnostics
            if diagnostics is None:
                continue
            rows.append(
                {
                    "sheet_name": diagnostics.sheet_name,
                    "ploidy": diagnostics.ploidy,
                    "dose_uM": diagnostics.dose_uM,
                    "analyte": diagnostics.analyte,
                    "n_timepoints_total": diagnostics.n_timepoints_total,
                    "n_timepoints_positive_finite": diagnostics.n_timepoints_positive_finite,
                    "peak_idx": diagnostics.peak_idx,
                    "peak_time_days": diagnostics.peak_time_days,
                    "peak_signal_uM": diagnostics.peak_signal_uM,
                    "tail_n_points": diagnostics.tail_n_points,
                    "tail_time_days": diagnostics.tail_time_days,
                    "tail_signal_uM": diagnostics.tail_signal_uM,
                    "tail_fit_used": diagnostics.tail_fit_used,
                    "fallback_used": diagnostics.fallback_used,
                    "tail_decay_per_day": diagnostics.tail_decay_per_day,
                    "tail_half_life_days": diagnostics.tail_half_life_days,
                    "tail_fit_r2": diagnostics.tail_fit_r2,
                    "warning_message": diagnostics.warning_message,
                }
            )
    return pd.DataFrame(rows)


def save_pk_tail_diagnostics(
    curves_by_ploidy: Dict[str, DfdctpSignalSurface],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_df = collect_pk_tail_diagnostics(curves_by_ploidy)
    output_path = output_dir / "pk_tail_fit_diagnostics.tsv"
    diagnostics_df.to_csv(output_path, sep="\t", index=False)
    return output_path

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
        tail_fit_diagnostics = PKTailFitDiagnostics(
            sheet_name=sheet_name,
            ploidy=sheet_name.split("_")[0],
            dose_uM=float(reference_dose_uM),
            analyte=analyte,
            n_timepoints_total=int(len(time_days)),
            n_timepoints_positive_finite=int(np.count_nonzero(np.isfinite(induced_signal_uM) & (induced_signal_uM > 0))),
            peak_idx=int(peak_idx),
            peak_time_days=time_of_peak_days,
            peak_signal_uM=peak_signal_uM,
            tail_n_points=int(tail_fit["tail_n_points"]),
            tail_time_days=tuple(float(x) for x in tail_fit["tail_time_days"]),
            tail_signal_uM=tuple(float(x) for x in tail_fit["tail_conc"]),
            tail_fit_used=True,
            fallback_used=False,
            tail_decay_per_day=tail_decay_per_day,
            tail_half_life_days=tail_half_life_days,
            tail_fit_r2=fit_r2,
            warning_message=None,
        )
    else:
        tail_half_life_days = float(fallback_half_life_days)
        tail_decay_per_day = float(decay_rate_from_half_life_days(tail_half_life_days))
        fit_r2 = None
        warning_message = (
            f"Using fallback dFdCTP tail half-life of {tail_half_life_days:.3f} days "
            f"for {sheet_name} because the terminal fit was not identifiable."
        )
        positive_mask = np.isfinite(induced_signal_uM) & (induced_signal_uM > 0)
        tail_fit_diagnostics = PKTailFitDiagnostics(
            sheet_name=sheet_name,
            ploidy=sheet_name.split("_")[0],
            dose_uM=float(reference_dose_uM),
            analyte=analyte,
            n_timepoints_total=int(len(time_days)),
            n_timepoints_positive_finite=int(np.count_nonzero(positive_mask)),
            peak_idx=int(peak_idx),
            peak_time_days=time_of_peak_days,
            peak_signal_uM=peak_signal_uM,
            tail_n_points=int(np.count_nonzero(positive_mask[peak_idx:])),
            tail_time_days=tuple(float(x) for x in time_days[peak_idx:][positive_mask[peak_idx:]]),
            tail_signal_uM=tuple(float(x) for x in induced_signal_uM[peak_idx:][positive_mask[peak_idx:]]),
            tail_fit_used=False,
            fallback_used=True,
            tail_decay_per_day=tail_decay_per_day,
            tail_half_life_days=tail_half_life_days,
            tail_fit_r2=None,
            warning_message=warning_message,
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
        tail_fit_diagnostics=tail_fit_diagnostics,
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
    model_variant="immediate_cytostasis_delayed_death",
    k_cyto=None,
):
    """
    ODE system driven directly by intracellular dFdCTP signal.

    State order:
    - A: Alive cells
    - Z_1 ... Z_n: transit compartments smoothing/delaying dFdCTP signal
    - D_obs: observable dead cells
    """
    model_variant = validate_model_variant(model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    A = y[0]
    D_obs = y[-1]
    c_now = float(np.asarray(dfdctp_signal_curve(t, dose_muM), dtype=float).reshape(-1)[0])

    dZ = np.zeros(n_tr)
    if n_tr > 0:
        delayed_signal = y[1:1+n_tr]
        dZ[0] = k_tr * (c_now - delayed_signal[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (delayed_signal[i - 1] - delayed_signal[i])
        c_delayed = delayed_signal[-1]
    else:
        c_delayed = c_now

    kappa_death = k_kill * c_delayed
    if model_variant == "delayed_death_only":
        growth_multiplier = 1.0
    else:
        if k_cyto is None or not np.isfinite(k_cyto) or k_cyto < 0:
            raise ValueError(
                f"k_cyto must be finite and >= 0 for model_variant={model_variant}, got {k_cyto}"
            )
        cyto_signal = c_now if model_variant == "immediate_cytostasis_delayed_death" else c_delayed
        growth_multiplier = cytostasis_multiplier(cyto_signal, k_cyto)

    dA = r * growth_multiplier * A * (1 - A / K) - kappa_death * A
    dD_obs = kappa_death * A - k_clear * D_obs
    return [dA] + dZ.tolist() + [dD_obs]


def simulate_joint_dfdctp(
    t,
    params,
    N0,
    D0,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve,
    n_tr,
    model_variant=JointFitConfig().model_variant,
):
    """
    Simulates the live/dead ODE with dFdCTP supplied directly as the drug-driver
    curve instead of integrating a separate intracellular drug state.
    """
    model_variant = validate_model_variant(model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    treatment_params = unpack_treatment_params(params, model_variant)
    y0 = [N0] + [0.0] * n_tr + [D0]
    sol = odeint(
        single_clone_dfdctp_ode_joint,
        y0,
        t,
        args=(
            r,
            K,
            dose_muM,
            dfdctp_signal_curve,
            treatment_params["k_tr"],
            treatment_params["k_kill"],
            treatment_params["k_clear"],
            n_tr,
            model_variant,
            treatment_params.get("k_cyto"),
        ),
    )
    return sol[:, 0], sol[:, -1]


def simulate_joint_dfdctp_safe(
    t: np.ndarray,
    params: Sequence[float],
    N0: float,
    D0: float,
    r: float,
    K: float,
    dose_muM: float,
    dfdctp_signal_curve: DfdctpSignalSurface,
    n_tr: int,
    fit_config: JointFitConfig,
    model_variant: str = JointFitConfig().model_variant,
) -> SimulationResult:
    model_variant = validate_model_variant(model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    t_arr = np.asarray(t, dtype=float)
    if t_arr.ndim != 1 or len(t_arr) == 0:
        return SimulationResult(np.array([]), np.array([]), False, "Empty or non-1D time grid")
    if not np.all(np.isfinite(t_arr)):
        return SimulationResult(np.array([]), np.array([]), False, "Non-finite time grid")
    if np.any(np.diff(t_arr) < 0):
        return SimulationResult(np.array([]), np.array([]), False, "Time grid must be sorted ascending")

    params_arr = np.asarray(params, dtype=float)
    expected_params = 3 if model_variant == "delayed_death_only" else 4
    if len(params_arr) != expected_params or not np.all(np.isfinite(params_arr)):
        return SimulationResult(np.array([]), np.array([]), False, "Invalid treatment params")
    if not np.all(np.isfinite([N0, D0, r, K, dose_muM])) or not np.isfinite(n_tr):
        return SimulationResult(np.array([]), np.array([]), False, "Non-finite simulation inputs")
    if N0 < 0 or D0 < 0 or r < 0 or K <= 0 or dose_muM < 0 or int(n_tr) < 0:
        return SimulationResult(np.array([]), np.array([]), False, "Inputs outside valid ranges")
    try:
        treatment_params = unpack_treatment_params(params_arr, model_variant)
    except ValueError as exc:
        return SimulationResult(np.array([]), np.array([]), False, str(exc))

    bounds_to_check = {
        "k_tr": treatment_params["k_tr"],
        "k_kill": treatment_params["k_kill"],
        "k_clear": treatment_params["k_clear"],
        "r": float(r),
        "K": float(K),
    }
    if "k_cyto" in treatment_params:
        bounds_to_check["k_cyto"] = treatment_params["k_cyto"]
    for name, value in bounds_to_check.items():
        lower = fit_config.lower_bounds.get(name, -np.inf)
        upper = fit_config.upper_bounds.get(name, np.inf)
        if value < (lower - 1e-12) or value > (upper + 1e-12):
            return SimulationResult(np.array([]), np.array([]), False, f"{name}={value} outside bounds")

    y0 = np.array([N0] + [0.0] * int(n_tr) + [D0], dtype=float)

    def rhs(t_eval, state):
        return np.asarray(
            single_clone_dfdctp_ode_joint(
                state,
                t_eval,
                r,
                K,
                dose_muM,
                dfdctp_signal_curve,
                treatment_params["k_tr"],
                treatment_params["k_kill"],
                treatment_params["k_clear"],
                int(n_tr),
                model_variant,
                treatment_params.get("k_cyto"),
            ),
            dtype=float,
        )

    try:
        result = solve_ivp(
            rhs,
            t_span=(float(t_arr[0]), float(t_arr[-1])),
            y0=y0,
            t_eval=t_arr,
            method=fit_config.solver_method,
            rtol=fit_config.solver_rtol,
            atol=fit_config.solver_atol,
        )
    except Exception as exc:
        return SimulationResult(np.array([]), np.array([]), False, f"solve_ivp exception: {exc}")

    if not result.success:
        return SimulationResult(np.array([]), np.array([]), False, result.message)
    if result.y.shape[1] != len(t_arr):
        return SimulationResult(np.array([]), np.array([]), False, "Solver returned unexpected shape")

    alive = np.asarray(result.y[0], dtype=float)
    dead = np.asarray(result.y[-1], dtype=float)
    if not np.all(np.isfinite(alive)) or not np.all(np.isfinite(dead)):
        return SimulationResult(alive, dead, False, "Non-finite predictions")
    if np.any(alive < -1e-8) or np.any(dead < -1e-8):
        return SimulationResult(alive, dead, False, "Materially negative predictions")
    alive = np.clip(alive, 0.0, np.inf)
    dead = np.clip(dead, 0.0, np.inf)
    return SimulationResult(alive, dead, True, result.message)


def safe_objective_value(
    objective_fn: Callable[[np.ndarray], float],
    x: np.ndarray,
    penalty: float,
) -> float:
    try:
        value = objective_fn(x)
    except Exception:
        return penalty
    if not np.isfinite(value):
        return penalty
    return float(value)
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


def _standardize_platemap_frame(plate_long: pd.DataFrame) -> pd.DataFrame:
    required = ["plate_row", "plate_col", "ploidy", "gem", "condition_raw"]
    missing = [col for col in required if col not in plate_long.columns]
    if missing:
        raise ValueError(f"Platemap parser did not produce required columns: {missing}")

    standardized = plate_long[required].copy()
    standardized["plate_row"] = standardized["plate_row"].astype(str).str.upper().str.strip()
    standardized["plate_col"] = pd.to_numeric(standardized["plate_col"], errors="coerce").astype("Int64")
    standardized["ploidy"] = standardized["ploidy"].apply(parse_ploidy_label)
    standardized["gem"] = standardized["gem"].apply(parse_gem_label)

    invalid_coords = standardized[
        ~standardized["plate_row"].isin(list("ABCDEFGH")) |
        standardized["plate_col"].isna() |
        (standardized["plate_col"] < 1) |
        (standardized["plate_col"] > 12)
    ]
    if not invalid_coords.empty:
        raise ValueError(
            "Invalid platemap coordinates encountered: "
            f"{invalid_coords[['plate_row', 'plate_col']].drop_duplicates().to_dict(orient='records')}"
        )

    invalid_ploidy = standardized["ploidy"].notna() & ~standardized["ploidy"].isin(["2N", "4N"])
    if invalid_ploidy.any():
        raise ValueError(
            "Invalid platemap ploidy values encountered: "
            f"{standardized.loc[invalid_ploidy, 'ploidy'].dropna().unique().tolist()}"
        )

    duplicate_mask = standardized.duplicated(subset=["plate_row", "plate_col"], keep=False)
    if duplicate_mask.any():
        conflicting = []
        for key, group in standardized.loc[duplicate_mask].groupby(["plate_row", "plate_col"], dropna=False):
            if len(group.drop_duplicates()) > 1:
                conflicting.append((key, group.to_dict(orient="records")))
        if conflicting:
            raise ValueError(f"Conflicting duplicate platemap entries: {conflicting}")
        standardized = standardized.drop_duplicates()

    return standardized.reset_index(drop=True)

def load_and_clean_counts(
    max_days: float = 5.0,
    paths: Optional[ExperimentPaths] = None,
) -> pd.DataFrame:
    """Loads the aggregated parquet counts data, cleans columns, and filters time."""
    paths = paths or default_experiment_paths()
    validate_paths(paths, required=("counts_agg",))
    df = pd.read_parquet(paths.counts_agg)
    
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

def load_platemap(paths: Optional[ExperimentPaths] = None) -> pd.DataFrame:
    """Loads and parses supported platemap formats into a standardized layout."""
    paths = paths or default_experiment_paths()
    validate_paths(paths, required=("platemap",))
    pm = pd.read_excel(paths.platemap, sheet_name=0)
    pm = pm.rename(columns={col: normalize_column_name(col) for col in pm.columns})
    normalized_cols = list(pm.columns)

    well_col = next((c for c in normalized_cols if c in {"well", "wellid"}), None)
    row_col = next((c for c in normalized_cols if c == "row"), None)
    col_col = next((c for c in normalized_cols if c in {"col", "column"}), None)
    condition_col = next((c for c in normalized_cols if c in {"condition", "condition_raw", "treatment", "label"}), None)
    ploidy_col = next((c for c in normalized_cols if c == "ploidy"), None)
    gem_col = next((c for c in normalized_cols if c in {"gem", "dose"}), None)

    if well_col is not None or (row_col is not None and col_col is not None):
        plate_long = pm.copy()
        if well_col is not None:
            parsed = plate_long[well_col].apply(parse_well_id)
            plate_long["plate_row"] = parsed.apply(lambda item: item[0])
            plate_long["plate_col"] = parsed.apply(lambda item: item[1])
        else:
            plate_long["plate_row"] = plate_long[row_col].astype(str).str.upper().str.strip()
            plate_long["plate_col"] = pd.to_numeric(plate_long[col_col], errors="coerce")

        if condition_col is not None:
            plate_long["condition_raw"] = plate_long[condition_col].astype(str)
            plate_long["ploidy"] = plate_long["condition_raw"].apply(parse_ploidy_label)
            plate_long["gem"] = plate_long["condition_raw"].apply(parse_gem_label)
        elif ploidy_col is not None and gem_col is not None:
            plate_long["ploidy"] = plate_long[ploidy_col]
            plate_long["gem"] = plate_long[gem_col]
            plate_long["condition_raw"] = (
                plate_long[ploidy_col].astype(str).fillna("") + " | " +
                plate_long[gem_col].astype(str).fillna("")
            )
        else:
            raise ValueError(
                "Unsupported long-format platemap. Detected columns "
                f"{normalized_cols} but missing either a combined condition column or separate ploidy/gem columns."
            )
        return _standardize_platemap_frame(plate_long)

    first_col = pm.columns[0]
    row_labels = pm[first_col].astype(str).str.upper().str.strip()
    value_cols = [c for c in pm.columns[1:] if str(c).isdigit() and 1 <= int(str(c)) <= 12]
    if row_labels.isin(list("ABCDEFGH")).any() and value_cols:
        wide = pm.copy()
        wide = wide.rename(columns={first_col: "row"})
        wide["row"] = wide["row"].astype(str).str.upper().str.strip()
        wide = wide[wide["row"].isin(list("ABCDEFGH"))].copy()
        plate_long = pd.melt(
            wide,
            id_vars="row",
            value_vars=value_cols,
            var_name="plate_col",
            value_name="condition_raw",
        )
        plate_long["plate_row"] = plate_long["row"]
        plate_long["ploidy"] = plate_long["condition_raw"].apply(parse_ploidy_label)
        plate_long["gem"] = plate_long["condition_raw"].apply(parse_gem_label)
        return _standardize_platemap_frame(plate_long)

    raise ValueError(
        "Unsupported platemap format. Detected normalized columns "
        f"{normalized_cols}. Expected wide row/1-12 layout or long format with well or row+column fields."
    )

def assemble_modeling_dataset(
    paths: Optional[ExperimentPaths] = None,
    fit_config: Optional[JointFitConfig] = None,
) -> pd.DataFrame:
    """Combines counts and platemap into a single dataframe."""
    paths = paths or default_experiment_paths()
    fit_config = fit_config or JointFitConfig()
    counts_df = load_and_clean_counts(max_days=fit_config.max_days, paths=paths)
    platemap_df = load_platemap(paths=paths)
    
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
# - `k_cyto` is in uM^-1 and acts through a bounded proliferation multiplier

def decay_rate_from_half_life_days(half_life_days: float) -> float:
    """Converts an extracellular half-life in days to a decay rate in days^-1."""
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days


def cytostasis_multiplier(c_dfdctp_uM: float, k_cyto: float) -> float:
    if k_cyto < 0:
        raise ValueError("k_cyto must be >= 0")
    c = max(float(c_dfdctp_uM), 0.0)
    return 1.0 / (1.0 + float(k_cyto) * c)


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
    """Legacy two-stage fit helper: baseline-only r/K fit with replicate-specific N0."""
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
    model_variant=JointFitConfig().model_variant,
    k_cyto=None,
):
    """
    Backward-compatible name for the dFdCTP-driven live/dead ODE.

    The default model variant is now immediate cytostasis plus delayed death,
    not the older delayed-death-only behavior.
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
        model_variant=model_variant,
        k_cyto=k_cyto,
    )


def simulate_joint_ext(
    t,
    params,
    N0,
    D0,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve,
    n_tr,
    model_variant=JointFitConfig().model_variant,
):
    """
    Backward-compatible wrapper for the dFdCTP-driven live/dead simulation.

    The default model variant is now immediate cytostasis plus delayed death.
    """
    return simulate_joint_dfdctp(
        t=t,
        params=params,
        N0=N0,
        D0=D0,
        r=r,
        K=K,
        dose_muM=dose_muM,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr=n_tr,
        model_variant=model_variant,
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
    treatment_params,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
):
    """
    Collects observed/predicted alive and dead counts for objective evaluation.

    `high_dose_weight` is an optional caller-specified dose weight. It is not a
    normalization term.
    """
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )
    fit_config = fit_config or JointFitConfig()
    model_variant = fit_config.model_variant
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    records = []

    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']

        sim_result = simulate_joint_dfdctp_safe(
            t=np.asarray(t_data, dtype=float),
            params=treatment_params,
            N0=N0,
            D0=D0,
            r=r_fixed,
            K=K_fixed,
            dose_muM=dose_muM,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=n_tr_test,
            model_variant=model_variant,
            fit_config=fit_config,
        )
        if not sim_result.success:
            raise RuntimeError(f"Simulation failed for dose {dose_muM}: {sim_result.message}")
        y_alive_pred, y_dead_pred = sim_result.alive, sim_result.dead

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


def residuals_global_joint(treatment_params, dose_data_list, r_fixed, K_fixed, dfdctp_signal_curve, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0,
                           observation_channels="alive_only",
                           model_variant=JointFitConfig().model_variant):
    """
    Residuals for the joint live/dead fit in raw observed cell-count units.

    `high_dose_weight` remains an optional caller-specified dose weight. It is
    not a normalization term and does not rescale residuals by within-dose
    maxima.
    """
    all_residuals = []

    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=n_tr_test,
            fit_config=JointFitConfig(model_variant=model_variant),
            observation_channels=observation_channels,
            fit_means_only=fit_means_only,
            high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        all_residuals.extend((record["alive_obs"] - record["alive_pred"]) * weight)
        if observation_channels == "alive_dead":
            all_residuals.extend((record["dead_obs"] - record["dead_pred"]) * weight)

    return np.array(all_residuals)


def live_dead_objective_nll(
    treatment_params,
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    objective="negative_binomial",
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
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
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )

    if objective == "negative_binomial":
        if theta_alive is None:
            raise ValueError("theta_alive is required for negative_binomial objective")
        if theta_alive <= 0:
            raise ValueError("theta_alive must be > 0")
        if observation_channels == "alive_dead":
            if theta_dead is None:
                raise ValueError("theta_dead is required for alive_dead negative_binomial objective")
            if theta_dead <= 0:
                raise ValueError("theta_dead must be > 0")

    total_nll = 0.0
    fit_config = fit_config or JointFitConfig()
    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr_test,
        fit_config=fit_config,
        observation_channels=observation_channels,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        if objective == "poisson":
            total_nll += weight * poisson_nll(record["alive_obs"], record["alive_pred"])
            if observation_channels == "alive_dead":
                total_nll += weight * poisson_nll(record["dead_obs"], record["dead_pred"])
        else:
            total_nll += weight * negative_binomial_nll(record["alive_obs"], record["alive_pred"], theta_alive)
            if observation_channels == "alive_dead":
                total_nll += weight * negative_binomial_nll(record["dead_obs"], record["dead_pred"], theta_dead)
    return float(total_nll)


def summarize_live_dead_fit(
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr,
    treatment_params,
    objective,
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
    theta_alive=None,
    theta_dead=None,
):
    fit_config = fit_config or JointFitConfig()
    records = collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr,
        fit_config=fit_config,
        observation_channels=observation_channels,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    )
    residual_vector = np.concatenate(
        [
            (record["alive_obs"] - record["alive_pred"]) * record["weight"]
            for record in records
        ] + (
            [
                (record["dead_obs"] - record["dead_pred"]) * record["weight"]
                for record in records
            ] if observation_channels == "alive_dead" else []
        )
    ) if records else np.array([], dtype=float)

    n_obs = int(sum(
        len(record["alive_obs"]) + (len(record["dead_obs"]) if observation_channels == "alive_dead" else 0)
        for record in records
    ))
    rmse = float(np.sqrt(np.mean(residual_vector ** 2))) if residual_vector.size > 0 else np.nan
    summary = {
        "objective": objective,
        "observation_channels": observation_channels,
        "model_variant": fit_config.model_variant,
        "n_observations": n_obs,
        "rmse": rmse,
        "theta_alive": theta_alive,
        "theta_dead": theta_dead if observation_channels == "alive_dead" else None,
    }
    summary.update(unpack_treatment_params(treatment_params, fit_config.model_variant))
    if objective == "least_squares":
        summary["objective_value"] = float(0.5 * np.sum(residual_vector ** 2))
        summary["nll"] = None
        summary["aic"] = None
        summary["bic"] = None
        summary["n_parameters"] = len(get_treatment_parameter_names(fit_config.model_variant))
    else:
        nll = live_dead_objective_nll(
            treatment_params=treatment_params,
            dose_data_list=dose_data_list,
            r_fixed=r_fixed,
            K_fixed=K_fixed,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=n_tr,
            objective=objective,
            fit_config=fit_config,
            observation_channels=observation_channels,
            fit_means_only=fit_means_only,
            high_dose_weight=high_dose_weight,
            theta_alive=theta_alive,
            theta_dead=theta_dead,
        )
        if objective == "negative_binomial":
            n_parameters = len(get_treatment_parameter_names(fit_config.model_variant)) + (
                2 if observation_channels == "alive_dead" else 1
            )
        else:
            n_parameters = len(get_treatment_parameter_names(fit_config.model_variant))
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
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
    max_nfev=3000,
):
    """
    Legacy two-stage fit: fits treatment parameters with baseline r/K held fixed.

    Available objectives:
    - `negative_binomial`
    - `poisson`
    - `least_squares` (legacy raw-residual objective)
    """
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )
    fit_config = fit_config or JointFitConfig()
    model_variant = fit_config.model_variant
    if objective == "least_squares":
        if model_variant == "delayed_death_only":
            guess_grid = [
                [0.2, 10.0, 0.2],
                [0.5, 25.0, 0.5],
                [1.0, 50.0, 1.0],
                [2.0, 100.0, 1.0],
                [5.0, 50.0, 0.5],
                [10.0, 20.0, 0.2],
            ]
        else:
            guess_grid = [
                [0.2, 10.0, 0.2, 1.0],
                [0.5, 25.0, 0.5, 10.0],
                [1.0, 50.0, 1.0, 100.0],
                [2.0, 100.0, 1.0, 10.0],
                [5.0, 50.0, 0.5, 1.0],
                [10.0, 20.0, 0.2, 100.0],
            ]
        treatment_names = get_treatment_parameter_names(model_variant)
        bounds = (
            [fit_config.lower_bounds[name] for name in treatment_names],
            [fit_config.upper_bounds[name] for name in treatment_names],
        )
        best_result = None
        best_summary = None
        best_cost = np.inf
        for guess in guess_grid:
            res = least_squares(
                residuals_global_joint,
                guess,
                bounds=bounds,
                args=(dose_data_list, r_opt, K_opt, dfdctp_signal_curve, n_tr, fit_means_only, high_dose_weight, observation_channels, model_variant),
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
                    treatment_params=res.x,
                    objective=objective,
                    fit_config=fit_config,
                    observation_channels=observation_channels,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                )
        if best_result is None or best_summary is None:
            return None
        return {
            "success": bool(best_result.success),
            "treatment_params": np.asarray(best_result.x, dtype=float),
            "objective": objective,
            "summary": best_summary,
            "optimizer_result": best_result,
        }

    if objective not in {"negative_binomial", "poisson"}:
        raise ValueError(f"Unsupported objective {objective}")

    if model_variant == "delayed_death_only":
        guess_grid = [
            [0.2, 10.0, 0.2],
            [0.5, 25.0, 0.5],
            [1.0, 50.0, 1.0],
            [2.0, 100.0, 1.0],
            [5.0, 50.0, 0.5],
            [10.0, 20.0, 0.2],
        ]
    else:
        guess_grid = [
            [0.2, 10.0, 0.2, 1.0],
            [0.5, 25.0, 0.5, 10.0],
            [1.0, 50.0, 1.0, 100.0],
            [2.0, 100.0, 1.0, 10.0],
            [5.0, 50.0, 0.5, 1.0],
            [10.0, 20.0, 0.2, 100.0],
        ]
    theta_guess_grid = [(20.0, 20.0), (50.0, 50.0), (100.0, 30.0)]
    treatment_names = get_treatment_parameter_names(model_variant)
    param_lower = np.array([fit_config.lower_bounds[name] for name in treatment_names], dtype=float)
    param_upper = np.array([fit_config.upper_bounds[name] for name in treatment_names], dtype=float)
    theta_lower = np.array([1e-3, 1e-3], dtype=float)
    theta_upper = np.array([1e6, 1e6], dtype=float)

    best_result = None
    best_summary = None
    best_nll = np.inf

    for guess in guess_grid:
        theta_guesses = theta_guess_grid if objective == "negative_binomial" else [(None, None)]
        for theta_guess_alive, theta_guess_dead in theta_guesses:
            if objective == "negative_binomial":
                if observation_channels == "alive_dead":
                    x0 = np.log(np.array([guess[0], guess[1], max(guess[2], 1e-8), theta_guess_alive, theta_guess_dead], dtype=float))
                    lower = np.log(np.concatenate([param_lower, theta_lower]))
                    upper = np.log(np.concatenate([param_upper, theta_upper]))
                else:
                    x0 = np.log(np.array([guess[0], guess[1], max(guess[2], 1e-8), theta_guess_alive], dtype=float))
                    lower = np.log(np.concatenate([param_lower, theta_lower[:1]]))
                    upper = np.log(np.concatenate([param_upper, theta_upper[:1]]))
            else:
                x0 = np.log(np.array([guess[0], guess[1], max(guess[2], 1e-8)], dtype=float))
                lower = np.log(param_lower)
                upper = np.log(param_upper)

            def objective_fn(log_params):
                natural_params = np.exp(np.asarray(log_params, dtype=float))
                treatment_params = natural_params[: len(treatment_names)]
                if objective == "negative_binomial":
                    theta_alive_idx = len(treatment_names)
                    theta_alive = float(natural_params[theta_alive_idx])
                    theta_dead = float(natural_params[theta_alive_idx + 1]) if observation_channels == "alive_dead" else None
                else:
                    theta_alive = None
                    theta_dead = None
                return live_dead_objective_nll(
                    treatment_params=treatment_params,
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr_test=n_tr,
                    objective=objective,
                    fit_config=fit_config,
                    observation_channels=observation_channels,
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
                treatment_params = natural_params[: len(treatment_names)]
                theta_alive = float(natural_params[len(treatment_names)]) if objective == "negative_binomial" else None
                theta_dead = (
                    float(natural_params[len(treatment_names) + 1])
                    if objective == "negative_binomial" and observation_channels == "alive_dead"
                    else None
                )
                best_nll = float(res.fun)
                best_result = res
                best_summary = summarize_live_dead_fit(
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_tr,
                    treatment_params=treatment_params,
                    objective=objective,
                    fit_config=fit_config,
                    observation_channels=observation_channels,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                    theta_alive=theta_alive,
                    theta_dead=theta_dead,
                )

    if best_result is None or best_summary is None:
        return None

    return {
        "success": bool(best_result.success),
        "treatment_params": np.asarray(np.exp(np.asarray(best_result.x, dtype=float))[: len(treatment_names)], dtype=float),
        "objective": objective,
        "observation_channels": observation_channels,
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
    k_cyto=None,
    model_variant=JointFitConfig().model_variant,
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
        
    treatment_params = [k_tr, k_kill, k_clear]
    if validate_model_variant(model_variant) != "delayed_death_only":
        treatment_params.append(1e-8 if k_cyto is None else k_cyto)
    
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
                t_data, treatment_params, N0, D0, r, K, dose_muM, dfdctp_signal_curve, n_tr,
                model_variant=model_variant,
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
        f"Model: {fit_summary.get('model_variant', model_variant) if fit_summary is not None else model_variant}",
        f"$k_{{tr}}$ = {k_tr:.3f} | $k_{{kill}}$ = {k_kill:.3f} | $k_{{clear}}$ = {k_clear:.3f}",
    ]
    if validate_model_variant(model_variant) != "delayed_death_only" and k_cyto is not None:
        param_lines.append(f"$k_{{cyto}}$ = {k_cyto:.3f}")
    if fit_summary is not None:
        param_lines.append(f"Fit objective: {fit_summary['objective']}")
        if fit_summary.get("observation_channels") is not None:
            param_lines.append(f"Observation channels: {fit_summary['observation_channels']}")
        if fit_summary.get("theta_alive") is not None and fit_summary.get("theta_dead") is not None:
            param_lines.append(
                f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f} | "
                f"$\\theta_{{dead}}$ = {fit_summary['theta_dead']:.3f}"
            )
        elif fit_summary.get("theta_alive") is not None:
            param_lines.append(f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f}")
        if fit_summary.get("nll") is not None and fit_summary.get("aic") is not None and fit_summary.get("bic") is not None:
            param_lines.append(
                f"NLL = {fit_summary['nll']:.3f} | AIC = {fit_summary['aic']:.3f} | "
                f"BIC = {fit_summary['bic']:.3f}"
            )
        elif fit_summary.get("nll") is not None:
            param_lines.append(f"Data NLL = {fit_summary['nll']:.3f}")
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
    observation_channels="alive_only",
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
            observation_channels=observation_channels,
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

    fitted_treatment = unpack_treatment_params(best_fit["treatment_params"], best_fit["summary"]["model_variant"])
    k_tr_opt = fitted_treatment["k_tr"]
    k_kill_opt = fitted_treatment["k_kill"]
    k_clear_opt = fitted_treatment["k_clear"]
    fit_summary = best_fit["summary"]
    
    print(f"\n--- {ploidy} Replicate {rep_idx} ---")
    print(f"Optimal n_tr: {best_n_tr}")
    print(f"Optimal k_tr: {k_tr_opt:.4f}")
    print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
    print(f"Optimal k_clear: {k_clear_opt:.4f}")
    if "k_cyto" in fitted_treatment:
        print(f"Optimal k_cyto: {fitted_treatment['k_cyto']:.4f}")
    print(f"Fit objective: {objective}")
    print(f"Observation channels: {observation_channels}")
    if fit_summary["nll"] is not None:
        print(f"Best negative log-likelihood: {fit_summary['nll']:.4f}")
        print(f"AIC: {fit_summary['aic']:.4f}")
        print(f"BIC: {fit_summary['bic']:.4f}")
        print(f"Alive dispersion theta: {fit_summary['theta_alive']:.4f}")
        if fit_summary["theta_dead"] is not None:
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
        'k_cyto': fitted_treatment.get('k_cyto', np.nan),
        'cost': best_objective_value,
        'objective': objective,
        'observation_channels': observation_channels,
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


def build_joint_fit_trajectories(
    df: pd.DataFrame,
    ploidies: Sequence[str],
    gem_doses: Sequence[str],
    fit_t_max: Optional[float],
) -> List[ReplicateTrajectory]:
    trajectories: List[ReplicateTrajectory] = []
    total_dropped_nonfinite = 0
    for ploidy in ploidies:
        for gem_dose in gem_doses:
            aligned = get_aligned_live_dead_data(df, gem_dose=gem_dose, ploidy=ploidy, t_max=fit_t_max)
            for rep_idx, replicate_id in enumerate(aligned["replicate_columns"]):
                t_rep, alive, dead, dropped_nonfinite = trim_finite_live_dead_observations(
                    aligned["t"],
                    aligned["y_alive"][:, rep_idx],
                    aligned["y_dead"][:, rep_idx],
                )
                total_dropped_nonfinite += dropped_nonfinite
                trajectories.append(
                    ReplicateTrajectory(
                        ploidy=ploidy,
                        dose_label=gem_dose,
                        dose_uM=float(gem_dose.split()[0]) / 1000.0,
                        replicate_id=(str(replicate_id[0]), int(replicate_id[1])),
                        t=t_rep,
                        alive=alive,
                        dead=dead,
                        N0=float(alive[0]),
                        D0=float(dead[0]),
                    )
                )
    if total_dropped_nonfinite > 0:
        print(
            f"Joint fit trajectories: dropped {total_dropped_nonfinite} non-finite paired "
            "alive/dead observations before likelihood evaluation."
        )
    return trajectories


def fit_joint_partial_pooling_model(
    trajectories: Sequence[ReplicateTrajectory],
    dfdctp_signal_curve_by_ploidy: Dict[str, DfdctpSignalSurface],
    fit_config: JointFitConfig,
    n_tr: int,
) -> Dict[str, Any]:
    """
    Primary MAP fit with log-scale partial pooling across ploidies.
    """
    if len(trajectories) == 0:
        return {
            "success": False,
            "posterior_objective": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "optimizer_message": "No trajectories provided",
            "optimizer_attempts": pd.DataFrame(),
            "model_variant": fit_config.model_variant,
        }
    if fit_config.objective in {"negative_binomial", "poisson"} and fit_config.fit_means_only:
        raise ValueError("fit_means_only=True is not allowed with count-likelihood objectives")
    if fit_config.observation_channels == "alive_only":
        raise ValueError("Partial-pooling primary fit requires observation_channels='alive_dead' to identify k_clear")
    parameter_names = get_parameter_names(fit_config.model_variant)

    ploidies = sorted({traj.ploidy for traj in trajectories})
    if any(ploidy not in dfdctp_signal_curve_by_ploidy for ploidy in ploidies):
        missing = [ploidy for ploidy in ploidies if ploidy not in dfdctp_signal_curve_by_ploidy]
        raise ValueError(f"Missing dFdCTP surfaces for ploidies: {missing}")

    alive_max = max(float(np.nanmax(traj.alive)) for traj in trajectories if len(traj.alive) > 0)
    base_start_grid = [
        {"r": 0.7, "K": max(alive_max, 1.0), "k_tr": 0.5, "k_kill": 25.0, "k_clear": 0.5, "theta_alive": 20.0, "theta_dead": 20.0},
        {"r": 1.2, "K": max(alive_max * 1.5, 1.0), "k_tr": 1.0, "k_kill": 50.0, "k_clear": 1.0, "theta_alive": 50.0, "theta_dead": 50.0},
        {"r": 0.3, "K": max(alive_max * 2.0, 1.0), "k_tr": 2.0, "k_kill": 100.0, "k_clear": 0.2, "theta_alive": 100.0, "theta_dead": 30.0},
    ]
    if fit_config.model_variant == "delayed_death_only":
        mu_start_grid = base_start_grid
    else:
        mu_start_grid = []
        for start_values in base_start_grid:
            for k_cyto in (1.0, 10.0, 100.0):
                expanded = dict(start_values)
                expanded["k_cyto"] = k_cyto
                mu_start_grid.append(expanded)

    prior_sd_by_param = {
        "r": fit_config.prior_sd_log_r,
        "K": fit_config.prior_sd_log_K,
        "k_tr": fit_config.prior_sd_log_k_tr,
        "k_kill": fit_config.prior_sd_log_k_kill,
        "k_clear": fit_config.prior_sd_log_k_clear,
    }
    if "k_cyto" in parameter_names:
        prior_sd_by_param["k_cyto"] = fit_config.prior_sd_log_k_cyto

    def unpack_parameters(log_x: np.ndarray) -> Dict[str, Any]:
        idx = 0
        mu_logs = {}
        for name in parameter_names:
            mu_logs[name] = float(log_x[idx])
            idx += 1
        delta_logs = {ploidy: {} for ploidy in ploidies}
        for ploidy in ploidies:
            for name in parameter_names:
                delta_logs[ploidy][name] = float(log_x[idx])
                idx += 1
        theta_alive = None
        theta_dead = None
        if fit_config.objective == "negative_binomial":
            theta_alive = float(np.exp(log_x[idx]))
            idx += 1
            theta_dead = float(np.exp(log_x[idx]))
            idx += 1

        ploidy_params = {}
        for ploidy in ploidies:
            ploidy_params[ploidy] = {
                name: float(np.exp(mu_logs[name] + delta_logs[ploidy][name]))
                for name in parameter_names
            }
        return {
            "mu_logs": mu_logs,
            "delta_logs": delta_logs,
            "ploidy_params": ploidy_params,
            "theta_alive": theta_alive,
            "theta_dead": theta_dead,
        }

    def data_nll_and_penalty(log_x: np.ndarray) -> Tuple[float, float]:
        unpacked = unpack_parameters(log_x)
        prior_penalty = 0.0
        for ploidy in ploidies:
            for name in parameter_names:
                prior_penalty += 0.5 * (unpacked["delta_logs"][ploidy][name] / prior_sd_by_param[name]) ** 2

        data_nll = 0.0
        for traj in trajectories:
            params = unpacked["ploidy_params"][traj.ploidy]
            sim = simulate_joint_dfdctp_safe(
                t=traj.t,
                params=[params[name] for name in get_treatment_parameter_names(fit_config.model_variant)],
                N0=traj.N0,
                D0=traj.D0,
                r=params["r"],
                K=params["K"],
                dose_muM=traj.dose_uM,
                dfdctp_signal_curve=dfdctp_signal_curve_by_ploidy[traj.ploidy],
                n_tr=n_tr,
                model_variant=fit_config.model_variant,
                fit_config=fit_config,
            )
            if not sim.success:
                raise RuntimeError(f"Simulation failed for {traj.ploidy} {traj.dose_label} {traj.replicate_id}: {sim.message}")
            alive_mask = np.isfinite(traj.alive) & np.isfinite(sim.alive)
            dead_mask = np.isfinite(traj.dead) & np.isfinite(sim.dead)
            if not np.any(alive_mask):
                raise RuntimeError(
                    f"No finite alive observations remain for {traj.ploidy} {traj.dose_label} {traj.replicate_id}"
                )
            if not np.any(dead_mask):
                raise RuntimeError(
                    f"No finite dead observations remain for {traj.ploidy} {traj.dose_label} {traj.replicate_id}"
                )
            if fit_config.objective == "negative_binomial":
                data_nll += negative_binomial_nll(
                    traj.alive[alive_mask],
                    sim.alive[alive_mask],
                    unpacked["theta_alive"],
                )
                data_nll += negative_binomial_nll(
                    traj.dead[dead_mask],
                    sim.dead[dead_mask],
                    unpacked["theta_dead"],
                )
            elif fit_config.objective == "poisson":
                data_nll += poisson_nll(traj.alive[alive_mask], sim.alive[alive_mask])
                data_nll += poisson_nll(traj.dead[dead_mask], sim.dead[dead_mask])
            else:
                data_nll += 0.5 * float(
                    np.sum((traj.alive[alive_mask] - sim.alive[alive_mask]) ** 2)
                    + np.sum((traj.dead[dead_mask] - sim.dead[dead_mask]) ** 2)
                )
        return float(data_nll), float(prior_penalty)

    def posterior_objective(log_x: np.ndarray) -> float:
        data_nll, prior_penalty = data_nll_and_penalty(log_x)
        return data_nll + prior_penalty

    def build_start_vector(start_values: Dict[str, float]) -> np.ndarray:
        values: List[float] = []
        for name in parameter_names:
            values.append(np.log(start_values[name]))
        for _ploidy in ploidies:
            for _name in parameter_names:
                values.append(0.0)
        if fit_config.objective == "negative_binomial":
            values.extend([np.log(start_values["theta_alive"]), np.log(start_values["theta_dead"])])
        return np.asarray(values, dtype=float)

    lower_bounds = []
    upper_bounds = []
    for name in parameter_names:
        lower_bounds.append(np.log(fit_config.lower_bounds[name]))
        upper_bounds.append(np.log(fit_config.upper_bounds[name]))
    for _ploidy in ploidies:
        for _name in parameter_names:
            lower_bounds.append(-5.0)
            upper_bounds.append(5.0)
    if fit_config.objective == "negative_binomial":
        lower_bounds.extend([np.log(fit_config.lower_bounds["theta_alive"]), np.log(fit_config.lower_bounds["theta_dead"])])
        upper_bounds.extend([np.log(fit_config.upper_bounds["theta_alive"]), np.log(fit_config.upper_bounds["theta_dead"])])
    bounds = list(zip(lower_bounds, upper_bounds))

    attempt_rows: List[Dict[str, Any]] = []
    best_result = None
    best_unpacked = None
    best_data_nll = np.nan
    best_prior_penalty = np.nan
    best_value = fit_config.large_objective_penalty

    for start_idx, start_values in enumerate(mu_start_grid):
        x0 = build_start_vector(start_values)

        def safe_posterior(x):
            return safe_objective_value(posterior_objective, np.asarray(x, dtype=float), fit_config.large_objective_penalty)

        try:
            result = minimize(
                safe_posterior,
                x0=x0,
                method=fit_config.optimizer_method,
                bounds=bounds,
                options={"maxiter": fit_config.max_nfev},
            )
            final_value = safe_objective_value(posterior_objective, result.x, fit_config.large_objective_penalty)
            if np.isfinite(final_value) and final_value < fit_config.large_objective_penalty:
                data_nll, prior_penalty = data_nll_and_penalty(result.x)
            else:
                data_nll, prior_penalty = np.nan, np.nan
            attempt_rows.append(
                {
                    "start_idx": start_idx,
                    "start_values": str(start_values),
                    "final_objective": final_value,
                    "success": bool(result.success and np.isfinite(final_value) and final_value < fit_config.large_objective_penalty),
                    "message": result.message,
                }
            )
            if bool(result.success) and np.isfinite(final_value) and final_value < best_value:
                best_value = float(final_value)
                best_result = result
                best_unpacked = unpack_parameters(result.x)
                best_data_nll = data_nll
                best_prior_penalty = prior_penalty
        except Exception as exc:
            attempt_rows.append(
                {
                    "start_idx": start_idx,
                    "start_values": str(start_values),
                    "final_objective": fit_config.large_objective_penalty,
                    "success": False,
                    "message": str(exc),
                }
            )

    attempts_df = pd.DataFrame(attempt_rows)
    n_observations = int(
        sum(np.count_nonzero(np.isfinite(traj.alive)) + np.count_nonzero(np.isfinite(traj.dead)) for traj in trajectories)
    )
    if best_result is None or best_unpacked is None:
        return {
            "success": False,
            "posterior_objective": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "n_observations": n_observations,
            "objective": fit_config.objective,
            "observation_channels": fit_config.observation_channels,
            "model_variant": fit_config.model_variant,
            "theta_alive": np.nan,
            "theta_dead": np.nan,
            "optimizer_message": "All optimizer starts failed",
            "optimizer_attempts": attempts_df,
            "population_parameters": {},
            "ploidy_parameters": {},
        }

    return {
        "success": True,
        "posterior_objective": float(best_value),
        "data_nll": float(best_data_nll),
        "prior_penalty": float(best_prior_penalty),
        "n_observations": n_observations,
        "objective": fit_config.objective,
        "observation_channels": fit_config.observation_channels,
        "model_variant": fit_config.model_variant,
        "theta_alive": best_unpacked["theta_alive"],
        "theta_dead": best_unpacked["theta_dead"],
        "optimizer_message": best_result.message,
        "optimizer_attempts": attempts_df,
        "population_parameters": best_unpacked["mu_logs"],
        "ploidy_parameters": best_unpacked["ploidy_params"],
    }

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

def main(
    paths: Optional[ExperimentPaths] = None,
    pk_config: Optional[PKConfig] = None,
    fit_config: Optional[JointFitConfig] = None,
) -> Dict[str, Any]:
    paths = paths or default_experiment_paths()
    pk_config = pk_config or PKConfig()
    fit_config = fit_config or JointFitConfig()
    validate_paths(paths, required=("counts_agg", "platemap", "pkpd_constants"))
    paths.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving figures to: {paths.output_dir}")

    pk_sheets = import_and_clean_pkpd(paths=paths, pk_config=pk_config)
    print_pk_workbook_summary(pk_sheets)
    dfdctp_signal_curve_by_ploidy = build_preferred_dfdctp_signal_surfaces(
        pk_sheets,
        ploidy_keys=("2N", "4N"),
        fallback_half_life_days=pk_config.fallback_half_life_days,
    )
    save_pk_tail_diagnostics(dfdctp_signal_curve_by_ploidy, paths.output_dir)

    df = assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    ploidy_options = [p for p in df['ploidy'].unique() if pd.notna(p)]
    gem_doses = sorted([d for d in df['gem'].unique() if pd.notna(d)], key=lambda x: float(x.split()[0]))
    live_dead_dose_uM_values = [float(dose.split()[0]) / 1000.0 for dose in gem_doses if dose != "0 nM"]

    for ploidy_key in ["2N", "4N"]:
        print_dfdctp_signal_curve_summary(
            ploidy_key,
            dfdctp_signal_curve_by_ploidy[ploidy_key],
            live_dead_dose_uM_values=live_dead_dose_uM_values,
        )
        plot_dfdctp_signal_curve(ploidy_key, dfdctp_signal_curve_by_ploidy[ploidy_key], output_dir=paths.output_dir)
        plot_dfdctp_amplitude_scaling(ploidy_key, dfdctp_signal_curve_by_ploidy[ploidy_key], output_dir=paths.output_dir)

    trajectories = build_joint_fit_trajectories(
        df=df,
        ploidies=ploidy_options,
        gem_doses=gem_doses,
        fit_t_max=fit_config.fit_t_max,
    )

    summary_rows: List[Dict[str, Any]] = []
    attempt_frames: List[pd.DataFrame] = []
    best_fit = None
    best_n_tr = None
    best_value = np.inf
    for n_tr in fit_config.n_tr_values:
        result = fit_joint_partial_pooling_model(
            trajectories=trajectories,
            dfdctp_signal_curve_by_ploidy=dfdctp_signal_curve_by_ploidy,
            fit_config=fit_config,
            n_tr=n_tr,
        )
        population = result.get("population_parameters", {})
        ploidy_params = result.get("ploidy_parameters", {})
        summary_row = {
            "n_tr": n_tr,
            "success": result["success"],
            "posterior_objective": result["posterior_objective"],
            "data_nll": result["data_nll"],
            "prior_penalty": result["prior_penalty"],
            "n_observations": result["n_observations"],
            "objective": result["objective"],
            "observation_channels": result["observation_channels"],
            "model_variant": result["model_variant"],
            "theta_alive": result["theta_alive"],
            "theta_dead": result["theta_dead"],
            "mu_log_r": population.get("r"),
            "mu_log_K": population.get("K"),
            "mu_log_k_tr": population.get("k_tr"),
            "mu_log_k_kill": population.get("k_kill"),
            "mu_log_k_clear": population.get("k_clear"),
            "mu_log_k_cyto": population.get("k_cyto"),
            "optimizer_message": result["optimizer_message"],
        }
        for ploidy in ["2N", "4N"]:
            params = ploidy_params.get(ploidy, {})
            for name in BASE_PARAMETER_NAMES + ("k_cyto",):
                summary_row[f"{ploidy}_{name}"] = params.get(name, np.nan)
        summary_rows.append(summary_row)
        attempts = result["optimizer_attempts"].copy()
        if not attempts.empty:
            attempts["n_tr"] = n_tr
            attempt_frames.append(attempts)
        if result["success"] and result["posterior_objective"] < best_value:
            best_value = float(result["posterior_objective"])
            best_fit = result
            best_n_tr = n_tr

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(paths.output_dir / "joint_fit_summary.tsv", sep="\t", index=False)
    attempts_df = pd.concat(attempt_frames, ignore_index=True) if attempt_frames else pd.DataFrame()
    attempts_df.to_csv(paths.output_dir / "optimizer_attempts.tsv", sep="\t", index=False)

    if best_fit is None or best_n_tr is None:
        print("No successful joint partial-pooling fit found.")
        return {"success": False, "joint_fit_summary": summary_df, "optimizer_attempts": attempts_df}

    print(f"Best joint partial-pooling fit used n_tr={best_n_tr}")
    print(f"Posterior objective = data_nll + prior_penalty = {best_fit['posterior_objective']:.4f}")
    print(f"Data NLL: {best_fit['data_nll']:.4f}")
    print(f"Prior penalty: {best_fit['prior_penalty']:.4f}")
    print(f"Objective: {best_fit['objective']}")
    print(f"Observation channels: {best_fit['observation_channels']}")
    print(f"Model variant: {best_fit['model_variant']}")

    for ploidy in ["2N", "4N"]:
        params = best_fit["ploidy_parameters"][ploidy]
        param_line = (
            f"{ploidy}: r={params['r']:.4f}, K={params['K']:.4f}, "
            f"k_tr={params['k_tr']:.4f}, k_kill={params['k_kill']:.4f}, k_clear={params['k_clear']:.4f}"
        )
        if "k_cyto" in params:
            param_line += f", k_cyto={params['k_cyto']:.4f}"
        print(param_line)

        dose_data_list = []
        for gem_dose in [dose for dose in gem_doses if dose != "0 nM"]:
            try:
                aligned = get_aligned_live_dead_data(df, gem_dose=gem_dose, ploidy=ploidy, t_max=fit_config.fit_t_max)
                dose_data_list.append(
                    {
                        "dose_label": gem_dose,
                        "dose_muM": float(gem_dose.split()[0]) / 1000.0,
                        "t": aligned["t"],
                        "y_alive": aligned["y_alive"],
                        "y_dead": aligned["y_dead"],
                        "N0": float(np.nanmean(aligned["y_alive"][0, :])),
                        "D0": float(np.nanmean(aligned["y_dead"][0, :])),
                    }
                )
            except ValueError:
                continue
        if dose_data_list:
            plot_global_fit_subplots_joint(
                dose_data_list=dose_data_list,
                ploidy_label=ploidy,
                r=params["r"],
                K=params["K"],
                n_tr=best_n_tr,
                dfdctp_signal_curve=dfdctp_signal_curve_by_ploidy[ploidy],
                k_tr=params["k_tr"],
                k_kill=params["k_kill"],
                k_clear=params["k_clear"],
                k_cyto=params.get("k_cyto"),
                model_variant=fit_config.model_variant,
                fit_summary={
                    "objective": best_fit["objective"],
                    "observation_channels": best_fit["observation_channels"],
                    "model_variant": best_fit["model_variant"],
                    "theta_alive": best_fit["theta_alive"],
                    "theta_dead": best_fit["theta_dead"],
                    "nll": best_fit["data_nll"],
                    "aic": None,
                    "bic": None,
                    "rmse": np.nan,
                },
                output_dir=paths.output_dir,
            )

    return {
        "success": True,
        "best_fit": best_fit,
        "best_n_tr": best_n_tr,
        "joint_fit_summary": summary_df,
        "optimizer_attempts": attempts_df,
    }


if __name__ == "__main__":
    main()
        
