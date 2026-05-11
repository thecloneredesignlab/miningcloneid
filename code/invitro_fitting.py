
import math
import re
import sys
import warnings
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

for name, path in PATHS.items():
    if not path.exists():
        print(f"Warning: {name} not found at {path}")


#################### eta, K_decay fitting utilities ####################

def import_and_clean_pkpd():
    """Imports PK sheets and cleans dFdCTP data (handles BDL/BQL automatically)."""
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
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0.0)
            
        clean_data_dict[sheet_name] = df
    return clean_data_dict

def get_r_eta_parameters(sheets, gem_molec_wt=503.1):
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
            eta_val = extract_eta(sheets[cohort], dose_muM=1.0, gem_molec_wt=gem_molec_wt)
            
            results[cohort] = {
                'eta': eta_val,
                'k_decay': k_val * 24 if k_val else 0.0, # Convert 1/hr to 1/day
                'hl_hours': t_half,
                'r2_confidence': r2
            }
    return results

def extract_half_life(df: pd.DataFrame, analyte: str = 'dFdCTP (ng/mL)') -> Tuple[Optional[float], Optional[float], Optional[float]]:
    mean_data = df.groupby('Timepoint')[analyte].mean().reset_index()
    mean_data = mean_data[mean_data[analyte] > 0]
    if len(mean_data) < 2: return None, None, None

    peak_idx = mean_data[analyte].idxmax()
    decay_phase = mean_data[mean_data['Timepoint'] >= mean_data.loc[peak_idx, 'Timepoint']]
    if len(decay_phase) < 2: return None, None, None

    slope, _, r_val, _, _ = linregress(decay_phase['Timepoint'], np.log(decay_phase[analyte]))
    k = -slope
    if k <= 0: return None, None, None 
    
    return round(k, 4), round(np.log(2)/k, 2), round(r_val**2, 4)

def extract_eta(df: pd.DataFrame, dose_muM: float = 1.0, analyte: str = 'dFdCTP (ng/mL)', gem_molec_wt: float = 503.1) -> Optional[float]:
    mean_data = df.groupby('Timepoint')[analyte].mean()
    if 0.0 not in mean_data.index or 1.0 not in mean_data.index: return None

    delta_c_hourly = mean_data[1.0] - mean_data[0.0]
    eta_daily = (delta_c_hourly / (0.00971 * dose_muM) ) * 24 / gem_molec_wt
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


def residuals_fixed_N0(params, t, y_true, N0_fixed):
    """Objective function that only fits r and K, keeping N0 fixed."""
    r, K = params
    y_pred = logistic_growth(t, N0_fixed, r, K)
    return y_true - y_pred

def fit_baseline_locked_N0(t_data, y_alive, y_dead, ploidy_label, t_max=None):
    
    if y_alive.ndim > 1:
        t_flat = np.repeat(t_data, y_alive.shape[1]) 
        y_alive_flat = y_alive.flatten()
        y_dead_flat = y_dead.flatten()
    else:
        t_flat = t_data
        y_alive_flat = y_alive
        y_dead_flat = y_dead
        
    valid_mask = ~np.isnan(y_alive_flat)
    
    if t_max is not None:
        time_mask = t_flat <= t_max
        combined_mask = valid_mask & time_mask
    else:
        combined_mask = valid_mask

    t_fit = t_flat[combined_mask]
    y_alive_fit = y_alive_flat[combined_mask]
    
    N0_fixed = np.nanmean(y_alive[0, :]) if y_alive.ndim > 1 else y_alive[0]
    
    K_guess = np.nanmax(y_alive_fit)
    r_guess = 0.5
    
    res = least_squares(residuals_fixed_N0, [r_guess, K_guess], 
                        args=(t_fit, y_alive_fit, N0_fixed), 
                        bounds=([0, 0], [np.inf, np.inf]))
    r_opt, K_opt = res.x
    print(f"\n{'='*60}")
    print(f"{ploidy_label} (LOCKED N0) ")
    if t_max is not None:
        print(f"Fitted only on t <= {t_max} Days")
    print(f"Fixed N0 = {N0_fixed:.1f}")
    print(f"r = {r_opt:.4f} days^-1")
    print(f"K = {K_opt:.1f}")
    print("=" * 60 + "\n")
    
    fig, ax1 = plt.subplots(figsize=(8, 5))
    color1 = 'tab:purple'
    ax1.set_xlabel('Time (Days)')
    ax1.set_ylabel('Alive Cell Count', color=color1)
    
    if y_alive.ndim > 1:
        ax1.plot(t_flat, y_alive_flat, 'o', color=color1, alpha=0.15, label='Alive (Replicates)')
        ax1.plot(t_data, np.nanmean(y_alive, axis=1), 's', color='indigo', label='Alive (Mean)')
    else:
        ax1.plot(t_data, y_alive, 's', color='indigo', label='Alive')
    
    t_smooth = np.linspace(min(t_data), max(t_data), 100)
    y_smooth = logistic_growth(t_smooth, N0_fixed, r_opt, K_opt)
    ax1.plot(t_smooth, y_smooth, '-', color='fuchsia', linewidth=2, label=f'Fit (r={r_opt:.3f})')
    
    if t_max is not None:
        ax1.axvline(x=t_max, color='gray', linestyle='--', alpha=0.7, label=f'Fit Cutoff (t={t_max})')
        
    ax1.tick_params(axis='y', labelcolor=color1)
    
    plt.title(f"Phase 1 Control Fit: {ploidy_label} (Locked N0)")
    plt.legend(loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return [N0_fixed, r_opt, K_opt]



#################### Cohort fitting utilities ####################

def single_clone_signal_ode_joint(y, t, r, K, dose_muM, eta, k_decay, k_tr, k_kill, k_clear, n_tr):
    """
    ODE system tracking:
    - A: Alive cells (y[0])
    - C_dna: DNA-incorporated active drug (y[1])
    - Z_1 ... Z_n: Transit compartments (y[2 : 2+n_tr])
    - D_obs: Observable Dead cells (y[-1])
    """
    A = y[0]
    C_dna = y[1]
    D_obs = y[-1]
    
    C_t =  0.00971 * dose_muM 
    
    # Active drug kinetics
    dC_dna = eta * C_t - k_decay * C_dna
    
    # Transit compartments
    dZ = np.zeros(n_tr)
    if n_tr > 0:
        Z_comps = y[2:2+n_tr]
        dZ[0] = k_tr * (C_dna - Z_comps[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (Z_comps[i-1] - Z_comps[i])
            
        kappa = k_kill * Z_comps[-1]
    else:
        kappa = k_kill * C_dna
        
    # Alive cells: Logistic growth minus drug-induced death
    dA = r * A * (1 - A / K) - kappa * A
    
    # Observable dead cells: Accumulate from death, deplete via lysis/clearance
    dD_obs = kappa * A - k_clear * D_obs
    
    return [dA, dC_dna] + dZ.tolist() + [dD_obs]


def simulate_joint_ext(t, params_3, N0, D0, r, K, dose_muM, eta_fixed, k_decay_fixed, n_tr):
    """
    Wrapper to simulate the ODE system. 
    params_3 = [k_tr, k_kill, k_clear]
    """
    k_tr, k_kill, k_clear = params_3 
    
    # Initial conditions: [A0, C_dna_0, Z_1_0, ..., Z_n_0, D_obs_0]
    y0 = [N0, 0.0] + [0.0] * n_tr + [D0] 
    
    sol = odeint(single_clone_signal_ode_joint, y0, t, 
                 args=(r, K, dose_muM, eta_fixed, k_decay_fixed, k_tr, k_kill, k_clear, n_tr))
    
    return sol[:, 0], sol[:, -1]


def residuals_global_joint(params_3, dose_data_list, r_fixed, K_fixed, eta_fixed, k_decay_fixed, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0):
    all_residuals = []
    
    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']
        
        y_alive_pred, y_dead_pred = simulate_joint_ext(
            t_data, params_3, N0, D0, r_fixed, K_fixed, dose_muM, eta_fixed, k_decay_fixed, n_tr_test
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


def plot_global_fit_subplots_joint(dose_data_list, ploidy_label, r, K, n_tr, eta, k_decay_fixed, k_tr, k_kill, k_clear):
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
            A_sim, D_sim = simulate_joint_ext(t_data, ode_params_3, N0, D0, r, K, dose_muM, eta, k_decay_fixed, n_tr)
            
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
        
    param_text = (f"Optima: n_tr = {n_tr} | $\eta$ (fixed) = {eta:.3f} | $k_{{decay}}$ (fixed) = {k_decay_fixed:.3f}\n"
                  f"$k_{{tr}}$ = {k_tr:.3f} | $k_{{kill}}$ = {k_kill:.3f} | $k_{{clear}}$ = {k_clear:.3f}")
    fig.suptitle(f"{ploidy_label} Cohort - Joint Live/Dead Kinetic Fit\n{param_text}", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85) 

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

def fit_joint_one_replicate(dose_data_list, r_opt, K_opt, eta_fixed, k_decay_fixed, ploidy, rep_idx):
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
                                      n_test, False, 1.0),
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

def plot_replicate_parameter_summary(replicate_fit_df):
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

if __name__ == "__main__":
    gemcitabine_mw = 503.1
    pk_sheets = import_and_clean_pkpd()
    fitted_map = get_r_eta_parameters(pk_sheets, gem_molec_wt=gemcitabine_mw)

    try:
        eta_2N, k_decay_2N = fitted_map['2N']['eta'], fitted_map['2N']['k_decay']
        eta_4N, k_decay_4N = fitted_map['4N']['eta'], fitted_map['4N']['k_decay']
        print(f"\n{'='*60}")
        print(f"Successfully extracted PK:")
        print(f"  2N -> eta: {eta_2N}, k_decay: {k_decay_2N:.4f} days⁻¹")
        print(f"  4N -> eta: {eta_4N}, k_decay: {k_decay_4N:.4f} days⁻¹")
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
    
    params_2N_locked = fit_baseline_locked_N0(t_2N, y_alive_2N, y_dead_2N, "2N Cohort", t_max=3.8)
    params_4N_locked = fit_baseline_locked_N0(t_4N, y_alive_4N, y_dead_4N, "4N Cohort", t_max=3.8)
    
    ploidy_options = [p for p in df['ploidy'].unique() if pd.notna(p)]
    gem_doses = sorted([d for d in df['gem'].unique() if pd.notna(d) and d != "0 nM"], key=lambda x: float(x.split()[0]))
    
    t_max = 5 

    for ploidy in ploidy_options:
        print(f"\n{'='*60}")
        print(f"3-PARAM FIT | COHORT: {ploidy}")
        print(f"{'='*60}")
        
        if ploidy == "2N":
            r_opt, K_opt = params_2N_locked[1], params_2N_locked[2]
            eta_fixed = eta_2N
            k_decay_fixed = k_decay_2N
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_locked[1], params_4N_locked[2]
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
                                            n_test, True, 1.0),
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
                eta_fixed, k_decay_fixed, k_tr_opt, k_kill_opt, k_clear_opt
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
            r_opt, K_opt = params_2N_locked[1], params_2N_locked[2]
            eta_fixed = eta_2N
            k_decay_fixed = k_decay_2N
        elif ploidy == "4N":
            r_opt, K_opt = params_4N_locked[1], params_4N_locked[2]
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

    plot_replicate_parameter_summary(replicate_fit_df)
    plt.show()    
        
