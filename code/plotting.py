"""
simulate_treatment.py
─────────────────────
Forward simulation of a fixed treatment schedule.
No beam search. No MCMC / Gibbs sampling.

Configure the sections marked  ◄ EDIT HERE ►  then run.
"""

import numpy as np
import matplotlib.pyplot as plt
from PujanEarlyVersionModel import ploidy_forcast

# Initial ploidy distribution (cells per ploidy class at t = 0)
INITIAL_PLOIDY = {2.0: 1.1343e9, 3.0: 0.1e9, 4.0: 0.1e9}

# Global model parameters
R_BASE = 0.28
K_CAP  = 6.0e10

# R_BASE = 0.2222
# K_CAP  = 1.565e+10


N_SIMS = 1000

# Treatment schedule ──────────────────────────────────────────────────────────
# List of (start_day, end_day, drug_name) in chronological order.
# Gaps between entries are automatically filled with "none" (no drug).
# drug_name must match a key in DRUG_PK below (or be "none").

TREATMENT_SCHEDULE = [
    (  0,   7, "none"),
    (  7,  14, "none"),
    ( 14,  21, "none"),
    ( 21,  28, "none"),
    ( 28,  31, "none"),
]

DRUG_PK = {
    # drug            C_peak   half_life (days)
    "gemcitabine": {"C_peak": 0.032, "half_life": 0.05},
    "bay1895344":  {"C_peak": 0.5,   "half_life": 0.50},
    "alisertib":   {"C_peak": 1.53,  "half_life": 19.0},
    "ispinesib":   {"C_peak": 0.09,  "half_life": 1.04},
    "none":        {},   # no PK params needed
}

_CELLS_PER_CM3 = 1e7

_OBSERVED_TUMOR_BURDENS_CM3 = {
     0:  133.43,
     3: 379.91,
     7: 459.47,
    10: 567.09,
    14: 958.81,
    17: 932.14,
    21: 766.32,
    24: 1441.37,
    28: 1902.76,
    31: 2622.36,
}

OBSERVED_TUMOR_BURDENS = {
    day: vol * _CELLS_PER_CM3
    for day, vol in _OBSERVED_TUMOR_BURDENS_CM3.items()
}

# Plot colours per drug ───────────────────────────────────────────────────────
DRUG_COLORS = {
    "gemcitabine": "orange",
    "bay1895344":  "red",
    "alisertib":   "green",
    "ispinesib":   "blue",
    "volasertib":  "purple",
    "cytarabine":  "brown",
    "none":        "lightgrey",
}

def _fill_gaps(schedule):
    """Insert 'none' cycles for any uncovered gaps in the schedule."""
    filled = []
    for start, end, drug in schedule:
        if filled and filled[-1][1] < start:
            filled.append((filled[-1][1], start, "none"))
        filled.append((start, end, drug))
    return filled


def simulate_cycle(ploidy_status, drug, cycle_len, r_base, k_cap, pk_overrides):
    """Run one treatment cycle; return (new_ploidy_dict, segment_trajectory).

    segment_trajectory shape: (n_timepoints-1, n_ploidies)  — first point
    is omitted so it does not duplicate the previous cycle's endpoint.
    """
    overrides = pk_overrides or {}
    ploidies, _, _, _, Tpaths = ploidy_forcast(
        ploidy_status, drug, cycle_len,
        N_SIMS=N_SIMS, R_BASE=r_base, K_CAP=k_cap, **overrides
    )
    mean_traj = np.mean(Tpaths, axis=0).T        # (n_tp, n_ploidies)
    final     = mean_traj[-1]
    new_status = {ploidies[k]: float(final[k]) for k in range(len(ploidies))}
    return new_status, mean_traj[1:]              # skip index 0 (= start of cycle)


def run_simulation():
    schedule = _fill_gaps(TREATMENT_SCHEDULE)

    ploidy    = dict(INITIAL_PLOIDY)
    all_times = [0.0]
    all_traj  = [np.array(list(INITIAL_PLOIDY.values()))]   # shape (n_ploidies,)
    path_info = []   # list of (drug, start_day, end_day)

    for start_day, end_day, drug in schedule:
        cycle_len = end_day - start_day
        pk        = DRUG_PK.get(drug, {})

        print(f"  Day {start_day:>5.1f} → {end_day:>5.1f}  drug={drug}  "
              f"C_peak={pk.get('C_peak', 0):.4g}  "
              f"half_life={pk.get('half_life', '—')}")

        new_ploidy, seg_traj = simulate_cycle(
            ploidy, drug, cycle_len, R_BASE, K_CAP, pk
        )

        n_pts     = len(seg_traj)
        step      = cycle_len / n_pts
        seg_times = start_day + (np.arange(1, n_pts + 1) * step)

        all_times.extend(seg_times.tolist())
        all_traj.extend(seg_traj.tolist())
        path_info.append((drug, start_day, end_day))

        ploidy = new_ploidy

    time_vec = np.array(all_times)
    traj_arr = np.array(all_traj)          # (n_timepoints, n_ploidies)
    total_burden = traj_arr.sum(axis=1)

    return time_vec, traj_arr, total_burden, path_info


# =============================================================================
#  PLOTTING
# =============================================================================

def plot_results(time_vec, traj_arr, total_burden, path_info):
    ploidy_keys  = list(INITIAL_PLOIDY.keys())
    ploidy_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)

    # ── Panel 1: per-ploidy burden ────────────────────────────────────────────
    ax1 = axes[0]
    _shade_background(ax1, path_info)

    for i, pk in enumerate(ploidy_keys):
        if i < traj_arr.shape[1]:
            ax1.plot(time_vec, traj_arr[:, i],
                     label=f"{int(pk) if pk == int(pk) else pk}n",
                     color=ploidy_colors[i % len(ploidy_colors)],
                     linewidth=2)

    ax1.set_ylabel("Cell Count")
    ax1.set_title("Tumor Burden per Ploidy")
    ax1.legend(loc="upper left")
    ax1.grid(True, ls="--", alpha=0.3)
    ax1.set_yscale("symlog", linthresh=1e5)

    # ── Panel 2: total burden + observations ──────────────────────────────────
    ax2 = axes[1]
    _shade_background(ax2, path_info, add_legend=True)

    ax2.plot(time_vec, total_burden,
             label="Simulated total", color="black", linewidth=2.5)

    if OBSERVED_TUMOR_BURDENS:
        obs_days    = np.array(sorted(OBSERVED_TUMOR_BURDENS))
        obs_burdens = np.array([OBSERVED_TUMOR_BURDENS[d] for d in obs_days])
        ax2.scatter(obs_days, obs_burdens,
                    color="crimson", zorder=5, s=60,
                    label="Observed", marker="D")

    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Total Cell Count")
    ax2.set_title("Total Tumor Burden")
    ax2.legend(loc="upper left")
    ax2.grid(True, ls="--", alpha=0.3)
    ax2.set_yscale("symlog", linthresh=1e5)

    plt.tight_layout()
    plt.savefig("simulation_output.png", dpi=150, bbox_inches="tight")
    print("\nPlot saved → simulation_output.png")
    plt.show()


def _shade_background(ax, path_info, add_legend=False):
    """Shade each treatment cycle on a given axes."""
    shaded = set()
    for drug, start_day, end_day in path_info:
        color = DRUG_COLORS.get(drug, "lightgrey")
        label = drug if (add_legend and drug not in shaded) else None
        ax.axvspan(start_day, end_day, color=color, alpha=0.15,
                   linewidth=0, label=label)
        shaded.add(drug)


# =============================================================================
#  ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    print("Running forward simulation...")
    time_vec, traj_arr, total_burden, path_info = run_simulation()
    print(f"\nFinal total burden: {total_burden[-1]:.3e} cells")
    plot_results(time_vec, traj_arr, total_burden, path_info)