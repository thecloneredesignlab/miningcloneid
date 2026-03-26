import numpy as np
import matplotlib.pyplot as plt
from Missegregation_Model import (
    simulate_karyotype_ode_piecewise,
    get_concentration_curve,
    build_times_with_doses,
    f,
)

# =============================================================================
#  CONFIGURATION — edit these to change the simulation
# =============================================================================

# Initial ploidy distribution (cells per chromosome count at t = 0).
# Keys are chromosome numbers: 46 = diploid (2n), 69 = triploid (3n), 92 = tetraploid (4n).
INITIAL_PLOIDY = {46: 0.94e6, 69: 0.05e6, 92: 0.01e6}
# INITIAL_PLOIDY = {46: 0.01e6, 69: 0.05e6, 92: 0.94e6}  # for 4n start

HAPLOID_N: int = 23   # chromosomes in a haploid complement

# Global model parameters
R_BASE    = 0.2528
K_CAP     = 1.822e+10
BETA_INIT = 0.01138
#MAP  R=0.2528  K=1.822e+10  beta=0.01138

ODE_STEP = 0.05       # ODE integration step (days)

# Minimum viable tumor burden. If total cell count drops below this threshold
# the simulation is halted, the crossing day is reported, and only data up to
# that point is plotted.
MIN_SIZE: float = 1e5

# Reference dose for C_peak scaling (mg/kg).
# C_peak in DRUG_PK is calibrated at this dose.
# For any other dose d:  C_peak_effective = C_peak_ref * (d / DOSE_REFERENCE_MG_KG).
DOSE_REFERENCE_MG_KG: float = 120.0

# Treatment schedule ──────────────────────────────────────────────────────────
# Each entry is either:
#   (start_day, end_day, drug_name)              ← uses default C_peak from DRUG_PK
#   (start_day, end_day, drug_name, dose_mg_kg)  ← scales C_peak by dose/DOSE_REFERENCE_MG_KG
#
# Gaps between entries are automatically filled with "none" segments.
# drug_name must match a key in DRUG_PK below (or be "none").

TREATMENT_SCHEDULE = [
    ( 0,  63, "bay1895344")
]

# Default PK parameters per drug.
# C_peak is the reference peak plasma concentration at DOSE_REFERENCE_MG_KG.
DRUG_PK = {
    "gemcitabine": {"C_peak": 0.032,  "half_life": 0.05,  "period": 7.0},
    "bay1895344":  {"C_peak": 0.5,    "half_life": 0.50,  "period": 0.5},
    "alisertib":   {"C_peak": 1.53,   "half_life": 19.0,  "period": 7.0},
    "ispinesib":   {"C_peak": 0.09,   "half_life": 1.04,  "period": 7.0},
    "none":        {},
}

# Plot colours per drug ───────────────────────────────────────────────────────
DRUG_COLORS = {
    "gemcitabine": "orange",
    "bay1895344":  "red",
    "alisertib":   "green",
    "ispinesib":   "blue",
    "volasertib":  "purple",
    "cytarabine":  "brown",
    "none":        "yellow",
}


# =============================================================================
#  SIMULATION HELPERS
# =============================================================================

def _normalise_schedule(schedule):
    """Ensure every entry is a 4-tuple (start, end, drug, dose_or_None).

    A 3-tuple entry means "use default dose" → dose is stored as None.
    A 4-tuple entry carries an explicit dose (mg/kg).
    """
    out = []
    for entry in schedule:
        if len(entry) == 3:
            start, end, drug = entry
            out.append((start, end, drug, None))
        else:
            start, end, drug, dose = entry
            out.append((start, end, drug, dose))
    return out


def _fill_gaps(schedule):
    """Insert 'none' segments to cover any gaps in the normalised schedule."""
    filled = []
    for start, end, drug, dose in schedule:
        if filled and filled[-1][1] < start:
            filled.append((filled[-1][1], start, "none", None))
        filled.append((start, end, drug, dose))
    return filled


def _pk_overrides_for(drug: str, dose_mg_kg=None) -> dict:
    """Return PK parameters for *drug*, optionally scaling C_peak by dose.

    If *dose_mg_kg* is None the default C_peak from DRUG_PK is used unchanged.
    Otherwise C_peak is scaled linearly:
        C_peak_effective = C_peak_ref * (dose_mg_kg / DOSE_REFERENCE_MG_KG)
    """
    base = dict(DRUG_PK.get(drug, {}))
    if dose_mg_kg is not None and "C_peak" in base and DOSE_REFERENCE_MG_KG > 0:
        base["C_peak"] = base["C_peak"] * (dose_mg_kg / DOSE_REFERENCE_MG_KG)
    return base


def _run_ode_cycle(ploidy_status, drug, cycle_len, r_base, k_cap,
                   pk_overrides, beta=BETA_INIT):
    """Run one ODE cycle.

    Returns
    -------
    new_status : dict {chr_count: cell_count}  — state at end of cycle
    Ns         : array of chromosome counts tracked by the ODE
    traj       : ndarray (n_timepoints-1, n_chr) — time-series excluding t=0
    times_rel  : ndarray (n_timepoints-1,)       — times relative to cycle start
    """
    overrides = pk_overrides or {}
    C_fn  = get_concentration_curve(drug, **overrides)
    TIMES = build_times_with_doses(
        (0.0, cycle_len), ODE_STEP,
        drug_name=drug, include_days=True, eps=1e-8,
    )
    _t, Ns, T_mat, _T_total, _M = simulate_karyotype_ode_piecewise(
        ploidy_status, drug,
        t_span=(0.0, cycle_len),
        r=r_base, Kcap=k_cap, beta=beta,
        N_min=10, N_max=90,
        C_fn=C_fn, f_param_fn=f, t_eval=TIMES, max_step=ODE_STEP,
        renormalize_M=True,
    )
    # T_mat: (n_chr, n_timepoints) → transpose to (n_tp, n_chr)
    traj       = T_mat.T
    final      = traj[-1]
    new_status = {int(N): float(c) for N, c in zip(Ns, final) if c > 0}
    return new_status, np.asarray(Ns), traj[1:], np.asarray(TIMES[1:])


# =============================================================================
#  MAIN SIMULATION
# =============================================================================

def run_simulation():
    """Run the full treatment schedule, stopping early if total burden drops
    below MIN_SIZE.

    Returns
    -------
    burden_timeline  : list of (day, total_cells)
    ploidy_snapshots : list of (day, {chr_count: cell_count})
    path_info        : list of (drug, start_day, end_day)
    extinction_day   : float or None — day burden crossed MIN_SIZE, else None
    """
    schedule = _fill_gaps(_normalise_schedule(TREATMENT_SCHEDULE))
    ploidy   = dict(INITIAL_PLOIDY)
    path_info     = []
    extinction_day = None

    # Seed with t=0 state
    burden_timeline  = [(0.0, float(sum(ploidy.values())))]
    ploidy_snapshots = [(0.0, dict(ploidy))]

    outer_done = False  # flag to break out of the segment loop early

    for start_day, end_day, drug, dose in schedule:
        if outer_done:
            break

        cycle_len = end_day - start_day
        pk        = _pk_overrides_for(drug, dose)

        dose_str = f"{dose:.4g} mg/kg" if dose is not None else "default"
        print(f"  Day {start_day:>5.1f} → {end_day:>5.1f}  drug={drug}  "
              f"dose={dose_str}  "
              f"C_peak={pk.get('C_peak', 0):.4g}  "
              f"half_life={pk.get('half_life', '—')}")

        new_ploidy, Ns, traj, times_rel = _run_ode_cycle(
            ploidy, drug, cycle_len, R_BASE, K_CAP, pk
        )

        for t_rel, row in zip(times_rel, traj):
            day   = start_day + t_rel
            total = float(row.sum())
            burden_timeline.append((day, total))
            snap  = {int(N): float(c) for N, c in zip(Ns, row) if c > 0}
            ploidy_snapshots.append((day, snap))

            # ── MIN_SIZE check ──────────────────────────────────────────────
            if total < MIN_SIZE:
                extinction_day = day
                print(
                    f"\n  *** Tumor burden dropped below MIN_SIZE "
                    f"({MIN_SIZE:.2e}) on day {extinction_day:.2f} "
                    f"(burden = {total:.3e} cells). Halting simulation. ***\n"
                )
                path_info.append((drug, start_day, extinction_day))
                outer_done = True
                break  # exit inner timepoint loop

        if not outer_done:
            path_info.append((drug, start_day, end_day))
            ploidy = new_ploidy

    return burden_timeline, ploidy_snapshots, path_info, extinction_day


# =============================================================================
#  PLOIDY BINNING  (chromosomes / HAPLOID_N, rounded to nearest 0.1)
# =============================================================================

def _build_ploidy_bins(ploidy_snapshots):
    """Convert per-chromosome snapshots into binned count/fraction matrices.

    Returns
    -------
    sorted_bins  : list of float bin labels (e.g. 2.0, 3.0, 4.0)
    snap_days    : ndarray of day values
    count_matrix : ndarray (n_bins, n_snaps) — raw cell counts per bin
    frac_matrix  : ndarray (n_bins, n_snaps) — fraction of total per bin
    """
    all_bins: set[float] = set()
    binned_data = []   # list of (day, {bin: count}, {bin: frac})

    for day, ploidy_dict in ploidy_snapshots:
        total = sum(ploidy_dict.values())
        if total <= 0:
            binned_data.append((day, {}, {}))
            continue
        counts: dict[float, float] = {}
        for chr_count, cell_count in ploidy_dict.items():
            ratio = chr_count / HAPLOID_N
            b     = round(round(ratio * 10) / 10, 1)
            counts[b] = counts.get(b, 0.0) + cell_count
        fracs = {b: v / total for b, v in counts.items()}
        binned_data.append((day, counts, fracs))
        all_bins.update(counts.keys())

    sorted_bins  = sorted(all_bins)
    snap_days    = np.array([d[0] for d in binned_data])
    count_matrix = np.zeros((len(sorted_bins), len(binned_data)))
    frac_matrix  = np.zeros((len(sorted_bins), len(binned_data)))

    for col_i, (_, counts, fracs) in enumerate(binned_data):
        for row_i, b in enumerate(sorted_bins):
            count_matrix[row_i, col_i] = counts.get(b, 0.0)
            frac_matrix[row_i,  col_i] = fracs.get(b,  0.0)

    return sorted_bins, snap_days, count_matrix, frac_matrix


# =============================================================================
#  PLOTTING
# =============================================================================

def _add_drug_shading(ax, path_info, add_legend=False):
    """Shade each treatment cycle on *ax*."""
    shaded: set[str] = set()
    for drug, start_day, end_day in path_info:
        color = DRUG_COLORS.get(drug, "lightgrey")
        label = drug if (add_legend and drug not in shaded) else None
        ax.axvspan(start_day, end_day, color=color, alpha=0.15,
                   linewidth=0, label=label)
        shaded.add(drug)


def plot_results(burden_timeline, ploidy_snapshots, path_info):
    burden_days   = np.array([t for t, _ in burden_timeline])
    burden_values = np.array([b for _, b in burden_timeline])

    sorted_bins, snap_days, count_matrix, frac_matrix = _build_ploidy_bins(
        ploidy_snapshots
    )

    cmap = plt.get_cmap("tab20")

    # ── Plot 1: Total tumor burden ────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    _add_drug_shading(ax1, path_info, add_legend=True)
    ax1.plot(burden_days, burden_values,
             color="steelblue", linewidth=1.8, label="Predicted (total)")
    ax1.set_xlabel("Day")
    ax1.set_ylabel("Number of cells")
    ax1.set_title("Tumor Burden — Total")
    ax1.legend(title="Treatment", loc="upper left", bbox_to_anchor=(1, 1))
    ax1.set_yscale("log")
    ax1.grid(True, alpha=0.3)
    fig1.tight_layout()
    fig1.savefig("plot1_tumor_burden.png", dpi=150)
    print("  Saved: plot1_tumor_burden.png")

    # ── Plot 1b: Total burden + per-ploidy cell counts (same axis) ────────────
    fig1b, ax1b = plt.subplots(figsize=(12, 6))
    _add_drug_shading(ax1b, path_info)
    ax1b.plot(burden_days, burden_values,
              color="steelblue", linewidth=2.5, zorder=3,
              label="Total burden (predicted)")
    for row_i, b in enumerate(sorted_bins):
        counts = count_matrix[row_i, :]
        if counts.max() > 1e-6:
            ax1b.plot(snap_days, counts,
                      color=cmap(row_i % 20),
                      linewidth=2.0,
                      label=f"Ploidy {b:.1f}×")
    ax1b.set_xlabel("Day")
    ax1b.set_ylabel("Number of cells")
    ax1b.legend(title="Treatment / Ploidy", loc="upper left",
                fontsize=7, ncol=2, bbox_to_anchor=(1, 1))
    ax1b.set_yscale("log")
    ax1b.grid(True, alpha=0.3)
    fig1b.suptitle(
        f"Tumor Burden & Ploidy Cell Counts\n"
        f"(ploidy = chromosomes / {HAPLOID_N}, rounded to 0.1)",
        fontsize=10,
    )
    fig1b.tight_layout()
    fig1b.savefig("plot1b_burden_and_ploidy_counts.png", dpi=150)
    print("  Saved: plot1b_burden_and_ploidy_counts.png")

    # ── Plot 2: Ploidy fraction distribution over time ────────────────────────
    fig2, ax2 = plt.subplots(figsize=(11, 6))
    _add_drug_shading(ax2, path_info)
    for row_i, b in enumerate(sorted_bins):
        fracs = frac_matrix[row_i, :]
        if fracs.max() > 1e-6:
            ax2.plot(snap_days, fracs,
                     color=cmap(row_i % 20),
                     linewidth=2.5,
                     label=f"Ploidy {b:.1f}×")
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Fraction of cells")
    ax2.set_title(
        f"Ploidy Distribution Over Time\n"
        f"(chromosomes / {HAPLOID_N}, rounded to nearest 0.1)"
    )
    ax2.legend(title="Treatment / Ploidy", loc="upper right",
               fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig("plot2_ploidy_distribution.png", dpi=150)
    print("  Saved: plot2_ploidy_distribution.png")

    plt.show()


# =============================================================================
#  ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    print("Running forward simulation...")
    burden_timeline, ploidy_snapshots, path_info, extinction_day = run_simulation()
    final_burden = burden_timeline[-1][1]

    if extinction_day is not None:
        print(f"Simulation halted: burden crossed MIN_SIZE ({MIN_SIZE:.2e}) "
              f"on day {extinction_day:.2f}  "
              f"(final recorded burden = {final_burden:.3e} cells)")
    else:
        print(f"Final total burden: {final_burden:.3e} cells")

    plot_results(burden_timeline, ploidy_snapshots, path_info)