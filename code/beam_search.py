import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from time import time
from PujanEarlyVersionModel import ploidy_forcast

# --- Configuration & Hyperparameters ---
BEAM_WIDTH = 100
MAX_DEPTH = 100
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
min_size = 1e5
max_size = 2e10
default_len = 7.0
cycle_lengths = {
    "volasertib": default_len, "alisertib": 21.0, "cytarabine": default_len,
    "gemcitabine": 28.0, "ispinesib": default_len, "bay1895344": default_len,
    "none": default_len
}


def get_cycle_length(drug):
    return cycle_lengths.get(drug, default_len)


def simulate_next_state(ploidy_status, drug):
    """Core simulation function."""
    ploidies, t_ode, T_mat_ode, t_sde, Tpaths = ploidy_forcast(
        ploidy_status, drug, get_cycle_length(drug), N_SIMS=1000, R_BASE=0.25, K_CAP = 2e10
    )
    final_per_ploidy = Tpaths[:, :, -1]
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)
    # Get trajectory. Note: we keep all points including the overlap to ensure index continuity
    mean_trajectory = np.mean(Tpaths, axis=0).T
    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k]) for k in range(len(ploidies))}

    # Returning the trajectory starting from index 1 to avoid duplicating the 'initial' state point
    return new_status, mean_trajectory[1:]


def simulate_next_state_wrapper(ploidy, drug, path, traj):
    """Helper for parallel execution to keep track of metadata."""
    new_status, segment_trajectory = simulate_next_state(ploidy, drug)
    return new_status, segment_trajectory, path, traj, drug


def beam_search_step(current_beams, executor):
    """Processes one level of the tree using the provided executor."""
    next_candidates = []
    futures = []

    for burden, ploidy, path, traj, extinct in current_beams:
        for drug in drugs:
            futures.append(executor.submit(simulate_next_state_wrapper, ploidy, drug, path, traj))

    for future in as_completed(futures):
        next_ploidy, segment_traj, old_path, old_traj, drug = future.result()
        new_burden = sum(next_ploidy.values())

        # Track the drug AND the actual number of points in this segment
        segment_info = (drug, len(segment_traj))

        # Success condition: Extinction
        if new_burden < min_size:
            return [(new_burden, next_ploidy, old_path + [segment_info], old_traj + list(segment_traj), True)]

        # Viability check
        if new_burden <= max_size:
            next_candidates.append(
                (new_burden, next_ploidy, old_path + [segment_info], old_traj + list(segment_traj), False)
            )

    # Sort by lowest burden (greedy optimization) and truncate to BEAM_WIDTH
    next_candidates.sort(key=lambda x: x[0])
    return next_candidates[:BEAM_WIDTH]


if __name__ == "__main__":
    start_time = time()
    initial_ploidy = {2.0: 1.5 * 1e9, 3.0: 0.3 * 1e9, 4.0: 0.25 * 1e9}
    initial_burden = sum(initial_ploidy.values())

    # beam format: (burden, ploidy_dict, path_list, trajectory_list, is_extinct)
    beam = [(initial_burden, initial_ploidy, [], [np.array(list(initial_ploidy.values()))], False)]

    print(f"Starting Beam Search (Width: {BEAM_WIDTH}, Max Depth: {MAX_DEPTH})...")

    with ProcessPoolExecutor(max_workers=8) as executor:
        for d in range(MAX_DEPTH):
            print(f"Step {d + 1}/{MAX_DEPTH} | Current Best Burden: {beam[0][0]:.2e}")
            beam = beam_search_step(beam, executor)

            if not beam:
                print("All paths led to tumor escape (max_size exceeded).")
                break

            if beam[0][4]:  # is_extinct
                print(f"Extinction target reached at depth {d + 1}!")
                break

    drug_colors = {
        "gemcitabine": "orange", "bay1895344": "red",
        "alisertib": "green", "ispinesib": "blue", "none": "yellow"
    }
    # --- Corrected Plotting Logic ---
    if beam:
        best_result = beam[0]
        winning_path_info = best_result[2]  # List of (drug, segment_length)
        winning_trajectory = np.array(best_result[3])

        # 1. Calculate the actual time for every point in the trajectory
        # We start at time 0.0
        time_points = [0.0]
        current_time = 0.0

        for drug_name, seg_len in winning_path_info:
            duration = get_cycle_length(drug_name)
            # Distribute the duration across the number of steps in this segment
            # segment_steps = np.linspace(current_time, current_time + duration, seg_len + 1)[1:]
            step_size = duration / seg_len
            segment_times = current_time + (np.arange(1, seg_len + 1) * step_size)
            time_points.extend(segment_times.tolist())
            current_time += duration

        time_vec = np.array(time_points)
        TB = np.sum(winning_trajectory, axis=1)

        # Ensure time_vec and TB match in length (handling floating point rounding)
        if len(time_vec) > len(TB):
            time_vec = time_vec[:len(TB)]

        # --- Standard Plotting Commands ---
        plt.figure(figsize=(12, 6))
        plt.xlabel("Time (Days)")
        plt.ylabel("Ploidy Tumor Volume (Cell Count)")
        plt.title(f"Beam Search Strategy (Depth: {len(winning_path_info)} N={BEAM_WIDTH})")

        # Dynamic Shading
        current_idx = 1
        current_time_mark = 0.0
        shaded_labels = set()

        for drug_name, seg_len in winning_path_info:
            duration = get_cycle_length(drug_name)
            end_time = current_time_mark + duration

            plt.axvspan(
                current_time_mark,
                end_time,
                color=drug_colors.get(drug_name, "gray"),
                linewidth=0,
                alpha=0.15,
                label=drug_name if drug_name not in shaded_labels else None
            )

            shaded_labels.add(drug_name)
            current_time_mark = end_time
            current_idx += seg_len

        # Plot Data
        plt.plot(time_vec, TB, label="Total", color="black", linewidth=3, alpha=0.7)
        ploidy_labels = [2.0, 3.0, 4.0]
        for i, label in enumerate(ploidy_labels):
            if i < winning_trajectory.shape[1]:
                plt.plot(time_vec, winning_trajectory[:, i], label=f"{int(label)}n", linewidth=2)

        plt.legend(title="Treatment Strategy", loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(True, which="both", ls="-", alpha=0.2)
        plt.tight_layout()
        plt.show()