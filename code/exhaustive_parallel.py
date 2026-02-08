import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from time import time
# Assuming PujanEarlyVersionModel is in your path
from PujanEarlyVersionModel import ploidy_forcast

# --- Configuration & Hyperparameters ---
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
MAX_DEPTH = 12 # Warning: Exhaustive search is O(len(drugs)^MAX_DEPTH)
# At a depth of 13, extinction will be reached
# Question is if extinction can be reached earlier; only able to validate with exhaustive search at depth 12
# Should take ~10 days on the cluster with 6 CPUs, ~32GB memory
min_size = 1e5
max_size = 2e10
default_len = 7.0
cycle_lengths = {
    "volasertib": default_len, "alisertib": 21.0, "cytarabine": default_len,
    "gemcitabine": 28.0, "ispinesib": default_len, "none": 1.0
}


def get_cycle_length(drug):
    return cycle_lengths.get(drug, default_len)


def simulate_next_state(ploidy_status, drug):
    # This remains the same as your original code
    ploidies, t_ode, T_mat_ode, t_sde, Tpaths = ploidy_forcast(
        ploidy_status, drug, get_cycle_length(drug), N_SIMS=1000, R_BASE=0.575
    )
    final_per_ploidy = Tpaths[:, :, -1]
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)
    mean_trajectory = np.mean(Tpaths, axis=0).T
    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k]) for k in range(len(ploidies))}
    return new_status, mean_trajectory[1:]


def exhaustive_search(current_ploidy, depth, path_so_far, trajectory_so_far):
    current_burden = sum(current_ploidy.values())
    if current_burden < min_size:
        return path_so_far, trajectory_so_far, True
    if depth >= MAX_DEPTH:
        return path_so_far, trajectory_so_far, False
    if current_burden > max_size:
        return None, None, False

    best_final_path = None
    best_final_trajectory = None

    for drug in drugs:
        next_ploidy, segment_trajectory = simulate_next_state(current_ploidy, drug)
        res_path, res_traj, success = exhaustive_search(
            next_ploidy, depth + 1, path_so_far + [drug], trajectory_so_far + list(segment_trajectory)
        )
        if success:
            return res_path, res_traj, True
        if res_path is not None:
            if best_final_path is None or sum(next_ploidy.values()) < sum(
                    best_final_path[-1] if isinstance(best_final_path[-1], dict) else {}):
                best_final_path, best_final_trajectory = res_path, res_traj
    return best_final_path, best_final_trajectory, False


# --- Parallel Wrapper ---
def run_branch(initial_drug, initial_ploidy):
    """Worker function to process one top-level drug branch."""
    print(f"Starting branch: {initial_drug}")
    next_ploidy, segment_trajectory = simulate_next_state(initial_ploidy, initial_drug)
    return exhaustive_search(
        next_ploidy, 1, [initial_drug], [np.array(list(initial_ploidy.values()))] + list(segment_trajectory)
    )


if __name__ == "__main__":
    start_time = time()
    initial_ploidy = {2.0: 1.5 * 1e9, 3.0: 0, 4.0: 0.55 * 1e9}

    print(f"Starting Parallel Exhaustive Search (CPUs: {os.cpu_count()})...")

    winning_path = None
    winning_trajectory = None
    is_extinct = False

    # Using ProcessPoolExecutor to run branches in parallel
    with ProcessPoolExecutor() as executor:
        # Submit each initial drug as a separate task
        futures = {executor.submit(run_branch, drug, initial_ploidy): drug for drug in drugs}

        for future in as_completed(futures):
            res_path, res_traj, success = future.result()
            if success:
                # If we found extinction, we can stop and take this result
                winning_path, winning_trajectory, is_extinct = res_path, res_traj, True
                print(f"Extinction found by branch: {futures[future]}")
                executor.shutdown(wait=False, cancel_futures=True)  # Kill other branches
                break
            else:
                # Keep track of the best result if no extinction yet
                if winning_path is None:
                    winning_path, winning_trajectory = res_path, res_traj

    end_time = time()

    # --- Plotting (Only in main process) ---
    if winning_path:
        print(f"\nSearch Complete in {end_time - start_time:.2f}s")
        if winning_path:
            print(f"\nSearch Complete in {end_time - start_time:.2f}s")
            print(f"Extinction Achieved: {is_extinct}")
            print(f"Best Strategy: {winning_path}")

            # Process trajectory for plotting
            ploidy_over_time = np.array(winning_trajectory)
            total_burden = np.sum(ploidy_over_time, axis=1)
            time_vec = np.arange(len(total_burden)) * 0.1

            drug_colors = {
                "gemcitabine": "orange", "bay1895344": "red",
                "alisertib": "green", "ispinesib": "blue", "none": "yellow"
            }

            plt.figure(figsize=(12, 6))
            plt.yscale('log')
            plt.xlabel("Time (Days)")
            plt.ylabel("Cell Count")
            plt.title(f"Exhaustive Search Results (Extinction: {is_extinct})")

            # Shading segments
            current_t = 0
            for drug in winning_path:
                duration = get_cycle_length(drug)
                plt.axvspan(current_t, current_t + duration, color=drug_colors.get(drug, "gray"), alpha=0.15,
                            label=drug)
                current_t += duration

            # Plot lines
            plt.plot(time_vec, total_burden, label="Total Burden", color="black", linewidth=3)
            for i, p_label in enumerate([2.0, 3.0, 4.0]):
                plt.plot(time_vec, ploidy_over_time[:, i], label=f"{int(p_label)}n", alpha=0.7)

            # Clean up legend (remove duplicates)
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), loc='upper left', bbox_to_anchor=(1, 1))

            plt.tight_layout()
            plt.show()
        else:
            print("No viable paths found (all paths exceeded max_size or search space exhausted).")
        plt.show()