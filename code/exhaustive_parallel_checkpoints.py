import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import json
import hashlib
from time import time
# Assuming PujanEarlyVersionModel is in your path
from PujanEarlyVersionModel import ploidy_forcast

# --- Configuration ---
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
MAX_DEPTH = 8
min_size = 1e5
max_size = 2e10
default_len = 7.0
cycle_lengths = {
    "volasertib": default_len, "alisertib": 21.0, "cytarabine": default_len,
    "gemcitabine": 28.0, "ispinesib": default_len, "none": 1.0
}

CACHE_DIR = "sim_cache"
if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR)


# --- Persistence Layer Fixes ---

def get_cache_key(ploidy_status, drug):
    """Creates a unique hash. Handles both float and string keys from JSON."""
    # FIX: Convert k to float before rounding because JSON loads keys as strings
    sorted_ploidy = sorted([(str(round(float(k), 1)), round(v, 2)) for k, v in ploidy_status.items()])
    key_string = f"{sorted_ploidy}_{drug}"
    return hashlib.md5(key_string.encode()).hexdigest()


def load_local_cache(branch_name):
    path = os.path.join(CACHE_DIR, f"cache_{branch_name}.json")
    if os.path.exists(path):
        with open(path, "r") as f:
            return json.load(f)
    return {}


def save_local_cache(branch_name, data):
    path = os.path.join(CACHE_DIR, f"cache_{branch_name}.json")
    with open(path, "w") as f:
        json.dump(data, f)


# --- Core Logic ---

def simulate_next_state(ploidy_status, drug, branch_cache=None, branch_name=None):
    cache_key = get_cache_key(ploidy_status, drug)

    # Check cache first
    if branch_cache is not None and cache_key in branch_cache:
        cached = branch_cache[cache_key]
        # FIX: Convert the string keys from JSON back to floats for the model
        new_status = {float(k): v for k, v in cached['new_status'].items()}
        return new_status, [np.array(x) for x in cached['trajectory']]

    # Run actual simulation if not in cache
    ploidies, t_ode, T_mat_ode, t_sde, Tpaths = ploidy_forcast(
        ploidy_status, drug, get_cycle_length(drug), N_SIMS=1000, R_BASE=0.575
    )

    final_per_ploidy = Tpaths[:, :, -1]
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)
    mean_trajectory = np.mean(Tpaths, axis=0).T

    new_status = {float(ploidies[k]): float(mean_sde_per_ploidy[k]) for k in range(len(ploidies))}
    trajectory_list = [x.tolist() for x in mean_trajectory[1:]]

    if branch_cache is not None:
        branch_cache[cache_key] = {
            'new_status': new_status,
            'trajectory': trajectory_list
        }
        save_local_cache(branch_name, branch_cache)

    return new_status, mean_trajectory[1:]


def exhaustive_search(current_ploidy, depth, path_so_far, trajectory_so_far, branch_cache, branch_name):
    current_burden = sum(current_ploidy.values())

    if depth <= 4:
        print(f"[{branch_name}] Depth {depth}: {' -> '.join(path_so_far)} | Burden: {current_burden:.2e}")

    if current_burden < min_size:
        return path_so_far, trajectory_so_far, True
    if depth >= MAX_DEPTH:
        return path_so_far, trajectory_so_far, False
    if current_burden > max_size:
        return None, None, False

    best_final_path = None
    best_final_trajectory = None

    for drug in drugs:
        next_ploidy, segment_trajectory = simulate_next_state(current_ploidy, drug, branch_cache, branch_name)
        res_path, res_traj, success = exhaustive_search(
            next_ploidy, depth + 1, path_so_far + [drug],
                         trajectory_so_far + list(segment_trajectory), branch_cache, branch_name
        )
        if success:
            return res_path, res_traj, True

        # Keep track of the path that resulted in the lowest burden if no extinction found
        if res_path is not None:
            if best_final_path is None or sum(next_ploidy.values()) < sum(
                    best_final_path[-1] if isinstance(best_final_path[-1], dict) else {}):
                best_final_path, best_final_trajectory = res_path, res_traj

    return best_final_path, best_final_trajectory, False


def run_branch(initial_drug, initial_ploidy):
    branch_name = initial_drug
    branch_cache = load_local_cache(branch_name)

    print(f"\n>>> INITIATING/RESUMING BRANCH: {branch_name}")
    next_ploidy, segment_trajectory = simulate_next_state(initial_ploidy, initial_drug, branch_cache, branch_name)

    return exhaustive_search(
        next_ploidy, 1, [initial_drug],
        [np.array(list(initial_ploidy.values()))] + list(segment_trajectory),
        branch_cache, branch_name
    )


def get_cycle_length(drug):
    return cycle_lengths.get(drug, default_len)


# --- Main Execution ---
if __name__ == "__main__":
    start_time = time()
    initial_ploidy = {2.0: 1.5 * 1e9, 3.0: 0, 4.0: 0.55 * 1e9}
    print(f"Starting Parallel Exhaustive Search with Persistent Sub-branch Caching...")

    winning_path = None
    winning_trajectory = None
    is_extinct = False

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(run_branch, drug, initial_ploidy): drug for drug in drugs}

        try:
            for future in as_completed(futures):
                res_path, res_traj, success = future.result()
                if success:
                    winning_path, winning_trajectory, is_extinct = res_path, res_traj, True
                    print(f"\nSUCCESS: Extinction found in branch {futures[future]}")
                    executor.shutdown(wait=False, cancel_futures=True)
                    break
                else:
                    if winning_path is None or (res_path and len(res_path) > len(winning_path)):
                        winning_path, winning_trajectory = res_path, res_traj
        except KeyboardInterrupt:
            print("\nExiting. All progress has been saved to the 'sim_cache' folder.")

    # --- Plotting ---
    if winning_path:
        if winning_path:
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