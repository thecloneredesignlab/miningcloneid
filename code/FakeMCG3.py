from math import sqrt
import numpy as np
import os
import matplotlib.pyplot as plt
import random
from Missegregation_Model import simulate_karyotype_ode_piecewise, get_concentration_curve, build_times_with_doses, f
from time import time

start = time()


class Node:
    def __init__(self, ploidy_status, cycle, parent=None):
        self.ploidy_status = ploidy_status
        self.cycle = cycle
        self.parent = parent
        self.children = {}
        self.N = 0
        self.W = 0.0  # reward

    def is_terminal(self):
        total = sum(self.ploidy_status.values())
        return self.cycle >= total_cycles or total < min_size or total > max_size

    def is_fully_expanded(self):
        return len(self.children) == len(drugs)


def simulate_next_state(ploidy_status, drug):
    C_fn = get_concentration_curve(drug)
    TIMES = build_times_with_doses((0.0, d_switch), 0.05, drug_name=drug, include_days=True, eps=1e-8)
    t, Ns, T_mat, T_total, M = simulate_karyotype_ode_piecewise(
        ploidy_status, drug, t_span=(0.0, d_switch), r=0.3, Kcap=1e12, beta=0.05,
        N_min=10, N_max=90,
        C_fn=C_fn,
        f_param_fn=f,
        t_eval=TIMES,
        max_step=0.05,
        renormalize_M=True
    )
    final_counts = T_mat[:, -1]
    new_status = {int(N): float(count) for N, count in zip(Ns, final_counts) if count > 0}
    trajectory = T_mat.T[1:]  # Skip first point (t=0), shape: (n_times-1, n_chromosomes)
    confidence = 1.0

    # Return Ns (chromosome numbers) so we can calculate average ploidy later
    return new_status, trajectory, confidence, Ns


def select_best_child_to_explore(node, c):
    best_score, best_child = -float('inf'), None
    for drug, child in node.children.items():
        Q = child.W / (child.N + 1e-6)
        U = c * np.sqrt(np.log(node.N + 1e-6) / (child.N + 1e-6))
        score = Q + U
        if score > best_score:
            best_score, best_child = score, child
    return best_child


def expand(node):
    untried = [d for d in drugs if d not in node.children]
    drug = random.choice(untried)
    new_ploidy, _, _, _ = simulate_next_state(node.ploidy_status, drug)
    child = Node(new_ploidy, node.cycle + 1, parent=node)
    node.children[drug] = child
    return child


def rollout(node, rollout_depth):
    ploidy = dict(node.ploidy_status)
    outcome = "normal"
    confidence = 1.0

    #parameters for reward calculation
    alpha = 0.01
    p_order = 3
    safe_size = 1e10
    beta_bonus = 1.0

    path_burdens = [sum(ploidy.values())]

    for step in range(rollout_depth):
        total = sum(ploidy.values())
        if total < min_size:
            break
        elif total > max_size:
            break
        drug = random.choice(drugs)
        ploidy, _, confidence, _ = simulate_next_state(ploidy, drug)
        path_burdens.append(sum(ploidy.values()))

    final_burden = sum(ploidy.values())

    # --- Reward logic ---
    # old: reward = confidence * ( (2e10 - final_burden) / 2e10)
    reward = 0.0
    for burden in path_burdens:
        reward -= (burden / max_size)
        if burden > safe_size:
            threshold_penalty = ((burden - safe_size) / (max_size - safe_size)) ** p_order
            reward -= alpha * threshold_penalty
    reward = reward / len(path_burdens)

    final_burden = path_burdens[-1]
    if final_burden < min_size:
        steps_to_extinct = next(i for i, b in enumerate(path_burdens) if b < min_size)
        bonus = beta_bonus * ((rollout_depth - steps_to_extinct) / max(1, rollout_depth))
        reward += bonus

    return reward * confidence


def backpropagate(node, reward):
    while node is not None:
        node.N += 1
        node.W += reward
        node = node.parent


# ---- Main loop ----
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
d_switch = 7
total_cycles = 12
min_size = 1e5
max_size = 2e10
depth = 30
num_rollouts = 100
c = sqrt(2)

cycle_counter = 0

# Initial ploidy status: map chromosome numbers to cell counts
ploidy_status = {46: 1.5e9, 69: 0.3e9, 92: 0.25e9}

# Store total tumor burden at each time point
tumor_burden_times = [sum(ploidy_status.values())]

# Calculate and store initial average ploidy
initial_chroms = np.array(list(ploidy_status.keys()))
initial_counts = np.array(list(ploidy_status.values()))
initial_avg_ploidy = np.sum(initial_chroms * initial_counts) / np.sum(initial_counts)
avg_ploidy_times = [initial_avg_ploidy]

best_drug_list = []

for decision in range(total_cycles):
    # Initialize a FRESH root node for the current state
    root = Node(ploidy_status, decision)

    # --- Run MCTS rollouts ---
    for _ in range(num_rollouts):
        node = root
        # Selection
        while node.is_fully_expanded() and not node.is_terminal():
            node = select_best_child_to_explore(node, c)
        # Expansion
        if not node.is_terminal():
            child = expand(node)
            reward = rollout(child, depth)
            backpropagate(child, reward)
        else:
            reward = rollout(node, 0)
            backpropagate(node, reward)


    def get_q(child):
        return child.W / (child.N + 1e-6)


    # Pick best drug
    best_drug = max(root.children.items(), key=lambda kv: get_q(kv[1]))[0]
    best_drug_list.append(best_drug)

    print(f"Cycle {decision + 1}: best drug is {best_drug} with tumor burden {sum(ploidy_status.values()):.2e} cells")

    cycle_counter += 1
    # Capture Ns (chromosome array) to calculate averages
    ploidy_status, trajectory, _, Ns = simulate_next_state(ploidy_status, best_drug)

    # Add trajectory points to tracking lists
    Ns = np.array(Ns)  # Ensure numpy array for vector math
    for time_point in trajectory:
        # time_point is an array of cell counts for each chromosome index
        current_total = np.sum(time_point)
        tumor_burden_times.append(current_total)

        # Calculate weighted average ploidy
        if current_total > 0:
            avg_p = np.sum(time_point * Ns) / current_total
        else:
            avg_p = 0
        avg_ploidy_times.append(avg_p)

    # Check for extinction or overflow
    if sum(ploidy_status.values()) < min_size:
        print("Tumor extinction detected.")
        break
    elif sum(ploidy_status.values()) > max_size:
        print("Tumor burden exceeded max size.")
        break

print(f"Final tumor burden: {sum(ploidy_status.values()):.2e}")
print(f"Total cycles: {cycle_counter}")
print(f"Best drug sequence: {best_drug_list}")
print(f"Final ploidy distribution: {ploidy_status}")

# --- Plotting ---
TB = np.array(tumor_burden_times)
AP = np.array(avg_ploidy_times)  # Average Ploidy array

drug_colors = {
    "gemcitabine": "black",
    "bay1895344": "red",
    "alisertib": "green",
    "ispinesib": "blue",
    "none": "yellow"
}

# Calculate time vector
points_per_cycle = int(d_switch / 0.05)
time_vec = np.arange(len(TB)) * 0.05

# --------------------------
# Plot 1: Tumor Burden
# --------------------------
plt.figure(figsize=(10, 6))

start_idx = 0
for i, drug in enumerate(best_drug_list):
    end_idx = min(start_idx + points_per_cycle, len(TB) - 1)

    # Color logic
    c = drug_colors.get(drug, "gray")
    lbl = drug if i == 0 or drug not in best_drug_list[:i] else None

    plt.plot(
        time_vec[start_idx:end_idx + 1],
        TB[start_idx:end_idx + 1],
        color=c,
        linewidth=2,
        label=lbl
    )
    start_idx = end_idx

plt.ylabel("Total Tumor Volume (Cell Count)", fontsize=12)
plt.xlabel("Time (Days)", fontsize=12)
plt.yscale('log')
plt.title("Total Tumor Burden Over Time", fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(title="Drug Applied", fontsize=10)
plt.tight_layout()
plt.show()  # Show first plot

# --------------------------
# Plot 2: Average Ploidy
# --------------------------
plt.figure(figsize=(10, 6))

start_idx = 0  # Reset start index for the second plot
for i, drug in enumerate(best_drug_list):
    end_idx = min(start_idx + points_per_cycle, len(AP) - 1)

    c = drug_colors.get(drug, "gray")

    # Logic to only label the first occurrence of each drug
    lbl = drug if i == 0 or drug not in best_drug_list[:i] else None

    plt.plot(
        time_vec[start_idx:end_idx + 1],
        AP[start_idx:end_idx + 1]/23,
        color=c,
        linewidth=2,
        label=lbl
    )
    start_idx = end_idx

plt.xlabel("Time (Days)", fontsize=12)
plt.ylabel("Average Ploidy", fontsize=12)
plt.title("Average Tumor Ploidy Over Time", fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(title="Drug Applied", fontsize=10)
plt.tight_layout()
plt.show()  # Show second plot

end = time()
print(f"\nTotal computation time: {end - start:.2f} seconds")