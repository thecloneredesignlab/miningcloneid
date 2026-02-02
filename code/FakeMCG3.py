# This code uses the original model from Pujan's early version and implements a Monte Carlo Tree Search (MCTS)
# This code will output the ploidy compositions of the initial ploidy populations (ie 2n, 3n, 4n) at each time point
from math import sqrt
import numpy as np
import os
import matplotlib.pyplot as plt
import random
from PujanEarlyVersionModel import ploidy_forcast
from time import time

start = time()

class Node:
    def __init__(self, ploidy_status, cycle, parent=None):
        self.ploidy_status = ploidy_status
        self.cycle = cycle
        self.parent = parent
        self.children = {}
        self.N = 0
        self.W = 0.0  # reward (negative tumor burden)

    def is_terminal(self):
        total = sum(self.ploidy_status.values())
        return self.cycle >= total_cycles or total < min_size or total > max_size

    def is_fully_expanded(self):
        return len(self.children) == len(drugs)

def simulate_next_state(ploidy_status, drug):
    ploidy_cell_count = ploidy_status
    ploidies, t_ode, T_mat_ode, t_sde, Tpaths = ploidy_forcast(ploidy_cell_count, drug, d_switch, N_SIMS=1000, R_BASE=0.15)
    final_total_burdens = Tpaths[:, :, -1].sum(axis=1)

    # Freedman-Diaconis Rule for Binning
    n = len(final_total_burdens)
    if n > 1:
        q75, q25 = np.percentile(final_total_burdens, [75, 25])
        iqr = q75 - q25
        # Calculate optimal bin width: 2 * IQR * n^(-1/3)
        if iqr > 0:
            bin_width = 2 * iqr * (n ** (-1 / 3))
            data_range = final_total_burdens.max() - final_total_burdens.min()
            num_bins = int(np.ceil(data_range / bin_width))
        else:
            num_bins = 1  # Data has zero variance
    else:
        num_bins = 1

    # Calculate Shannon Entropy
    counts, _ = np.histogram(final_total_burdens, bins=num_bins)
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    entropy = -np.sum(probs * np.log(probs))
    max_entropy = np.log(num_bins)
    if max_entropy > 1e-9:  # Avoid division by zero if num_bins=1
        uncertainty_score = entropy / max_entropy
    else:
        uncertainty_score = 0.0


    final_per_ploidy = Tpaths[:, :, -1]  # shape: (num_paths, num_ploidies)
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)

    mean_trajectory = np.mean(Tpaths, axis=0)
    y = mean_trajectory.T

    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k]) for k in range(len(ploidies))}
    confidence = 1 - uncertainty_score
    return new_status, y[1:], confidence

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
    new_ploidy, _, _ = simulate_next_state(node.ploidy_status, drug)
    child = Node(new_ploidy, node.cycle + 1, parent=node)
    node.children[drug] = child
    return child

def rollout(node, rollout_depth):
    ploidy = dict(node.ploidy_status)
    outcome = "normal"
    confidence = 1.0

    # parameters for reward calculation
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
        ploidy, _, confidence = simulate_next_state(ploidy, drug)
        path_burdens.append(sum(ploidy.values()))

    final_burden = sum(ploidy.values())

    # --- Reward logic ---
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
num_rollouts = 500
c = sqrt(2)

cycle_counter = 0
ploidy_status = {2.0: 1.5*1e9, 3.0: 0.3*1e9, 4.0: 0.25*1e9}
tumor_burden_times = [np.array(list(ploidy_status.values()), dtype=float)]
best_drug_list = []

for decision in range(total_cycles):
    # Initialize a FRESH root node for the current state
    root = Node(ploidy_status, decision)

    # --- Run MCTS rollouts ---
    for _ in range(num_rollouts):
        node = root
        # Selection
        while node.is_fully_expanded() and not node.is_terminal(): # Traverse down tree acording to UCB1
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
    ploidy_status, ploidies, _ = simulate_next_state(ploidy_status, best_drug)
    tumor_burden_times.extend(ploidies)

    temp = np.array([np.sum(arr) for arr in tumor_burden_times])
    # print(min(temp), max(temp))

    if min(temp) < min_size:
        print("Tumor extinction detected in trajectory.")
        # Delete values in tumor_burden_times after extinction
        extinction_index = next(i for i, v in enumerate(temp) if v < min_size)
        tumor_burden_times = tumor_burden_times[:extinction_index + 1]
        break
    elif max(temp) > max_size:
        print("Tumor burden exceeded max size in trajectory.")
        # Delete values in tumor_burden_times after exceeding max size
        exceed_index = next(i for i, v in enumerate(temp) if v > max_size)
        tumor_burden_times = tumor_burden_times[:exceed_index + 1]
        break

print(f"Final tumor burden: {sum(ploidy_status.values())}")
print(cycle_counter)
print(best_drug_list)
print(ploidy_status)

# --- Plotting ---
TB = np.array([np.sum(arr) for arr in tumor_burden_times])
colors = ['black','red', 'green', 'blue', 'yellow']
plt.figure(figsize=(8, 6))
plt.xlabel("Time (Days)")
plt.ylabel("Total Tumor Volume (Cell Count)")
plt.yscale('log')

drug_colors = {
    "gemcitabine": "black",
    "bay1895344": "red",
    "alisertib": "green",
    "ispinesib": "blue",
    "none": "yellow"
}

time_vec = np.arange(len(TB)) * 0.1
start_idx = 0
for i, drug in enumerate(best_drug_list):
    end_idx = min(start_idx + (10*d_switch), len(TB) - 1)
    plt.plot(
        time_vec[start_idx:end_idx + 1],
        TB[start_idx:end_idx + 1],
        color=drug_colors.get(drug, "gray"),
        linewidth=2,
        label=drug if i == 0 or drug not in best_drug_list[:i] else None
    )
    start_idx = end_idx

plt.legend(title="Drug Applied")
plt.tight_layout()
plt.show()

plt.xlabel("Time (Days)")
plt.ylabel("Ploidy Tumor Volume (Cell Count)")
# plt.yscale('log')

ploidy_labels = [2.0, 3.0, 4.0]
ploidy_over_time = np.array(tumor_burden_times)  # shape: (time_steps, 3)
time_vec = np.arange(len(ploidy_over_time)) * 0.1
for i, label in enumerate(ploidy_labels):
    plt.plot(
        time_vec,
        ploidy_over_time[:, i],
        label=f"{int(label)}n",
        linewidth=2
    )

plt.legend(title="Ploidies")
plt.tight_layout()
plt.show()

end = time()
print(f"Total time taken: {end - start:.2f} seconds")