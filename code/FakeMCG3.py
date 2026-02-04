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

default_len = 7.0
cycle_lengths = {
    "volasertib": default_len, "alisertib": 21.0, "cytarabine": default_len,
    "gemcitabine": 28.0, "ispinesib": default_len, "umi-77": default_len,
    "navitoclax": default_len, "bay1895344": default_len, "Topotecan": 28.0,
    "Doxorubicin": 21.0, "osi-027": default_len, "abt-199": default_len,
    "abt-263": default_len, "ceralasertib": default_len, "adavosertib": default_len,
    "tas": default_len, "tegafur": default_len, "capecitabine": default_len,
    "5-azacytidine": default_len, "Elimusertib": 7.0, "none": 1.0
}

def get_cycle_length(drug):
    return cycle_lengths.get(drug, default_len)

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
        return self.cycle >= depth# or total < min_size or total > max_size

    def is_fully_expanded(self):
        return len(self.children) == len(drugs)

def simulate_next_state(ploidy_status, drug):
    ploidy_cell_count = ploidy_status
    ploidies, t_ode, T_mat_ode, t_sde, Tpaths = ploidy_forcast(ploidy_cell_count, drug, get_cycle_length(drug), N_SIMS=1000, R_BASE=0.575)
    final_per_ploidy = Tpaths[:, :, -1]  # shape: (num_paths, num_ploidies)
    mean_sde_per_ploidy = np.mean(final_per_ploidy, axis=0)

    mean_trajectory = np.mean(Tpaths, axis=0)
    y = mean_trajectory.T

    new_status = {ploidies[k]: float(mean_sde_per_ploidy[k]) for k in range(len(ploidies))}
    confidence = 1
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
    beta_bonus = 100.0

    path_burdens = [sum(ploidy.values())]

    for step in range(rollout_depth):
        total = sum(ploidy.values())
        # if total < min_size:
        #     break
        # elif total > max_size:
        #     break
        drug = random.choice(drugs)
        ploidy, _, confidence = simulate_next_state(ploidy, drug)
        path_burdens.append(sum(ploidy.values()))

    final_burden = sum(ploidy.values())

    # This reward method needs a robust version of confidence to work well
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

    # return (1 - (final_burden / max_size))

def backpropagate(node, reward):
    while node is not None:
        node.N += 1
        node.W += reward
        node = node.parent

# ---- Main loop ----
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
d_switch = 7
total_cycles = 10
min_size = 1e5
max_size = 2e10
depth = 30
num_rollouts = 10000
c = sqrt(2)

cycle_counter = 0
ploidy_status = {2.0: 1.5*1e9, 3.0: 0, 4.0: 0.55*1e9}
tumor_burden_times = [np.array(list(ploidy_status.values()), dtype=float)]
best_drug_list = []
actual_cycle_lengths = []

rigged_drugs = ["bay1895344", "ispinesib", "ispinesib"] * 100

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
    # best_drug = max(root.children.items(), key=lambda kv: get_q(kv[1]))[0]
    best_drug = rigged_drugs[decision]
    best_drug_list.append(best_drug)

    current_len = get_cycle_length(best_drug)
    actual_cycle_lengths.append(current_len)

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
ploidy_over_time = np.array(tumor_burden_times)
TB = np.sum(ploidy_over_time, axis=1)
time_vec = np.arange(len(TB)) * 0.1

drug_colors = {
    "gemcitabine": "orange", "bay1895344": "red",
    "alisertib": "green", "ispinesib": "blue", "none": "yellow"
}

plt.figure(figsize=(10, 6))
plt.yscale('log')
plt.xlabel("Time (Days)")
plt.ylabel("Ploidy Tumor Volume (Cell Count)")

# Dynamic Shading and Plotting
start_idx = 0
shaded_labels = set()

for i, drug in enumerate(best_drug_list):
    if start_idx >= len(time_vec): break

    # The length of this segment in time steps (0.1 day increments)
    seg_len = int(actual_cycle_lengths[i] * 10)
    end_idx = min(start_idx + seg_len, len(time_vec) - 1)

    plt.axvspan(
        time_vec[start_idx], time_vec[end_idx],
        color=drug_colors.get(drug, "gray"),
        alpha=0.15,
        label=drug if drug not in shaded_labels else None
    )
    shaded_labels.add(drug)
    start_idx = end_idx

# Plot Total and Individual Ploidies
plt.plot(time_vec, TB, label="Total", color="black", linewidth=2, alpha=0.6)
ploidy_labels = [2.0, 3.0, 4.0]
for i, label in enumerate(ploidy_labels):
    plt.plot(time_vec, ploidy_over_time[:, i], label=f"{int(label)}n", linewidth=2)

plt.legend(title="Treatment Strategy", loc='upper left', bbox_to_anchor=(1, 1))
plt.tight_layout()
plt.show()

end = time()
print(f"Total time taken: {end - start:.2f} seconds")