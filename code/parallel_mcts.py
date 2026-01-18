from math import sqrt
import numpy as np
import os
import matplotlib.pyplot as plt
import random
from PujanEarlyVersionModel import ploidy_forcast
from time import time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

epsilon = 1e-30

class Node:
    def __init__(self, ploidy_status, cycle, parent=None):
        self.ploidy_status = ploidy_status
        self.cycle = cycle
        self.parent = parent
        self.children = {}
        self.N = 0
        self.W = 0.0  # reward (negative tumor burden)

    def is_terminal(self, total_cycles, min_size, max_size):
        total = sum(self.ploidy_status.values())
        return self.cycle >= total_cycles or total < min_size or total > max_size

    def is_fully_expanded(self, drugs):
        return len(self.children) == len(drugs)


def simulate_next_state(ploidy_status, drug, d_switch):
    ploidy_cell_count = ploidy_status
    ploidies, t_ode, T_mat_ode, _, _ = ploidy_forcast(ploidy_cell_count, drug, d_switch)
    states = ploidies
    y = T_mat_ode.T
    new_status = {states[k]: float(y[-1][k]) for k in range(len(states))}
    return new_status, y[1:]


def select_best_child(node, c, current_cycle, decay_factor):
    best_score, best_child = -float('inf'), None
    for drug, child in node.children.items():
        age = max(0, current_cycle - child.cycle)
        effective_factor = decay_factor ** age
        eff_W = child.W * effective_factor
        eff_N = child.N * effective_factor

        Q = eff_W / (eff_N + epsilon)
        U = c * np.sqrt(np.log(node.N + 1) / (eff_N + epsilon))
        score = Q + U

        if score > best_score:
            best_score, best_child = score, child
    return best_child


def expand(node, drugs, d_switch):
    untried = [d for d in drugs if d not in node.children]
    if not untried:
        return None
    drug = random.choice(untried)
    new_ploidy, _ = simulate_next_state(node.ploidy_status, drug, d_switch)
    child = Node(new_ploidy, node.cycle + 1, parent=node)
    node.children[drug] = child
    return child


def rollout(ploidy_status, rollout_depth, drugs, d_switch, min_size, max_size):
    ploidy = dict(ploidy_status)
    extinct = False
    maxed_out = False
    extinction_step = 0
    outcome = "normal"
    step = 0

    for step in range(rollout_depth):
        total = sum(ploidy.values())
        if total < min_size:
            extinct = True
            extinction_step = step
            outcome = "extinct"
            break
        elif total > max_size:
            maxed_out = True
            outcome = "maxout"
            break
        drug = random.choice(drugs)
        ploidy, _ = simulate_next_state(ploidy, drug, d_switch)

    final_burden = sum(ploidy.values())

    if extinct:
        extinction_boost = (rollout_depth - extinction_step) / rollout_depth
        reward = 0.9 + extinction_boost * 0.1
    elif maxed_out:
        reward = 0.01 * step
    else:
        reward = 0.9 - final_burden / 2e10

    return reward, outcome


def backpropagate(node, reward, gamma):
    dist = 0
    while node is not None:
        node.N += 1
        node.W += reward * (gamma ** dist)
        node = node.parent
        dist += 1


# ---- Execution Guard ----
if __name__ == "__main__":
    start_time = time()

    # Parameters
    drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib", "none"]
    d_switch = 7
    total_cycles = 12
    min_size = 1e5
    max_size = 2e10
    depth = 30
    num_rollouts = 1000
    decay_factor = 0.1
    c_param = sqrt(2)
    gamma = 0.9

    ploidy_status = {2.0: 1.5 * 1e9, 3.0: 0.3 * 1e9, 4.0: 0.25 * 1e9}
    tumor_burden_times = [np.array(list(ploidy_status.values()), dtype=float)]
    best_drug_list = []

    root = Node(ploidy_status, 0)
    num_workers = mp.cpu_count()

    for decision in range(total_cycles):
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []

            for _ in range(num_rollouts):
                node = root

                # Selection
                while node.is_fully_expanded(drugs) and not node.is_terminal(total_cycles, min_size, max_size):
                    next_node = select_best_child(node, c_param, decision, decay_factor)
                    if next_node is None: break
                    node = next_node

                # Expansion
                if not node.is_terminal(total_cycles, min_size, max_size):
                    child = expand(node, drugs, d_switch)
                    target_node = child if child else node
                    futures.append((target_node,
                                    executor.submit(rollout, target_node.ploidy_status, depth, drugs, d_switch,
                                                    min_size, max_size)))
                else:
                    futures.append(
                        (node, executor.submit(rollout, node.ploidy_status, 0, drugs, d_switch, min_size, max_size)))

            for node_to_backprop, fut in futures:
                reward, rollout_outcome = fut.result()
                backpropagate(node_to_backprop, reward, gamma)


        def get_eff_q(child):
            age = max(0, decision - child.cycle)
            factor = decay_factor ** age
            return (child.W * factor) / (child.N * factor + epsilon)


        if not root.children:
            print("No children expanded. Exiting.")
            break

        best_drug = max(root.children.items(), key=lambda kv: get_eff_q(kv[1]))[0]
        best_drug_list.append(best_drug)

        print(
            f"Cycle {decision + 1}: best drug is {best_drug} with tumor burden {sum(ploidy_status.values()):.2e} cells")

        # Update state for next cycle
        ploidy_status, ploidies = simulate_next_state(ploidy_status, best_drug, d_switch)
        tumor_burden_times.extend(ploidies)

        temp_sum = np.array([np.sum(arr) for arr in tumor_burden_times])

        if temp_sum[-1] < min_size:
            print("Tumor extinction detected.")
            break
        elif temp_sum[-1] > max_size:
            print("Tumor burden exceeded max size.")
            break

        root = root.children[best_drug]
        root.parent = None

    # --- Plotting ---
    drug_colors = {
        "gemcitabine": "black",
        "bay1895344": "red",
        "alisertib": "green",
        "ispinesib": "blue",
        "none": "yellow"
    }

    TB = np.array([np.sum(arr) for arr in tumor_burden_times])
    plt.figure(figsize=(10, 6))
    plt.subplot(1, 2, 1)
    time_vec = np.arange(len(TB)) * 0.1
    start_idx = 0
    for i, drug in enumerate(best_drug_list):
        end_idx = min(start_idx + (10 * d_switch), len(TB) - 1)
        plt.plot(
            time_vec[start_idx:end_idx + 1],
            TB[start_idx:end_idx + 1],
            color=drug_colors.get(drug, "gray"),
            linewidth=2,
            label=drug if i == 0 or drug not in best_drug_list[:i] else None
        )
        start_idx = end_idx
    plt.xlabel("Time (Days)")
    plt.ylabel("Total Tumor Volume (Cell Count)")
    plt.yscale('log')
    plt.legend(title="Drug Applied")

    plt.subplot(1, 2, 2)
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
    plt.xlabel("Time (Days)")
    plt.ylabel("Ploidy Tumor Volume (Cell Count)")
    plt.tight_layout()
    plt.show()

    print(f"Total time taken: {time() - start_time:.2f} seconds")