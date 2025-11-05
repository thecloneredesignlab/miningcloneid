# 3 cycle prediction, max 30 cycle
import numpy as np
from matplotlib import pyplot as plt

from simulated_data import simulate_ploidy_dynamics_counts

# Have drugs
drugs = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10']

total_cycles = 30
d_target = 3 # Number of switches we look ahead to
d_switch = 7 # NUmber of days between times we are willing to switch drugs

# Define cycles
def permutations_with_replacement(items, length):
    if length == 1:
        return [[item] for item in items]
    else:
        result = []
        for item in items:
            for p in permutations_with_replacement(items, length - 1):
                result.append([item] + p)
        return result

# Generate all permutations
options = permutations_with_replacement(drugs, d_target)
ploidy_status = {2: 7000, 4: 3000}

for i in range(total_cycles):
    lowest_TB = float('inf')
    best_epicycle = None
    for epicycle in options:
        starting_size = sum(ploidy_status.values())
        starting_ploidy_status = ploidy_status.copy()
        epicycle_tumor_burdens = []
        not_working_prev_drug = None
        for j in range(len(epicycle)):
            drug = epicycle[j]
            if not_working_prev_drug == drug:
                continue
            _, y, states = simulate_ploidy_dynamics_counts(ploidy_status, ms_rate=0.01, drug_type=drug, drug_concentration=1.2, t_max=(d_switch+1))
            ploidy_status = {states[k]: float(y[d_switch][k]) for k in range(len(states))}
            tumor_burden_timepoints = y.sum(axis=1, keepdims=True)
            epicycle_tumor_burdens.extend(tumor_burden_timepoints.ravel().tolist()[1:])
            current_size = sum(ploidy_status.values())
            if (current_size - starting_size)/starting_size > 0.2: not_working_prev_drug = drug
            else: starting_size = current_size
        # End of treatment; calculate TB and cost
        TB = sum(ploidy_status.values())
        # print(f"Cycle: {i+1}, Epicycle: {epicycle}, Tumor Burden: {TB}")
        if TB < lowest_TB:
            lowest_TB = TB
            best_epicycle = epicycle
        # Reset ploidy status for next cycle
        ploidy_status = starting_ploidy_status.copy()
    # Apply first treatment in best epicycle
    best_drug = best_epicycle[0]
    print(f"Best epicycle for cycle {i+1}: {best_epicycle}, using Drug: {best_drug}")
    _, y, states = simulate_ploidy_dynamics_counts(ploidy_status, ms_rate=0.01, drug_type=best_drug, drug_concentration=1.2, t_max=8)
    ploidy_status = {states[k]: float(y[7][k]) for k in range(len(states))}
    if sum(ploidy_status.values()) < 1000: break
    elif sum(ploidy_status.values()) > 20000:
        print("Tumor burden exceeded 20000, stopping treatment.")
        break


print(f"Ploidy status at the end of treatment: {ploidy_status}")
print(f"Tumor burden at the end of treatment: {sum(ploidy_status.values())}")
