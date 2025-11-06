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
        return self.cycle >= total_cycles or sum(self.ploidy_status.values()) < min_size or sum(self.ploidy_status.values()) > max_size

    def is_fully_expanded(self):
        return len(self.children) == len(drugs)

def simulate_next_state(ploidy_status, drug):
    ploidy_cell_count = ploidy_status
    ploidies, t_ode,T_mat_ode, _, _ = ploidy_forcast(ploidy_cell_count, drug, d_switch)
    states=ploidies
    # t=t_ode
    y= T_mat_ode.T #do stuff on tmatode
    new_status = {states[k]: float(y[d_switch][k]) for k in range(len(states))}
    return new_status  # negative burden = reward
def select_best_child(node, c):
    best_score, best_child = -float('inf'), None
    for drug, child in node.children.items():
        Q = child.W / (child.N + 1e-6)
        U = c * np.sqrt(np.log(node.N + 1) / (child.N + 1e-6))
        score = Q + U
        if score > best_score:
            best_score, best_child = score, child
    return best_child


def expand(node):
    untried = [d for d in drugs if d not in node.children]
    drug = random.choice(untried)
    new_ploidy = simulate_next_state(node.ploidy_status, drug)
    child = Node(new_ploidy, node.cycle + 1, parent=node)
    node.children[drug] = child
    return child


def rollout(node, rollout_depth):
    ploidy = dict(node.ploidy_status)
    for _ in range(rollout_depth):
        if sum(ploidy.values()) < min_size or sum(ploidy.values()) > max_size:
            break
        drug = random.choice(drugs)
        ploidy = simulate_next_state(ploidy, drug)
    final_burden = sum(ploidy.values())
    return -final_burden
def backpropagate(node, reward):
    while node is not None:
        node.N += 1
        node.W += reward
        node = node.parent

# Recursively decay W and N of a node and all its children
def decay_node(node, factor):
    node.W *= factor
    node.N *= factor
    for child in node.children.values():
        decay_node(child, factor)


# ---- Main loop ----
drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib"]
d_switch = 7
total_cycles = 25
min_size = 1e8
max_size = 2e10
depth = 30
num_rollouts = 1000 #should be 100000
decay_factor = 0.2  # After a few cycles, old stats are negligible
c=1.4  # Exploration parame ter; higher values encourage exploration of less-visited nodes

#number of cycles counter
cycle_counter=0

# ploidy_status = {2: 0.7*1e9, 4: 0.3*1e9}

ploidy_status = {2.1: 0.9*1e9, 2.3: 0.05*1e9, 3.8: 0.05*1e9}

best_drug_list = []

root = Node(ploidy_status, 0)
for decision in range(total_cycles):
    for _ in range(num_rollouts):  # number of MCTS simulations per decision
        node = root
        # 1. Selection
        while node.is_fully_expanded() and not node.is_terminal():
            node = select_best_child(node, c)
        # 2. Expansion
        if not node.is_terminal():
            child = expand(node)  # expand returns child; no immediate reward
            # 3. Simulation: rollout from child for the remaining cycles
            reward = rollout(child, depth)
            # 4. Backpropagation
            backpropagate(child, reward)
        else:
            reward = -sum(node.ploidy_status.values())
            backpropagate(node, reward)

    # Pick best drug for this decision
    best_drug = max(root.children.items(), key=lambda kv: kv[1].W / (kv[1].N + 1e-6))[0]
    print(f"Cycle {decision + 1}: best drug is {best_drug}")
    best_drug_list.append(best_drug)

    cycle_counter=cycle_counter + 1

    # Move environment forward
    ploidy_status = simulate_next_state(ploidy_status, best_drug)
    root = root.children[best_drug]
    root.parent = None
    decay_node(root, decay_factor)

    if sum(ploidy_status.values()) < min_size:
        print("Tumor extinction.")
        break
    elif sum(ploidy_status.values()) > max_size:
        print("Tumor burden exceeded 2000 mm3, stopping treatment.")
        break
print(f"Final tumor burden: {sum(ploidy_status.values())}")
end = time()
print(f"Total time taken: {end - start:.2f} seconds")
print(cycle_counter)
print(best_drug_list)
print(ploidy_status)

# print(best_drug_list[3])

plt.figure(figsize=(8, 6)) # Create a single figure
plt.xlabel("Time (Days)")
plt.ylabel("Total Tumor Volume (Cell Count)")
# plt.title("Multiple Lines on One Plot")
plt.yscale('log')
for b in range(cycle_counter):
    #simulate fp for given ploidy and drug
    drug=best_drug_list[b]
    ploidies, t_ode,T_mat_ode, _, _ = ploidy_forcast(ploidy_status, drug, d_switch)
    # ploidy_status=ploidies
    total_vol = np.squeeze(np.sum(T_mat_ode,axis=0))
    tvec=t_ode+ b*d_switch
    # plt.yscale('log')
    plt.plot(tvec, total_vol, label=f"Line {b+1}")
    plt.text(tvec[1], total_vol[0], drug, rotation=55, fontsize=12, color='black', ha='left', va='bottom')
    # a=T_mat_ode[0,-1]/sum(T_mat_ode[:,-1])
    # c=T_mat_ode[1,-1]/sum(T_mat_ode[:,-1])
    ploidy_status={ploidies[0]:T_mat_ode[0,-1], ploidies[1] :T_mat_ode[1,-1]}

home_directory = os.path.expanduser("~")

# Define the filename for your figure
filename = "my_figure.png"

# Construct the full path to save the figure
full_path = os.path.join(home_directory, filename)

# Save the figure
plt.savefig(full_path)

print(f"Figure saved to: {full_path}")

ploidy_status = {2.1: 0.9*1e9, 2.3: 0.05*1e9, 3.8: 0.05*1e9}

plt.figure(figsize=(8, 6)) # Create a single figure
plt.xlabel("Time (Days)")
plt.ylabel("Tumor Volume by Ploidy (Cell Count)")
# plt.title("Multiple Lines on One Plot")
plt.yscale('log')
colors = ['black','red', 'green', 'blue', 'yellow', 'purple']
#plot individual ploidy curves
for b in range(cycle_counter):
    drug=best_drug_list[b]
    ploidies, t_ode,T_mat_ode, _, _ = ploidy_forcast(ploidy_status, drug, d_switch)
    print(len(ploidies))
    # ploidy_status=ploidies
    # total_vol = np.squeeze(np.sum(T_mat_ode,axis=0))
    tvec=t_ode+ b*d_switch
    # plt.yscale('log')
    for j in range(len(ploidies)):
        clr=colors[j]
        if b==0:
            val=ploidies[j]
            plt.plot(tvec, T_mat_ode[j,:],color=clr,label=f'Ploidy {val}')
        else:
            plt.plot(tvec, T_mat_ode[j,:],color=clr,label='_nolegend_')
    # plt.text(tvec[1], total_vol[0], drug, rotation=55, fontsize=12, color='black', ha='left', va='bottom')
    # a=T_mat_ode[0,-1]/sum(T_mat_ode[:,-1])
    # c=T_mat_ode[1,-1]/sum(T_mat_ode[:,-1])
    ploidy_status={ploidies[0]:T_mat_ode[0,-1], ploidies[1] :T_mat_ode[1,-1]}

# legendstrings=[]
# for j in range (len(ploidies)):
#     legendstr=f"Ploidy {ploidies(j)}"
#     legendstrings=legendstrings.append(legendstr)
plt.legend()
filename = "my_figure2.png"

# Construct the full path to save the figure
full_path = os.path.join(home_directory, filename)

# Save the figure
plt.savefig(full_path)

print(f"Figure saved to: {full_path}")
# print(ploidies)
print(T_mat_ode)
# # print(a)
# # print(c)
print(total_vol)
# print(t_ode)
# print(tvec)
# print(len(tvec))
# print(len(total_vol))
