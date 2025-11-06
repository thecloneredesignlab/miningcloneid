import numpy as np
import random
from copy import deepcopy
from PujanEarlyVersionModel import ploidy_forcast
# from FakeODEModelG3 import simulate_ploidy_dynamics_counts
from time import time

start = time()

drugs = ["gemcitabine", "bay1895344", "alisertib", "ispinesib"]
#drugs=['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10']
d_switch = 7
total_cycles = 11

# ploidy_status={2: 7000, 4: 3000}
# t, y, states=simulate_ploidy_dynamics_counts(ploidy_status, ms_rate=0.01, drug_type='D2', drug_concentration=1.2, t_max=(d_switch + 1))
# print(states)
# print(t)
# print(y)

# def simulate_next_state(ploidy_status, drug):

# ploidies, t_ode,T_mat_ode, _, _ = ploidy_forcast({2.0: 1e3, 4.0: 2e2, 6.0: 1.5e3}, "alisertib", d_switch)
# print("ploidies") #list of ploidy values
# print(ploidies)
# print("t_ode")
# print(t_ode) #vector of time steps between 0 and 7 days
# print("T_mat_ode")
# print(len(T_mat_ode))
# print(T_mat_ode) #array of arrays with cell count trajectory for the different ploidies




class Node:
    def __init__(self, ploidy_status, cycle, parent=None):
        self.ploidy_status = ploidy_status
        self.cycle = cycle
        self.parent = parent
        self.children = {}
        self.N = 0
        self.W = 0.0  # cumulative reward (negative tumor burden)

    def is_terminal(self):
        return self.cycle >= total_cycles or sum(self.ploidy_status.values()) < 1000

    def is_fully_expanded(self):
        return len(self.children) == len(drugs)

#ploidy_forcast(ploidy_cell_count, drug, T)
#ploidy_status = {2: 7000, 4: 3000}
#ploidy_cell_count = {2.0: 1e3, 4.0: 2e2, 6.0: 1.5e3}
#_, y, states=simulate_ploidy_dynamics_counts(ploidy_status, ms_rate=0.01, drug_type=drug, drug_concentration=1.2, t_max=(d_switch + 1))
def simulate_next_state(ploidy_status, drug):
    ploidy_cell_count = ploidy_status
    ploidies, t_ode,T_mat_ode, _, _ = ploidy_forcast(ploidy_cell_count, drug, d_switch)
    states=ploidies
    #new_status = {states[k]: float(y[d_switch][k]) for k in range(len(states))}
    #new_status= # do stuff w t_mat_ode??? over the ploidies
    # TB = sum(new_status.values())
    t=t_ode
    y= T_mat_ode.T #do stuff on tmatode
    new_status = {states[k]: float(y[d_switch][k]) for k in range(len(states))}
    TB = sum(new_status.values())

    return new_status, -TB  # negative burden = reward


def select_best_child(node, c=1.4):
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
    new_ploidy, reward = simulate_next_state(node.ploidy_status, drug)
    child = Node(new_ploidy, node.cycle + 1, parent=node)
    node.children[drug] = child
    return child, reward


def rollout(node, rollout_depth=5):
    ploidy = deepcopy(node.ploidy_status)
    reward = 0
    for _ in range(rollout_depth):
        if sum(ploidy.values()) < 1000:
            break
        drug = random.choice(drugs)
        ploidy, reward = simulate_next_state(ploidy, drug)
    return reward


def backpropagate(node, reward):
    while node is not None:
        node.N += 1
        node.W += reward
        node = node.parent


# ---- Main loop ----
ploidy_status = {2: 7000, 4: 3000}
root = Node(ploidy_status, 0)

for decision in range(total_cycles):
    for _ in range(1000):  # number of MCTS simulations per decision
        node = root
        # 1. Selection
        while node.is_fully_expanded() and not node.is_terminal():
            node = select_best_child(node)
        # 2. Expansion
        if not node.is_terminal():
            node, reward = expand(node)
            # 3. Simulation
            reward += rollout(node)
        else:
            reward = -sum(node.ploidy_status.values())
        # 4. Backpropagation
        backpropagate(node, reward)

    # Pick best drug for this decision
    best_drug = max(root.children.items(), key=lambda kv: kv[1].W / (kv[1].N + 1e-6))[0]
    print(f"Cycle {decision + 1}: best drug is {best_drug}")

    # Move environment forward
    ploidy_status, _ = simulate_next_state(ploidy_status, best_drug)
    root = Node(ploidy_status, decision + 1)  # reset root for next decision

    if sum(ploidy_status.values()) < 1000:
        print("Tumor eliminated early.")
        break
    elif sum(ploidy_status.values()) > 20000:
        print("Tumor burden exceeded 20000, stopping treatment.")
        break

print(f"Final tumor burden: {sum(ploidy_status.values())}")
end = time()
print(f"Total time taken: {end - start:.2f} seconds")



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OLD VERSION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# import numpy as np
# import random
# from copy import deepcopy
# from FakeODEModelG3 import simulate_ploidy_dynamics_counts
# from time import time

# start = time()

# drugs = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10']
# d_switch = 7
# total_cycles = 11


# class Node:
#     def __init__(self, ploidy_status, cycle, parent=None):
#         self.ploidy_status = ploidy_status
#         self.cycle = cycle
#         self.parent = parent
#         self.children = {}
#         self.N = 0
#         self.W = 0.0  # cumulative reward (negative tumor burden)

#     def is_terminal(self):
#         return self.cycle >= total_cycles or sum(self.ploidy_status.values()) < 1000

#     def is_fully_expanded(self):
#         return len(self.children) == len(drugs)


# def simulate_next_state(ploidy_status, drug):
#     _, y, states = simulate_ploidy_dynamics_counts(ploidy_status, ms_rate=0.01, drug_type=drug, drug_concentration=1.2, t_max=(d_switch + 1))
#     new_status = {states[k]: float(y[d_switch][k]) for k in range(len(states))}
#     TB = sum(new_status.values())
#     return new_status, -TB  # negative burden = reward


# def select_best_child(node, c=1.4):
#     best_score, best_child = -float('inf'), None
#     for drug, child in node.children.items():
#         Q = child.W / (child.N + 1e-6)
#         U = c * np.sqrt(np.log(node.N + 1) / (child.N + 1e-6))
#         score = Q + U
#         if score > best_score:
#             best_score, best_child = score, child
#     return best_child


# def expand(node):
#     untried = [d for d in drugs if d not in node.children]
#     drug = random.choice(untried)
#     new_ploidy, reward = simulate_next_state(node.ploidy_status, drug)
#     child = Node(new_ploidy, node.cycle + 1, parent=node)
#     node.children[drug] = child
#     return child, reward


# def rollout(node, rollout_depth=5):
#     ploidy = deepcopy(node.ploidy_status)
#     reward_sum = 0.0
#     for _ in range(rollout_depth):
#         if sum(ploidy.values()) < 1000:
#             break
#         drug = random.choice(drugs)
#         ploidy, reward = simulate_next_state(ploidy, drug)
#         reward_sum += reward
#     return reward_sum


# def backpropagate(node, reward):
#     while node is not None:
#         node.N += 1
#         node.W += reward
#         node = node.parent


# # ---- Main loop ----
# ploidy_status = {2: 7000, 4: 3000}
# root = Node(ploidy_status, 0)

# for decision in range(total_cycles):
#     for _ in range(1000):  # number of MCTS simulations per decision
#         node = root
#         # 1. Selection
#         while node.is_fully_expanded() and not node.is_terminal():
#             node = select_best_child(node)
#         # 2. Expansion
#         if not node.is_terminal():
#             node, reward = expand(node)
#             # 3. Simulation
#             reward += rollout(node)
#         else:
#             reward = -sum(node.ploidy_status.values())
#         # 4. Backpropagation
#         backpropagate(node, reward)

#     # Pick best drug for this decision
#     best_drug = max(root.children.items(), key=lambda kv: kv[1].W / (kv[1].N + 1e-6))[0]
#     print(f"Cycle {decision + 1}: best drug is {best_drug}")

#     # Move environment forward
#     ploidy_status, _ = simulate_next_state(ploidy_status, best_drug)
#     root = Node(ploidy_status, decision + 1)  # reset root for next decision

#     if sum(ploidy_status.values()) < 1000:
#         print("Tumor eliminated early.")
#         break
#     elif sum(ploidy_status.values()) > 20000:
#         print("Tumor burden exceeded 20000, stopping treatment.")
#         break

# print(f"Final tumor burden: {sum(ploidy_status.values())}")
# end = time()
# print(f"Total time taken: {end - start:.2f} seconds")
