import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def simulate_ploidy_dynamics_counts(ploidy_init, ms_rate, drug_type, drug_concentration,
                                   t_max=0, dt=1, base_growth_rate=0.05):
    """
    Simulates ploidy-specific population dynamics under drug selection.
    Ploidy_init: dict of absolute cell counts (e.g. {2: 7000, 4: 3000})
    """
    ploidy_states = sorted(ploidy_init.keys())
    y0 = np.array([ploidy_init[p] for p in ploidy_states], dtype=float)

    # Hard-coded EC50 vs ploidy relationship (linear)
    base_ec50 = {"D1": 1.0, "D2": 0.8, "D3": 1.2, "D4": 1.2, "D5": 1.2, "D6": 1.2, "D7": 1.2, "D8": 1.2, "D9": 1.0, "D10": 0.001}.get(drug_type, 1.0)
    ec50s = {p: base_ec50 * (1 + 0.2 * (p - 2)) for p in ploidy_states}

    def dydt(y, t):
        dydt = np.zeros_like(y)
        for i, p in enumerate(ploidy_states):
            # Growth rate (drug-free)
            r = base_growth_rate
            # Drug-induced death (Hill equation)
            drug_kill = drug_concentration / (drug_concentration + ec50s[p])
            dydt[i] = y[i] * (r * (1 - drug_kill))
        return dydt

    t = np.arange(0, t_max, dt)
    y = odeint(dydt, y0, t)
    return t, y, ploidy_states


# Example usage with absolute counts
ploidy_init = {2: 7000, 4: 3000}
t, y, states = simulate_ploidy_dynamics_counts(ploidy_init, ms_rate=0.01, drug_type="D1", drug_concentration=1.2)

# print(y)

# Add y lengthwise to calculate TB
# TB = y.sum(axis=1, keepdims=True)
# print(TB)


# Plot absolute counts
# plt.figure(figsize=(8, 5))
# for i, p in enumerate(states):
#     plt.plot(t, TB, label=f"{p}n (absolute)")
# plt.xlabel("Time (a.u.)")
# plt.ylabel("Cell count")
# plt.title("Absolute population dynamics under drug-induced selection")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()

# # Plot relative composition
# plt.figure(figsize=(8, 5))
# fractions = y / y.sum(axis=1, keepdims=True)
# for i, p in enumerate(states):
#     plt.plot(t, fractions[:, i], label=f"{p}n (fraction)")
# plt.xlabel("Time (a.u.)")
# plt.ylabel("Fraction of total population")
# plt.title("Relative ploidy composition over time")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()
