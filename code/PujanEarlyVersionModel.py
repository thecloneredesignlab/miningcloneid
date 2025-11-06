# Imports
import numpy as np
from scipy.integrate import solve_ivp
import math


def phi_Hill(C, EC50, n, Emax=1):
    """
    Hill-type kill rate function:
    phi(C) = Emax * C^n / (EC50^n + C^n)
    """
    return  (C**n) / (EC50**n + C**n)

def f(ploidy, drug):
    drug = drug.lower()
    def clamp_ec50(x):  return max(x, 1e-12)
    def clamp_n(x):     return max(x, 0.1)   # keep Hill slope sensible

    if drug == "bay1895344":
        n_out = clamp_n(3.85 * math.exp(-0.861 * ploidy) + 0.81)
        ec50  = clamp_ec50(1.04 * math.exp(0.35 * ploidy) - 2.05)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "alisertib":
        n_out = clamp_n(1.0)
        ec50  = clamp_ec50(51.02 * math.exp(-0.62 * ploidy) - 4.78)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "ispinesib":
        n_out = clamp_n(0.94 * math.exp(-0.303 * ploidy) - 0.73)
        ec50  = clamp_ec50(1.185 * math.exp(-0.21 * ploidy) - 0.56)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "gemcitabine":
        n_out = clamp_n(28.92 * math.exp(-0.94 * ploidy) + 0.92)
        ec50  = clamp_ec50(0.004 * math.exp(0.78 * ploidy) - 0.01)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    else:
        raise ValueError(f"Unknown drug: {drug}")
        
# Drug Dependent Dosing Functions

def pulsed_dose(C_peak=5.0, half_life=2.0, period=7.0):
    """
    IV bolus at fixed interval (days); exponential decay between doses.
    C(t) is in the same arbitrary units as C_peak; t in days.
    """
    lam = np.log(2) / max(half_life, 1e-12)
    def C_t(t):
        t = np.asarray(t, dtype=float)
        modt = np.mod(t, period)
        return C_peak * np.exp(-lam * modt)
    return C_t

def oral_pulsed_ss_days(dose=100.0, F=0.7, Vd=70.0, ka_day=1.2, ke_day=0.3, period=1.0, tlag=0.0):
    """
    Steady-state oral dosing every `period` days.
    ka_day, ke_day are per-day rate constants.
    """
    if abs(ka_day - ke_day) < 1e-12:
        # safe limiting form if ka ≈ ke
        def C_t(t):
            t = np.asarray(t, dtype=float)
            tau = float(period)
            tstar = np.mod(t - tlag, tau)
            num = np.exp(-ke_day*tstar)
            den = max(1.0 - np.exp(-ke_day*tau), 1e-12)
            return (F*dose/Vd) * (ke_day * tstar) * num / den
        return C_t

    A = (F * dose * ka_day) / (Vd * (ka_day - ke_day))
    def C_t(t):
        t = np.asarray(t, dtype=float)
        tau = float(period)
        tstar = np.mod(t - tlag, tau)
        term_elim = np.exp(-ke_day*tstar) / max(1.0 - np.exp(-ke_day*tau), 1e-12)
        term_abs  = np.exp(-ka_day*tstar) / max(1.0 - np.exp(-ka_day*tau), 1e-12)
        return A * (term_elim - term_abs)
    return C_t


def get_concentration_curve(drug_name: str, **overrides):
    """
    Returns C(t) for the given drug, using route from `drug_dosing_schedules`.
    - For IV: uses pulsed_dose (days).
    - For Oral: uses steady-state Bateman (days).
    You can override any PK parameter via kwargs.
    """
    drug_name = drug_name.lower()
    route = drug_dosing_schedules.get(drug_name)
    if route is None:
        raise KeyError(f"Unknown drug: {drug_name}")

    # Gather per-drug params (if any), then apply overrides
    per_drug = PER_DRUG.get(drug_name, {})
    if route == "IV":
        params = {**IV_DEFAULTS, **per_drug, **overrides}
        return pulsed_dose(**params)
    elif route == "Oral":
        base = {**ORAL_DEFAULTS, **per_drug, **overrides}
        return oral_pulsed_ss_days(**base)
    else:
        raise ValueError(f"Unsupported route '{route}' for {drug_name}")
    
def dose_times(drug_name, t_span, tlag=0.0):
    drug_name = drug_name.lower()
    route = drug_dosing_schedules[drug_name]
    if route == "IV":
        period = (PER_DRUG.get(drug_name, IV_DEFAULTS)).get("period", IV_DEFAULTS["period"])
    else:  # Oral daily by default unless overridden
        period = (PER_DRUG.get(drug_name, ORAL_DEFAULTS)).get("period", ORAL_DEFAULTS["period"])
    t0, t1 = t_span
    # all multiples of period shifted by tlag within [t0, t1]
    start_k = math.ceil((t0 - tlag) / period)
    end_k   = math.floor((t1 - tlag) / period)
    return tlag + period * np.arange(start_k, end_k + 1)

def get_dose_period(drug_name: str) -> float:
    """
    Return the dosing period (days) for the drug based on route and PER_DRUG fallback.
    """
    dn = drug_name.lower()
    route = drug_dosing_schedules.get(dn)
    if route == "IV":
        return PER_DRUG.get(dn, IV_DEFAULTS).get("period", IV_DEFAULTS["period"])
    elif route == "Oral":
        return PER_DRUG.get(dn, ORAL_DEFAULTS).get("period", ORAL_DEFAULTS["period"])
    else:
        raise KeyError(f"Unknown route for {drug_name}")

def build_times_with_doses(t_span, dt, drug_name, tlag: float = 0.0,
                           include_days: bool = True, eps: float = 1e-8):
    """
    Build a monotonically increasing time grid that:
      - samples with base step dt,
      - includes all dose instants and a tiny +eps after each,
      - optionally includes integer day ticks.
    """
    t0, t1 = float(t_span[0]), float(t_span[1])
    tau = float(get_dose_period(drug_name))

    # Base uniform grid
    base = np.arange(t0, t1 + 0.5*dt, dt, dtype=float)

    # Dose instants within [t0, t1]
    # Find first k such that tlag + k*tau >= t0
    k_start = math.ceil((t0 - tlag) / tau)
    k_end   = math.floor((t1 - tlag) / tau)
    if k_end >= k_start:
        doses = (tlag + tau * np.arange(k_start, k_end + 1, dtype=float))
    else:
        doses = np.array([], dtype=float)

    # Optional whole-day ticks
    days = np.arange(math.ceil(t0), math.floor(t1) + 1, dtype=float) if include_days else np.array([], dtype=float)

    # Combine: base, doses, and post-dose samples
    all_times = np.concatenate([base, doses, np.minimum(doses + eps, t1), days])
    times = np.unique(all_times)
    # Ensure endpoints are present
    if times[0] > t0:   times = np.insert(times, 0, t0)
    if times[-1] < t1:  times = np.append(times, t1)
    return times

# System building and simulation functions

def build_multiclonal(initial_tumor_fraction, drug, r, K, C_fn=None, r_by_ploidy=None):
    ploidies = np.array(sorted(initial_tumor_fraction.keys()), dtype=float)
    T0 = np.array([initial_tumor_fraction[p] for p in ploidies], dtype=float)
    phi_params = [f(p, drug) for p in ploidies]
    r_vec = np.array([r_by_ploidy.get(p, r) if r_by_ploidy else r for p in ploidies], dtype=float)
    C_func = C_fn if C_fn is not None else (lambda t: 0.0)
    return ploidies, T0, phi_params, r_vec, C_func

def rhs_ode(t, T, r_vec, K, phi_params, C_func):
    Tsum = np.sum(T)
    C = C_func(t)
    phi_vals = np.array([phi_Hill(C, **pp) for pp in phi_params])
    return (r_vec * (1 - Tsum / K) - phi_vals) * T

def simulate_ode_piecewise(initial_tumor_fraction, drug, t_span, r, K, C_fn=None,
                           r_by_ploidy=None, max_step=None, atol=1e-9, rtol=1e-7,
                           tlag=0.0, include_points=None, t_eval=None, eps=1e-10):
    """
    Piecewise integrate across dose times, *and* honor an explicit t_eval grid.

    - include_points: extra times to force into output (e.g., whole days)
    - t_eval: global grid you want in the final output (e.g., the TIMES you also gave SDE)
    """
    ploidies, T0, phi_params, r_vec, C_func = build_multiclonal(
        initial_tumor_fraction, drug, r, K, C_fn, r_by_ploidy
    )
    t0, t1 = float(t_span[0]), float(t_span[1])
    if t0 >= t1:
        raise ValueError("t_span must have t0 < t1")

    # --- Build hard breaks at dosing times ---
    dtimes = dose_times(drug, (t0, t1), tlag=tlag)
    dtimes = dtimes[(dtimes > t0) & (dtimes < t1)]  # interior only

    breaks = np.array([t0, *dtimes.tolist(), t1], dtype=float)

    # --- Build the set of global evaluation points we must return ---
    add = []
    if include_points is not None:
        add.extend(include_points)
    if t_eval is not None:
        add.extend(t_eval)

    # tiny sample just after each dose time (right-continuous for C(t))
    post_dose = np.minimum(dtimes + eps, t1) if dtimes.size else np.array([], dtype=float)

    # global required samples
    required_global = np.unique(
        np.concatenate(([t0, t1], np.asarray(add, dtype=float) if len(add) else [], post_dose))
    )

    # --- Integrate segment by segment ---
    t_all = []
    Y_all = []
    T_init = T0.copy()

    for a, b in zip(breaks[:-1], breaks[1:]):
        if b <= a + 1e-14:
            continue  # degenerate segment (can happen with floating rounding)

        # All requested eval points that fall in [a, b]
        mask = (required_global >= a) & (required_global <= b)
        seg_points = required_global[mask]

        # Ensure we hit a, a+eps (if room), and b within the segment
        seg_core = np.array([a, min(a + eps, b), b], dtype=float)
        seg_eval = np.unique(np.concatenate((seg_points, seg_core)))

        # If seg_eval collapses to a single point (rare), add b to keep solver happy
        if seg_eval.size == 1:
            seg_eval = np.unique(np.array([seg_eval[0], b], dtype=float))

        sol = solve_ivp(
            rhs_ode, (a, b), T_init,
            t_eval=seg_eval, args=(r_vec, K, phi_params, C_func),
            rtol=rtol, atol=atol, max_step=max_step
        )

        # Append, avoiding duplicate of the left endpoint
        if t_all and abs(sol.t[0] - t_all[-1]) < 1e-12:
            t_all.extend(sol.t[1:].tolist())
            Y_all.append(sol.y[:, 1:])
        else:
            t_all.extend(sol.t.tolist())
            Y_all.append(sol.y)

        # Restart from the last state (state is continuous; C(t) jumped at a)
        T_init = sol.y[:, -1]

    # Concatenate segments
    t_all = np.array(t_all, dtype=float)
    T_mat = np.hstack(Y_all)
    T_total = T_mat.sum(axis=0)

    # Optionally, if the caller provided a t_eval, return exactly that ordering
    # by interpolating the piecewise solution to the requested grid.
    # Usually unnecessary because we already forced those times into seg_eval,
    # but this preserves *ordering* identical to t_eval if desired.
    if t_eval is not None:
        t_req = np.asarray(t_eval, dtype=float)
        # Ensure coverage
        if t_req[0] < t_all[0] - 1e-12 or t_req[-1] > t_all[-1] + 1e-12:
            raise ValueError("t_eval outside integration range.")
        # Interpolate each clone onto t_req
        T_mat_req = np.vstack([np.interp(t_req, t_all, T_mat[i, :]) for i in range(T_mat.shape[0])])
        return t_req, ploidies, T_mat_req, T_mat_req.sum(axis=0)

    return t_all, ploidies, T_mat, T_total

def simulate_sde(initial_tumor_fraction, drug, t_span, dt, r, K, n_sims,
                 C_fn=None, r_by_ploidy=None, beta_by_ploidy=None, rng=None,
                 times=None, rho=0.0):
    """
    If `times` is provided, it is used directly (dose-aware); otherwise build with dt.
    `rho`: fraction (0..1) of systemic (shared) noise per clone.
    """
    ploidies, T0, phi_params, r_vec, C_func = build_multiclonal(initial_tumor_fraction, drug, r, K, C_fn, r_by_ploidy)
    M = len(ploidies)

    if beta_by_ploidy is None:
        beta_by_ploidy = {float(p): 0.02 for p in ploidies}
    beta_vec = np.array([beta_by_ploidy[float(p)] for p in ploidies], dtype=float)

    if rng is None:
        rng = np.random.default_rng(12345)

    # --- shared dose-aware grid ---
    if times is None:
        times = build_times_with_doses(t_span, dt, drug_name=drug, tlag=0.0, include_days=True, eps=1e-8)
    times = np.asarray(times, dtype=float)
    N = times.size
    Tpaths = np.zeros((n_sims, M, N))
    Tpaths[:, :, 0] = T0

    # Exponential Euler update (positivity-preserving)
    for k in range(N - 1):
        t  = times[k]
        dtk = times[k+1] - t
        sqrt_dtk = math.sqrt(max(dtk, 0.0))

        C = C_func(t)
        phi_vals = np.array([phi_Hill(C, **pp) for pp in phi_params], dtype=float)  # (M,)

        T_curr = Tpaths[:, :, k]                             # (n_sims, M)
        Tsum   = T_curr.sum(axis=1, keepdims=True)           # (n_sims, 1)
        mu     = (r_vec * (1.0 - Tsum / K)) - phi_vals       # broadcast -> (n_sims, M)

        if rho > 0.0:
            beta_sys  = rho * beta_vec
            beta_idio = (1.0 - rho) * beta_vec
            Z_sys  = rng.normal(0.0, 1.0, size=(T_curr.shape[0], 1))
            Z_idio = rng.normal(0.0, 1.0, size=T_curr.shape)
            dW = beta_sys * sqrt_dtk * Z_sys + beta_idio * sqrt_dtk * Z_idio
            beta_eff_sq = (beta_sys**2) + (beta_idio**2)
        else:
            Z = rng.normal(0.0, 1.0, size=T_curr.shape)
            dW = beta_vec * sqrt_dtk * Z
            beta_eff_sq = beta_vec**2

        expo = (mu - 0.5 * beta_eff_sq) * dtk + dW
        T_next = T_curr * np.exp(expo)
        Tpaths[:, :, k+1] = T_next

    return times, ploidies, Tpaths

# Main forecasting function for ploidy-specific cell counts

def ploidy_forcast(ploidy_cell_count, drug, T, R_BASE=.4, K_CAP = 4e10, beta_by_ploidy=None, DT = 0.1, N_SIMS =200):
    ploidies = np.array(sorted(ploidy_cell_count.keys()), dtype=float)
    if beta_by_ploidy is None:
        beta_by_ploidy = {float(p): 0.02 for p in ploidies}
    C_fn = get_concentration_curve(drug)
    T_SPAN = (0.0, T)
    DRUG = drug.lower()
    initial_tumor_fraction = ploidy_cell_count
    # Build the shared grid once (the same one you pass to SDE)
    TIMES = build_times_with_doses(T_SPAN, DT, drug_name=DRUG, tlag=0.0, include_days=True, eps=1e-8)

    # ODE (piecewise) on that exact grid
    t_ode, ploidies, T_mat_ode, T_total_ode = simulate_ode_piecewise(
        initial_tumor_fraction, DRUG, t_span=T_SPAN, r=R_BASE, K=K_CAP, C_fn=C_fn,
        max_step=0.05, tlag=0.0, t_eval=TIMES
    )

    # SDE on the same grid (as you already set up)
    t_sde, ploidies_sde, Tpaths = simulate_sde(
        initial_tumor_fraction, DRUG, t_span=T_SPAN, dt=DT, r=R_BASE, K=K_CAP, n_sims=N_SIMS,
        C_fn=C_fn, beta_by_ploidy=beta_by_ploidy, times=TIMES, rho=0.0
    )
    assert np.allclose(ploidies, ploidies_sde)
        
    # Aggregates for SDE
    T_total_sde = Tpaths.sum(axis=1)  # (n_sims, N)
    # Extract final time and values
    return ploidies, t_ode,T_mat_ode,t_sde, Tpaths.mean(axis=0)

# Drug dosing schedules and parameters
drug_dosing_schedules = {
    "volasertib"    : "IV",
    "umi-77"        : "IV",
    "tegafur"       : "Oral",
    "tas"           : "Oral",
    "pf473776_chk"  : "IV",
    "osi-027"       : "Oral",
    "alisertib"     : "IV",
    "5-azacytidine" : "Oral",
    "abt-199"       : "Oral",
    "abt-263"       : "Oral",
    "capecitabine"  : "Oral",
    "ceralasertib"  : "Oral",
    "cytarabine"    : "IV",
    "gemcitabine"   : "IV",
    "bay1895344"    : "IV",
    "ispinesib"     : "IV",
    "navitoclax"    : "Oral",
    "adavosertib"   : "Oral",
}

IV_DEFAULTS = dict(C_peak=1.0, half_life=0.5, period=7.0)  # weekly IV as a baseline

ORAL_DEFAULTS = {
    "dose": 100.0, "F": 0.6, "Vd": 60.0,
    "ka_day": 2.0,                   # ~ ka = 2 / day (t1/2(abs) ≈ 0.35 d ≈ 8.3 h)
    "ke_day": 0.7,                   # ~ ke = 0.7 / day (t1/2 ≈ 1 d)
    "period": 1.0,                   # daily
    "tlag": 0.0
}
PER_DRUG = {
    # IV examples
    "volasertib":   {"C_peak": 1.0, "half_life": 4.0, "period": 7.0},
    "alisertib":    {"C_peak": 42.5, "half_life": 0.875, "period": 0.5},
    "cytarabine":   {"C_peak": 1.0, "half_life": 0.2, "period": 3.5},
    "gemcitabine":  {"C_peak": 239, "half_life": 0.05, "period": 7.0},
    "ispinesib":    {"C_peak": 2.1, "half_life": 1.04, "period": 7},
    "umi-77":       {"C_peak": 1.0, "half_life": 0.8, "period": 7.0},
    "navitoclax":   {"C_peak": 1.0, "half_life": 0.73, "period": 1},
    "bay1895344":   {"C_peak": 6.2, "half_life": 0.50, "period": 0.5},

    # Oral examples
    "osi-027":      {"dose": 50, "F": 0.5, "Vd": 80, "ka_day": 1.8, "ke_day": 0.5, "period": 1.0},
    "abt-199":      {"dose": 100, "F": 0.6, "Vd": 250, "ka_day": 1.0, "ke_day": 0.3, "period": 1.0},
    "abt-263":      {"dose": 100, "F": 0.6, "Vd": 120, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    # "navitoclax":   {"dose": 100, "F": 0.6, "Vd": 120, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "ceralasertib": {"dose": 80,  "F": 0.5, "Vd": 100,"ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    # "bay1895344":   {"dose": 60,  "F": 0.5, "Vd": 90, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "adavosertib":  {"dose": 100, "F": 0.6, "Vd": 65, "ka_day": 2.4, "ke_day": 0.6, "period": 1.0},
    "tas":          {"dose": 60,  "F": 0.5, "Vd": 40, "ka_day": 2.4, "ke_day": 0.7, "period": 1.0},
    "tegafur":      {"dose": 40, "F": 0.5,"Vd": 45, "ka_day": 1.6, "ke_day": 0.5, "period": 1.0},
    "capecitabine": {"dose": 100, "F": 0.8, "Vd": 40, "ka_day": 3.0, "ke_day": 0.6, "period": 1.0},
    "5-azacytidine":{"dose": 100, "F": 0.2, "Vd": 40, "ka_day": 3.0, "ke_day": 2.0, "period": 1.0},  # oral-like placeholder
}


if __name__ == "__main__":
    
    #sample usage 
    ploidy_cell_count = {2.0: 1e3, 4.0: 2e2, 6.0: 1.5e3} #dict{ploidy: cell count}
    drug = "gemcitabine"  
    T = 30
    ploidies, t_ode,T_mat_ode, t_sde, Tpaths = ploidy_forcast(ploidy_cell_count, drug, T)
    
