# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math
from scipy.special import comb

plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
plt.rcParams["text.latex.preamble"] += r"\usepackage{amssymb}"


# Drug effect functions
def phi_Hill(C, EC50, n, Emax=1.0):
    C = float(C)
    EC50 = max(float(EC50), 1e-12)
    n = max(float(n), 0.1)
    return Emax * (C ** n) / (EC50 ** n + C ** n)


def make_phi_vec_fn(Ns, drug_name, f_param_fn, baseline=23.0):
    """
    Returns phi_vec_fn(C, t) -> vector of kill rates for each N in Ns.
    """
    params = []
    for N in Ns:
        ploidy = N / baseline
        params.append(f_param_fn(ploidy, drug_name))

    def phi_vec_fn(C, t=None):
        return np.array([phi_Hill(C, **pp) for pp in params], dtype=float)

    return phi_vec_fn


def f(ploidy, drug):
    drug = drug.lower()

    def clamp_ec50(x):
        return max(x, 1e-12)

    def clamp_n(x):
        return max(x, 0.1)  # keep Hill slope sensible

    if drug == "bay1895344":
        n_out = clamp_n(3.85 * math.exp(-0.861 * ploidy) + 0.81)
        ec50 = clamp_ec50(1.04 * math.exp(0.35 * ploidy) - 2.05)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "alisertib":
        n_out = clamp_n(1.0)
        ec50 = clamp_ec50(51.02 * math.exp(-0.62 * ploidy) - 4.78)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "ispinesib":
        n_out = clamp_n(0.94 * math.exp(-0.303 * ploidy) - 0.73)
        ec50 = clamp_ec50(1.185 * math.exp(-0.21 * ploidy) - 0.56)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "gemcitabine":
        n_out = clamp_n(28.92 * math.exp(-0.94 * ploidy) + 0.92)
        ec50 = clamp_ec50(0.004 * math.exp(0.78 * ploidy) - 0.01)
        return dict(EC50=ec50, n=n_out, Emax=1.0)

    elif drug == "none":
        n_out = clamp_n(28.92 * math.exp(-0.94 * ploidy) + 0.92)
        ec50 = clamp_ec50(0.004 * math.exp(0.78 * ploidy) - 0.01)
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
            num = np.exp(-ke_day * tstar)
            den = max(1.0 - np.exp(-ke_day * tau), 1e-12)
            return (F * dose / Vd) * (ke_day * tstar) * num / den

        return C_t

    A = (F * dose * ka_day) / (Vd * (ka_day - ke_day))

    def C_t(t):
        t = np.asarray(t, dtype=float)
        tau = float(period)
        tstar = np.mod(t - tlag, tau)
        term_elim = np.exp(-ke_day * tstar) / max(1.0 - np.exp(-ke_day * tau), 1e-12)
        term_abs = np.exp(-ka_day * tstar) / max(1.0 - np.exp(-ka_day * tau), 1e-12)
        return A * (term_elim - term_abs)

    return C_t


# Drug dosing schedules and parameters
drug_dosing_schedules = {
    "volasertib": "IV",
    "umi-77": "IV",
    "tegafur": "Oral",
    "tas": "Oral",
    "pf473776_chk": "IV",
    "osi-027": "Oral",
    "alisertib": "IV",
    "5-azacytidine": "Oral",
    "abt-199": "Oral",
    "abt-263": "Oral",
    "capecitabine": "Oral",
    "ceralasertib": "Oral",
    "cytarabine": "IV",
    "gemcitabine": "IV",
    "none": "IV",
    "bay1895344": "IV",
    "ispinesib": "IV",
    "navitoclax": "Oral",
    "adavosertib": "Oral",
}

IV_DEFAULTS = dict(C_peak=1.0, half_life=0.5, period=7.0)  # weekly IV as a baseline

ORAL_DEFAULTS = {
    "dose": 100.0, "F": 0.6, "Vd": 60.0,
    "ka_day": 2.0,  # ~ ka = 2 / day (t1/2(abs) ≈ 0.35 d ≈ 8.3 h)
    "ke_day": 0.7,  # ~ ke = 0.7 / day (t1/2 ≈ 1 d)
    "period": 1.0,  # daily
    "tlag": 0.0
}
PER_DRUG = {
    # IV examples
    "volasertib": {"C_peak": 1.0, "half_life": 4.0, "period": 7.0},
    "alisertib": {"C_peak": 42.5, "half_life": 0.875, "period": 0.5},
    "cytarabine": {"C_peak": 1.0, "half_life": 0.2, "period": 3.5},
    "gemcitabine": {"C_peak": 239, "half_life": 0.05, "period": 7.0},
    "ispinesib": {"C_peak": 2.1, "half_life": 1.04, "period": 7},
    "umi-77": {"C_peak": 1.0, "half_life": 0.8, "period": 7.0},
    "navitoclax": {"C_peak": 1.0, "half_life": 0.73, "period": 1},
    "bay1895344": {"C_peak": 6.2, "half_life": 0.50, "period": 0.5},
    "none": {"C_peak": 0, "half_life": 0.50, "period": 7},

    # Oral examples
    "osi-027": {"dose": 50, "F": 0.5, "Vd": 80, "ka_day": 1.8, "ke_day": 0.5, "period": 1.0},
    "abt-199": {"dose": 100, "F": 0.6, "Vd": 250, "ka_day": 1.0, "ke_day": 0.3, "period": 1.0},
    "abt-263": {"dose": 100, "F": 0.6, "Vd": 120, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    # "navitoclax":   {"dose": 100, "F": 0.6, "Vd": 120, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "ceralasertib": {"dose": 80, "F": 0.5, "Vd": 100, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    # "bay1895344":   {"dose": 60,  "F": 0.5, "Vd": 90, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "adavosertib": {"dose": 100, "F": 0.6, "Vd": 65, "ka_day": 2.4, "ke_day": 0.6, "period": 1.0},
    "tas": {"dose": 60, "F": 0.5, "Vd": 40, "ka_day": 2.4, "ke_day": 0.7, "period": 1.0},
    "tegafur": {"dose": 40, "F": 0.5, "Vd": 45, "ka_day": 1.6, "ke_day": 0.5, "period": 1.0},
    "capecitabine": {"dose": 100, "F": 0.8, "Vd": 40, "ka_day": 3.0, "ke_day": 0.6, "period": 1.0},
    "5-azacytidine": {"dose": 100, "F": 0.2, "Vd": 40, "ka_day": 3.0, "ke_day": 2.0, "period": 1.0},
    # oral-like placeholder
}


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
    end_k = math.floor((t1 - tlag) / period)
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
    base = np.arange(t0, t1 + 0.5 * dt, dt, dtype=float)

    # Dose instants within [t0, t1]
    # Find first k such that tlag + k*tau >= t0
    k_start = math.ceil((t0 - tlag) / tau)
    k_end = math.floor((t1 - tlag) / tau)
    if k_end >= k_start:
        doses = (tlag + tau * np.arange(k_start, k_end + 1, dtype=float))
    else:
        doses = np.array([], dtype=float)

    # Optional whole-day ticks
    days = np.arange(math.ceil(t0), math.floor(t1) + 1, dtype=float) if include_days else np.array([], dtype=float)

    # Combine: base, doses, and post-dose samples
    all_times = np.concatenate([base, doses, np.minimum(doses, t1), days])
    times = np.unique(all_times)
    # Ensure endpoints are present
    if times[0] > t0:   times = np.insert(times, 0, t0)
    if times[-1] < t1:  times = np.append(times, t1)
    return times


# Missegregation Matrix Construction
def build_misseg_matrix(N_min: int, N_max: int, beta: float, renormalize=True):
    """
    Build M where M[i,j] = expected # of daughters with N=Ns[i]
    produced by one division of a mother with N=Ns[j].

    If renormalize=True: columns sum to 2 within [N_min,N_max]
    (mass that would go out of bounds is redistributed in-bounds).
    If renormalize=False: leakage out of bounds is allowed (col sum <= 2).
    """
    Ns = np.arange(N_min, N_max + 1, dtype=int)
    m = len(Ns)
    P = np.zeros((m, m), dtype=float)  # daughter-probabilities

    for j, q in enumerate(Ns):
        col = np.zeros(m, dtype=float)

        # K ~ Binomial(q, beta)
        for k in range(q + 1):
            pk = comb(q, k) * (beta ** k) * ((1.0 - beta) ** (q - k))
            if pk == 0.0:
                continue

            # X ~ Binomial(k, 0.5); Delta = 2X - k
            for x in range(k + 1):
                delta = 2 * x - k
                p = q + delta
                if p < N_min or p > N_max:
                    continue
                px = comb(k, x) * (0.5 ** k)
                col[p - N_min] += pk * px

        if renormalize:
            s = col.sum()
            if s > 0:
                col /= s
            else:
                # fallback: keep daughter at same state
                col[q - N_min] = 1.0

        P[:, j] = col

    M = 2.0 * P  # expected daughters
    return Ns, M


# Karyotype ODE RHS
def rhs_karyotype_ode(t, T, r, Kcap, M, phi_vec_fn, C_fn):
    Tsum = float(np.sum(T))
    crowd = 1.0 - Tsum / Kcap
    C = float(C_fn(t))
    phi = np.asarray(phi_vec_fn(C, t), dtype=float)  # shape (m,)

    # division operator: inflow - outflow
    div = r * (M @ T - T)  # (m,)
    dT = div * crowd - phi * T
    return dT


# Karyotype ODE Piecewise Simulator
def simulate_karyotype_ode_piecewise(T0_by_N, drug, t_span, r, Kcap, beta,
                                     N_min=30, N_max=90,
                                     C_fn=None, phi_vec_fn=None, f_param_fn=None,
                                     t_eval=None, include_points=None,
                                     max_step=None, atol=1e-9, rtol=1e-7,
                                     tlag=0.0, eps=1e-10, renormalize_M=True):
    """
    Matrix-based karyotype ODE with piecewise integration across dose times.
    """
    if C_fn is None:
        C_fn = lambda t: 0.0

    Ns, M = build_misseg_matrix(N_min, N_max, beta, renormalize=renormalize_M)

    if phi_vec_fn is None:
        if f_param_fn is None:
            raise ValueError("Provide phi_vec_fn or f_param_fn to construct drug kill.")
        phi_vec_fn = make_phi_vec_fn(Ns, drug, f_param_fn)

    # initial state vector
    T0 = np.zeros(len(Ns), dtype=float)
    for N, val in T0_by_N.items():
        if N_min <= int(N) <= N_max:
            T0[int(N) - N_min] = float(val)

    t0, t1 = float(t_span[0]), float(t_span[1])

    # dose breaks (reuse your existing dose_times())
    dtimes = dose_times(drug, (t0, t1), tlag=tlag)
    dtimes = dtimes[(dtimes > t0) & (dtimes < t1)]
    breaks = np.array([t0, *dtimes.tolist(), t1], dtype=float)

    add = []
    if include_points is not None:
        add.extend(include_points)
    if t_eval is not None:
        add.extend(t_eval)

    post_dose = np.minimum(dtimes + eps, t1) if dtimes.size else np.array([], dtype=float)
    required_global = np.unique(
        np.concatenate(([t0, t1],
                        np.asarray(add, dtype=float) if len(add) else [],
                        post_dose))
    )

    t_all = []
    Y_all = []
    T_init = T0.copy()

    for a, b in zip(breaks[:-1], breaks[1:]):
        if b <= a + 1e-14:
            continue

        seg_points = required_global[(required_global >= a) & (required_global <= b)]
        seg_core = np.array([a, min(a + eps, b), b], dtype=float)
        seg_eval = np.unique(np.concatenate((seg_points, seg_core)))

        sol = solve_ivp(
            rhs_karyotype_ode, (a, b), T_init,
            t_eval=seg_eval,
            args=(r, Kcap, M, phi_vec_fn, C_fn),
            rtol=rtol, atol=atol, max_step=max_step
        )

        if t_all and abs(sol.t[0] - t_all[-1]) < 1e-12:
            t_all.extend(sol.t[1:].tolist())
            Y_all.append(sol.y[:, 1:])
        else:
            t_all.extend(sol.t.tolist())
            Y_all.append(sol.y)

        T_init = sol.y[:, -1]

    t_all = np.array(t_all, dtype=float)
    T_mat = np.hstack(Y_all)
    T_total = T_mat.sum(axis=0)

    if t_eval is not None:
        t_req = np.asarray(t_eval, dtype=float)
        T_req = np.vstack([np.interp(t_req, t_all, T_mat[i, :]) for i in range(T_mat.shape[0])])
        return t_req, Ns, T_req, T_req.sum(axis=0), M

    return t_all, Ns, T_mat, T_total, M


if __name__ == "__main__":
    DRUG = "none"  # change drug here
    C_fn = get_concentration_curve(DRUG)  # your existing IV PK

    T_SPAN = (0.0, 30)
    DT = 0.05
    TIMES = build_times_with_doses(T_SPAN, DT, drug_name=DRUG, include_days=True, eps=1e-8)

    T0_by_N = {23: 3e3, 46: 3e3, 69: 3e3}

    beta = .05
    r = .3
    Kcap = 1e12

    t, Ns, T_mat, T_total, M = simulate_karyotype_ode_piecewise(
        T0_by_N, DRUG, t_span=T_SPAN, r=r, Kcap=Kcap, beta=beta,
        N_min=10, N_max=90,
        C_fn=C_fn,
        f_param_fn=f,
        t_eval=TIMES,
        max_step=0.05,
        renormalize_M=True
    )

    print(len(Ns))
    print(len(T_mat[:, -1]))
