# Imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import math
from functools import lru_cache
from scipy.special import comb
from scipy.stats import norm
plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
plt.rcParams["text.latex.preamble"] += r"\usepackage{amssymb}"

# Drug Sensitivities Functions


def phi_Hill(C, EC50, n, Emax=1.0):
    C = float(C)
    EC50 = max(float(EC50), 1e-12)
    n = max(float(n), 0.1)
    return Emax * (C**n) / (EC50**n + C**n)

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
    def clamp_ec50(x):  return max(x, 1e-12)
    def clamp_n(x):     return max(x, 0.1)

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
    
    elif drug == "none":
        n_out = clamp_n(28.92 * math.exp(-0.94 * ploidy) + 0.92)
        ec50  = clamp_ec50(0.004 * math.exp(0.78 * ploidy) - 0.01)
        return dict(EC50=ec50, n=n_out, Emax=1.0)
    else:
        raise ValueError(f"Unknown drug: {drug}")


# Dosing Functions

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
    "none"          : "IV",
    "bay1895344"    : "IV",
    "ispinesib"     : "IV",
    "navitoclax"    : "Oral",
    "adavosertib"   : "Oral",
}

IV_DEFAULTS = dict(C_peak=1.0, half_life=0.5, period=7.0)

ORAL_DEFAULTS = {
    "dose": 100.0, "F": 0.6, "Vd": 60.0,
    "ka_day": 2.0,
    "ke_day": 0.7,
    "period": 1.0,
    "tlag": 0.0
}

PER_DRUG = {
    "volasertib":   {"C_peak": 1.0, "half_life": 4.0, "period": 7.0},
    "alisertib":    {"C_peak": 42.5, "half_life": 0.875, "period": 0.5},
    "cytarabine":   {"C_peak": 1.0, "half_life": 0.2, "period": 3.5},
    "gemcitabine":  {"C_peak": 239, "half_life": 0.05, "period": 7.0},
    "ispinesib":    {"C_peak": 2.1, "half_life": 1.04, "period": 7},
    "umi-77":       {"C_peak": 1.0, "half_life": 0.8, "period": 7.0},
    "navitoclax":   {"C_peak": 1.0, "half_life": 0.73, "period": 1},
    "bay1895344":   {"C_peak": 6.2, "half_life": 0.50, "period": 0.5},
    "none":         {"C_peak": 0, "half_life": 0.50, "period": 7},
    "osi-027":      {"dose": 50, "F": 0.5, "Vd": 80, "ka_day": 1.8, "ke_day": 0.5, "period": 1.0},
    "abt-199":      {"dose": 100, "F": 0.6, "Vd": 250, "ka_day": 1.0, "ke_day": 0.3, "period": 1.0},
    "abt-263":      {"dose": 100, "F": 0.6, "Vd": 120, "ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "ceralasertib": {"dose": 80,  "F": 0.5, "Vd": 100,"ka_day": 2.0, "ke_day": 0.5, "period": 1.0},
    "adavosertib":  {"dose": 100, "F": 0.6, "Vd": 65, "ka_day": 2.4, "ke_day": 0.6, "period": 1.0},
    "tas":          {"dose": 60,  "F": 0.5, "Vd": 40, "ka_day": 2.4, "ke_day": 0.7, "period": 1.0},
    "tegafur":      {"dose": 40, "F": 0.5,"Vd": 45, "ka_day": 1.6, "ke_day": 0.5, "period": 1.0},
    "capecitabine": {"dose": 100, "F": 0.8, "Vd": 40, "ka_day": 3.0, "ke_day": 0.6, "period": 1.0},
    "5-azacytidine":{"dose": 100, "F": 0.2, "Vd": 40, "ka_day": 3.0, "ke_day": 2.0, "period": 1.0},
}

def get_concentration_curve(drug_name: str, **overrides):
    drug_name = drug_name.lower()
    route = drug_dosing_schedules.get(drug_name)
    if route is None:
        raise KeyError(f"Unknown drug: {drug_name}")

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
    else:
        period = (PER_DRUG.get(drug_name, ORAL_DEFAULTS)).get("period", ORAL_DEFAULTS["period"])
    t0, t1 = t_span
    start_k = math.ceil((t0 - tlag) / period)
    end_k   = math.floor((t1 - tlag) / period)
    return tlag + period * np.arange(start_k, end_k + 1)

def get_dose_period(drug_name: str) -> float:
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
    t0, t1 = float(t_span[0]), float(t_span[1])
    tau = float(get_dose_period(drug_name))

    base = np.arange(t0, t1 + 0.5*dt, dt, dtype=float)

    k_start = math.ceil((t0 - tlag) / tau)
    k_end   = math.floor((t1 - tlag) / tau)
    if k_end >= k_start:
        doses = (tlag + tau * np.arange(k_start, k_end + 1, dtype=float))
    else:
        doses = np.array([], dtype=float)

    days = np.arange(math.ceil(t0), math.floor(t1) + 1, dtype=float) if include_days else np.array([], dtype=float)

    all_times = np.concatenate([base, doses, np.minimum(doses , t1), days])
    times = np.unique(all_times)
    if times[0] > t0:   times = np.insert(times, 0, t0)
    if times[-1] < t1:  times = np.append(times, t1)
    return times

@lru_cache(maxsize=64)
def _build_ms_matrix_cached(
    N_min: int,
    N_max: int,
    beta: float,
    boundary: str,
    renormalize: bool,
    tail_truncation: bool,
    eps_tail: float,
):
    """Vectorised, cached builder for the expected-offspring MS matrix.

    Each chromosome independently contributes to the daughter count:
      0  with probability  beta/2   (misseg → other daughter)
      1  with probability  1-beta   (faithful segregation)
      2  with probability  beta/2   (misseg → this daughter)

    The daughter count for a mother with *q* chromosomes is the sum of *q*
    i.i.d. copies of the above — computed efficiently via FFT-based
    convolution (O(q log q) per column instead of the O(q²) triple loop).

    Results are cached on all arguments so repeated calls with the same
    (N_min, N_max, beta, …) return the same arrays in O(1).
    """
    if boundary not in ("drop", "absorbminmax"):
        raise ValueError("boundary must be 'drop' or 'absorbminmax'.")
    if not (0.0 <= beta <= 1.0):
        raise ValueError("beta must be in [0,1].")
    if tail_truncation and not (0.0 < eps_tail < 1.0):
        raise ValueError("eps_tail must be in (0,1).")

    Ns = np.arange(N_min, N_max + 1, dtype=int)
    m = len(Ns)
    P = np.zeros((m, m), dtype=float)

    z = norm.ppf(1.0 - 0.5 * eps_tail) if tail_truncation else None

    # Single-chromosome PMF: values {0, 1, 2}
    p1 = np.array([beta / 2.0, 1.0 - beta, beta / 2.0])

    absorb = boundary == "absorbminmax"

    for j, q in enumerate(Ns):
        # PMF of daughter chromosome count via FFT exponentiation.
        # Support: 0 … 2q  (length 2q+1).
        n_fft = 2 * q + 1
        pmf = np.real(np.fft.ifft(np.fft.fft(p1, n=n_fft) ** q))
        pmf = np.maximum(pmf, 0.0)          # remove tiny negative noise

        # ── Tail truncation ────────────────────────────────────────────
        if tail_truncation:
            Tq = min(float(q), z * math.sqrt(float(q) * float(beta)))
            Tq_int = int(math.floor(Tq))
            lo_keep = max(0, q - Tq_int)
            hi_keep = min(2 * q, q + Tq_int)
            if lo_keep > 0:
                pmf[:lo_keep] = 0.0
            if hi_keep < 2 * q:
                pmf[hi_keep + 1:] = 0.0

        # ── Map onto [N_min … N_max] grid ─────────────────────────────
        col = np.zeros(m, dtype=float)

        # Overlap between pmf support [0..2q] and grid [N_min..N_max]
        src_lo = max(N_min, 0)
        src_hi = min(N_max, 2 * q)
        if src_hi >= src_lo:
            col[src_lo - N_min: src_hi - N_min + 1] = pmf[src_lo: src_hi + 1]

        if absorb:
            # Accumulate out-of-range mass onto boundary bins
            if N_min > 0:
                col[0] += pmf[:N_min].sum()
            if N_max < 2 * q:
                col[-1] += pmf[N_max + 1:].sum()

        # ── Optional renormalisation ───────────────────────────────────
        if renormalize:
            s = col.sum()
            if s > 0:
                col /= s
            else:
                col[q - N_min] = 1.0

        P[:, j] = col

    B = 2.0 * P
    return Ns, B


def build_ms_expected_offspring_matrix(
    N_min: int,
    N_max: int,
    beta: float,
    *,
    boundary: str = "absorbminmax",
    renormalize: bool = False,
    tail_truncation: bool = True,
    eps_tail: float = 1e-8,
):
    """Build B (expected-offspring MS matrix).

    Thin wrapper that preserves the original keyword-only API while
    delegating to the cached, vectorised implementation.
    """
    return _build_ms_matrix_cached(
        N_min, N_max, beta, boundary, renormalize, tail_truncation, eps_tail,
    )


@lru_cache(maxsize=64)
def build_wgd_matrix(N_min: int, N_max: int, boundary: str = "absorbminmax"):
    """
    BW[i,j] = expected # of daughters in state Ns[i] from WGD division of mother Ns[j].
    (BW)_{p,q} = 2 if p=2q else 0.
    """
    if boundary not in ("drop", "absorbminmax"):
        raise ValueError("boundary must be 'drop' or 'absorbminmax'.")
    Ns = np.arange(N_min, N_max + 1, dtype=int)
    m = len(Ns)
    BW = np.zeros((m, m), dtype=float)

    for j, q in enumerate(Ns):
        p = 2 * q
        if boundary == "absorbminmax":
            p = max(N_min, min(p, N_max))
            BW[p - N_min, j] += 2.0
        else:
            if N_min <= p <= N_max:
                BW[p - N_min, j] += 2.0

    return Ns, BW


# Karyotype ODE

def rhs_karyotype_ode(
    t, Y,
    r, Kcap, p_wgd, B0, B1, BW,
    damage_vec_fn,  
    C_fn,
    n_tr: int,
    k_tr: float,
    k_kill: float,
):
    """
    State layout:
      Y = [ x0 (m), x1 (m), Z1 (m), Z2 (m), ..., Zn (m) ]
    where Z_i are transit-chain stages driven by damage input D(t).
    Kill hazard: kappa(t) = k_kill * Z_n(t).
    """
    m = B0.shape[0]
    x0 = Y[0:m]
    x1 = Y[m:2*m]

    Z = Y[2*m:].reshape((n_tr, m))  

    xsum = float(np.sum(x0) + np.sum(x1))
    crowd = max(0.0, 1.0 - xsum / Kcap)

    C = float(C_fn(t))
    D = np.asarray(damage_vec_fn(C, t), dtype=float)  # (m,)

    # Transit chain dynamics
    dZ = np.zeros_like(Z)
    dZ[0, :] = k_tr * (D - Z[0, :])
    for i in range(1, n_tr):
        dZ[i, :] = k_tr * (Z[i-1, :] - Z[i, :])

    # Commitment -> kill hazard
    kappa = k_kill * Z[-1, :]  # (m,)

    # Population dynamics with delayed kill
    dx0 = (r * (1.0 - p_wgd) * (B0 @ x0) - r * x0) * crowd - kappa * x0
    dx1 = (r * (B1 @ x1) - r * x1 + r * p_wgd * (BW @ x0)) * crowd - kappa * x1

    return np.concatenate([dx0, dx1, dZ.reshape(-1)])

def simulate_karyotype_ode_piecewise(
    T0_by_N, drug, t_span, r, Kcap, beta,
    N_min=30, N_max=90,
    p_wgd=0.0,
    boundary="drop",
    C_fn=None,
    f_param_fn=None,
    t_eval=None, include_points=None,
    max_step=None, atol=1e-9, rtol=1e-7,
    tlag=0.0, eps=1e-10, renormalize_M=False,
    # --- transit-chain kill params ---
    n_tr: int = 3,
    k_tr: float = 1.0,
    k_kill: float = 1.0,
):
    """
    Kill model:
      D_N(t) = Hill(C(t); params(N))
      Z chain: dZ1 = k_tr(D - Z1), dZi = k_tr(Z_{i-1} - Zi)
      kappa_N(t) = k_kill * Z_n(t)
      death term: -kappa_N(t) * x_N(t)

    Notes:
      - Setting k_kill=1 makes steady-state kappa match instantaneous Hill kill.
      - n_tr controls the shape/spread (gamma distributed delay).
      - k_tr sets timescale (mean delay ~ n_tr / k_tr).
    """
    
    if C_fn is None:
        raise ValueError("Concentration function C_fn must be provided.")
    if not (0.0 <= p_wgd <= 1.0):
        raise ValueError("p_wgd must be in [0,1].")
    if n_tr < 1:
        raise ValueError("n_tr must be >= 1.")
    if k_tr <= 0.0:
        raise ValueError("k_tr must be > 0.")
    if k_kill < 0.0:
        raise ValueError("k_kill must be >= 0.")
    if f_param_fn is None:
        raise ValueError("Provide f_param_fn to construct Hill damage input.")

    
    Ns, B0 = build_ms_expected_offspring_matrix(
        N_min, N_max, beta, boundary=boundary, renormalize=renormalize_M
    )
    B1 = B0                                       # B1 uses identical parameters
    _, BW = build_wgd_matrix(N_min, N_max, boundary=boundary)

    damage_vec_fn = make_phi_vec_fn(Ns, drug, f_param_fn)

    m = len(Ns)

    # Initial populations
    x0_0 = np.zeros(m, dtype=float)
    x1_0 = np.zeros(m, dtype=float)
    for N, val in T0_by_N.items():
        if N_min <= int(N) <= N_max:
            x0_0[int(N) - N_min] = float(val)

    # Initial transit-chain stages (start at 0 damage/commitment)
    Z0 = np.zeros((n_tr, m), dtype=float)

    Y_init = np.concatenate([x0_0, x1_0, Z0.reshape(-1)])

    t0, t1 = float(t_span[0]), float(t_span[1])

    # Piecewise breaks at dosing times
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

    for a, b in zip(breaks[:-1], breaks[1:]):
        if b <= a + 1e-14:
            continue

        seg_points = required_global[(required_global >= a) & (required_global <= b)]
        seg_core = np.array([a, min(a + eps, b), b], dtype=float)
        seg_eval = np.unique(np.concatenate((seg_points, seg_core)))

        sol = solve_ivp(
            rhs_karyotype_ode, (a, b), Y_init,
            method='LSODA',
            t_eval=seg_eval,
            args=(r, Kcap, p_wgd, B0, B1, BW, damage_vec_fn, C_fn, n_tr, k_tr, k_kill),
            rtol=rtol, atol=atol, max_step=max_step
        )

        if t_all and abs(sol.t[0] - t_all[-1]) < 1e-12:
            t_all.extend(sol.t[1:].tolist())
            Y_all.append(sol.y[:, 1:])
        else:
            t_all.extend(sol.t.tolist())
            Y_all.append(sol.y)

        Y_init = sol.y[:, -1]

    t_all = np.array(t_all, dtype=float)
    Y_mat = np.hstack(Y_all)  

    x0_mat = Y_mat[:m, :]
    x1_mat = Y_mat[m:2*m, :]
    Z_mat  = Y_mat[2*m:, :].reshape((n_tr, m, Y_mat.shape[1]))  

    x_total = x0_mat.sum(axis=0) + x1_mat.sum(axis=0)

    if t_eval is not None:
        t_req = np.asarray(t_eval, dtype=float)
        Y_req = np.vstack([np.interp(t_req, t_all, Y_mat[i, :]) for i in range(Y_mat.shape[0])])

        x0_req = Y_req[:m, :]
        x1_req = Y_req[m:2*m, :]
        Z_req  = Y_req[2*m:, :].reshape((n_tr, m, t_req.size))

        x_tot_req = x0_req.sum(axis=0) + x1_req.sum(axis=0)
        return t_req, Ns, x0_req, x1_req, x_tot_req, B0, B1, BW, Z_req

    return t_all, Ns, x0_mat, x1_mat, x_total, B0, B1, BW, Z_mat


# Plotting Functions

def plot_total_burden(t, T_total, title="Total tumor burden"):
    plt.figure()
    plt.plot(t, T_total)
    plt.xlabel("Time (days)")
    plt.ylabel("Total cells")
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_selected_states(t, Ns, T_mat, states, title="Selected karyotype trajectories"):
    plt.figure()
    for N in states:
        if N < Ns[0] or N > Ns[-1]:
            continue
        plt.plot(t, T_mat[N - Ns[0], :], label=f"N={N}")
    plt.xlabel("Time (days)")
    plt.ylabel("Cells in state")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plot_karyotype_heatmap(t, Ns, T_mat, title="Karyotype distribution over time"):
    plt.figure()
    Z = np.log1p(T_mat)
    plt.imshow(
        Z,
        aspect="auto",
        origin="lower",
        extent=[t[0], t[-1], Ns[0], Ns[-1]]
    )
    plt.colorbar(label=r"$\log(1 + T_N)$")
    plt.xlabel("Time (days)")
    plt.ylabel("Chromosome number N")
    plt.title(title)
    plt.tight_layout()
    plt.show()
    
def plot_dosing_curve(drug_name, t_span=(0.0, 14.0), n_points=2000, tlag=0.0,
                      mark_doses=True, overrides=None, title=None):
    overrides = overrides or {}
    dn = drug_name.lower()
    C_fn = get_concentration_curve(dn, **overrides)

    t0, t1 = float(t_span[0]), float(t_span[1])
    t = np.linspace(t0, t1, int(n_points))
    C = C_fn(t)

    plt.figure()
    plt.plot(t, C, label=f"{dn}")

    if mark_doses:
        try:
            dt = dose_times(dn, (t0, t1), tlag=tlag)
            dt = dt[(dt >= t0) & (dt <= t1)]
            ymax = float(np.nanmax(C)) if np.size(C) else 1.0
            plt.vlines(dt, ymin=0.0, ymax=ymax * 1.02, label="dose times")
        except Exception:
            pass

    plt.xlabel("Time (days)")
    plt.ylabel("Concentration C(t)")
    if title is None:
        route = drug_dosing_schedules.get(dn, "")
        title = f"Dosing curve: {dn}" + (f" ({route})" if route else "")
    plt.title(title)
    plt.tight_layout()
    plt.show()
    
def weighted_avg_karyotype(t, Ns, T_mat, eps=1e-30):
    Ns = np.asarray(Ns, dtype=float)
    T = np.asarray(T_mat, dtype=float)
    total = T.sum(axis=0)
    meanN = (Ns[:, None] * T).sum(axis=0) / (total + eps)
    return meanN, total

def plot_weighted_avg_karyotype(t, Ns, T_mat, title="Weighted average karyotype"):
    meanN, total = weighted_avg_karyotype(t, Ns, T_mat)
    plt.figure()
    plt.plot(t, meanN)
    plt.xlabel("Time (days)")
    plt.ylabel("Weighted average chromosome number  N̄(t)")
    plt.ylim(10, 90)
    plt.title(title)
    plt.tight_layout()
    plt.show()
    
def plot_pre_post_wgd_totals(t, x0_mat, x1_mat, title="Pre- vs Post-WGD burden"):
    plt.figure()
    plt.plot(t, x0_mat.sum(axis=0), label="Pre-WGD")
    plt.plot(t, x1_mat.sum(axis=0), label="Post-WGD")
    plt.plot(t, (x0_mat + x1_mat).sum(axis=0), "--", label="Total")
    plt.xlabel("Time (days)")
    plt.ylabel("Cells")
    plt.legend()
    plt.title(title)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # Sample Usage 
    DRUG =  "ispinesib"  # "alisertib", "gemcitabine", "bay1895344", "none"
    C_fn = get_concentration_curve(DRUG)

    T_SPAN = (0.0, 30)
    DT = 0.05
    TIMES = build_times_with_doses(T_SPAN, DT, drug_name=DRUG, include_days=True, eps=1e-8)

    T0_by_N = {23: 3e3, 46: 3e3, 69: 3e3} 

    beta = 0
    r = .5
    Kcap = 1e12
    pWGD = 0

    t, Ns, x0, x1, xtot, B0, B1, BW, Z = simulate_karyotype_ode_piecewise(
        T0_by_N, DRUG, t_span=T_SPAN, r=r, Kcap=Kcap, beta=beta,
        N_min=10, N_max=150,
        p_wgd=pWGD,
        boundary="drop",  #absorbminmax or drop
        C_fn=C_fn,
        f_param_fn=f,
        t_eval=TIMES,
        max_step=0.05,
        renormalize_M=True, 
        n_tr=3,
        k_tr=1.5,
        k_kill=1
    )

    # Combine compartments for karyotype-resolved plots
    x_combined = x0 + x1

    # Plot
    plot_total_burden(t, xtot, title=f"Total burden (drug={DRUG}, beta={beta}, p_WGD={pWGD})")
    plot_selected_states(t, Ns, x_combined, states=[23, 46, 69], title="Selected N trajectories (combined)")
    plot_karyotype_heatmap(t, Ns, x_combined, title="Karyotype heatmap (combined, log scale)")
    plot_weighted_avg_karyotype(t, Ns, x_combined, title="Mean karyotype vs time (combined)")
    plot_total_burden(t, x0.sum(axis=0), title="Pre-WGD total burden")
    plot_total_burden(t, x1.sum(axis=0), title="Post-WGD total burden")
