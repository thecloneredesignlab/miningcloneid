from __future__ import annotations

import math
import os
from multiprocessing import Pool
from time import time

import numpy as np
import torch
import importlib

N_CHAINS       = 32
N_BURNIN       = 6000
N_SAMPLES      = 14000
THIN           = 2
ADAPT_INTERVAL = 50
TARGET_ACCEPT  = 0.234  # https://www.sciencedirect.com/science/article/pii/S0304414919306982#sec3

# ── In-vivo live/dead model parameters (from invitro_fitting.py) ──────────────
# All 9 biological parameters from Table 1, plus log-noise sigma_B.
#
# `r` (net growth rate) is the one exception to log-space fitting: it is fit
# directly in natural, SIGNED units so the sampler can reach negative values
# (net tumor regression independent of the kill term), which a log
# parameterization can never represent since exp(x) > 0 for all x.  Every
# other biological parameter must stay strictly positive and is still fit in
# log-space as before.
#
# Every parameter — including `r` — carries a hard bound (PARAM_BOUNDS,
# below).  Bounds are wide enough to comfortably contain the posterior, but
# finite: without them, unconstrained MCMC proposals can wander into regions
# where the ODE (e.g. huge K combined with huge k_cyto/k_kill) becomes
# extremely stiff, and LSODA's internal step count — and therefore wall-clock
# time per iteration — blows up, which is what caused fitting to stall.
#
# Init values use the 2N in-vitro MAP as a warm start; MCMC will move them to
# the in-vivo posterior.
INIT_VALS = {
    "r":              1.25,
    # K_init = 1e8 cells (in-vivo scale).  The original value of 25358 came
    # from an in-vitro calibration where cell counts were much smaller; at
    # in-vivo scale that K lies below MIN_SIZE (1e5), making every simulation
    # return None and pinning all chains at energy = 1e12.  1e8 keeps
    # alive* = K×(1 − μ/r) ≫ MIN_SIZE.
    "log_K":          math.log(1e8),
    "log_k_kill":     math.log(217.3),
    "log_k_clear":    math.log(0.152),
    "log_k_cyto":     math.log(8041.0),
    "log_beta_dose":  math.log(0.253),
    "log_mu_base":    math.log(0.117),
    "log_mu_conf":    math.log(0.051),
    "log_sigma_B":    math.log(0.5),
}

INIT_STEP = {
    "r":              0.15,   # natural-space step (≈12% of the r=1.25 init)
    "log_K":          0.15,
    "log_k_kill":     0.20,
    "log_k_clear":    0.15,
    "log_k_cyto":     0.20,
    "log_beta_dose":  0.10,
    "log_mu_base":    0.15,
    "log_mu_conf":    0.15,
    "log_sigma_B":    0.05,
}

# Hard bounds on every fitted base parameter.  For the 8 log-space
# parameters these bound the LOG value (i.e. natural value is bounded to
# [exp(lo), exp(hi)]); for `r` they bound the natural, signed value directly.
# These are deliberately generous — wide enough that the bound should rarely
# bind for a well-fit mouse — but finite, which is what keeps MCMC proposals
# (and the resulting ODE solves) numerically tractable.
PARAM_BOUNDS = {
    # r: allow net regression (negative) as well as strong net growth. ±5/day
    # is already a biologically extreme daily rate, so this comfortably spans
    # anything a real tumor should show while still bounding the ODE.
    "r":              (-5.0,            5.0),
    "log_K":          (math.log(1e2),   math.log(1e10)),
    "log_k_kill":     (math.log(1e-5),  math.log(1e7)),
    "log_k_clear":    (math.log(1e-4),  math.log(20.0)),
    "log_k_cyto":     (math.log(1e-4),  math.log(1e7)),
    "log_beta_dose":  (math.log(0.05),  math.log(3.0)),
    "log_mu_base":    (math.log(1e-5),  math.log(2.0)),
    "log_mu_conf":    (math.log(1e-5),  math.log(2.0)),
    "log_sigma_B":    (math.log(0.01),  math.log(5.0)),
}

# ── Joint (hierarchical) model constants ──────────────────────────────────────
# 8 biological parameters receive a per-mouse epsilon offset in the same
# space they're fit in: additive in log-space for the 7 log parameters (i.e.
# a *multiplicative* deviation in natural units), and additive directly in
# natural units for `r` (since `r` has no log transform to be additive in).
# sigma_B and PK params remain global (shared across all mice).
#
#   param_m = param_global + epsilon_m        (in the parameter's own fit-space)
#   epsilon_m   ~ N(0, EPSILON_PRIOR_STD²), truncated to ±3σ (see eps bounds
#                 in run_mcmc_joint) for the same stiffness reasons as above.
#
# EPSILON_PRIOR_STD = 0.2        means each mouse may deviate ≈ ±22 % from the
N_BIO_PARAMS     = 8
BIO_PARAM_NAMES  = [
    "r", "log_K", "log_k_kill", "log_k_clear",
    "log_k_cyto", "log_beta_dose", "log_mu_base", "log_mu_conf",
]
# Display-only natural-unit labels, one per entry in BIO_PARAM_NAMES (e.g.
# "log_k_kill" -> "k_kill").  Kept as an explicit list rather than
# string-stripping a "log_" prefix so it stays correct for entries like "r"
# that have no prefix to strip.
BIO_PARAM_LABELS = ["r", "K", "k_kill", "k_clear",
                     "k_cyto", "beta_dose", "mu_base", "mu_conf"]
EPSILON_PRIOR_STD = 0.2   # prior σ on per-mouse epsilons
EPSILON_STEP      = 0.05  # initial MCMC proposal step size for epsilons

_W = {}


# ─────────────────────────────────────────────────────────────────────────────
# Single-mouse MCMC (original, kept for backward compatibility)
# ─────────────────────────────────────────────────────────────────────────────

def _init_worker(caller_module_name, state_dict):
    caller = importlib.import_module(caller_module_name)

    _W["obs"]            = state_dict["obs"]
    _W["end_ploidy"]     = state_dict["end_ploidy"]
    _W["drugs"]          = state_dict["drugs"]
    _W["pk_cfg"]         = state_dict["pk_cfg"]
    _W["sample_name"]    = state_dict["sample_name"]
    _W["initial_ploidy"] = state_dict["initial_ploidy"]

    # Override the caller module's globals so its functions use our data
    caller.INITIAL_PLOIDY                   = state_dict["initial_ploidy"]
    caller.OBSERVED_TUMOR_BURDENS           = state_dict["obs"]
    caller.OBSERVED_DRUGS_ADMINISTERED      = state_dict["drugs"]
    caller.OBSERVED_END_PLOIDY_DISTRIBUTION = state_dict["end_ploidy"]

    _W["fn_combined"] = caller._simulate_burden_and_ploidy


# Objective function; lower is better
def _log_posterior_single(param_vec, pk_param_names, pk_param_drugs):
    (r, log_K, log_k_kill, log_k_clear,
     log_k_cyto, log_beta_dose, log_mu_base, log_mu_conf,
     log_sB) = param_vec[:9]

    K          = math.exp(log_K)
    k_kill     = math.exp(log_k_kill)
    k_clear    = math.exp(log_k_clear)
    k_cyto     = math.exp(log_k_cyto)
    beta_dose  = math.exp(log_beta_dose)
    mu_base    = math.exp(log_mu_base)
    mu_conf    = math.exp(log_mu_conf)
    sigma_B    = math.exp(log_sB)

    # PK params follow the 9 base params (indices 9, 10, …)
    pk_state = {}
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_state.setdefault(drug, {})[pname] = math.exp(param_vec[9 + i])

    obs        = _W["obs"]
    obs_drugs  = _W["drugs"]
    pk_cfg     = _W["pk_cfg"]

    # ── Single ODE pass: get burden timeline ──────────────────────────────────
    nll = 0.0
    n_burden = 0

    need_burden = len(obs) > 1

    timeline = None
    if need_burden:
        timeline, _pred_ploidy = _W["fn_combined"](
            obs_drugs, r, K, k_kill, k_clear, k_cyto,
            beta_dose, mu_base, mu_conf,
            pk_state=pk_state)
        if timeline is None:
            return 1e12

    # L_B = Σ (log σ_B + 0.5 log(2π) + (log|v_obs| - log|v_hat|)² / (2σ_B²))
    if need_burden and timeline is not None:
        pred = {round(t, 2): v for t, v in timeline}
        obs_days = sorted(d for d in obs if d > 0)
        if obs_days:
            n_burden = len(obs_days)
            burden_nll = 0.0
            for d in obs_days:
                v_obs = obs[d]
                closest_t = min(pred.keys(), key=lambda t: abs(t - d))
                v_pred = pred[closest_t]
                if v_pred <= 0 or v_obs <= 0:
                    return 1e12
                residual = (math.log(v_obs) - math.log(v_pred)) ** 2
                burden_nll += (
                    math.log(sigma_B)
                    + 0.5 * math.log(2 * math.pi)
                    + residual / (2 * sigma_B ** 2)
                )
            nll += burden_nll

    # PK priors
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        cfg = pk_cfg[drug][pname]
        nll += 0.5 * ((param_vec[9 + i] - cfg["prior_log_mean"]) / cfg["prior_log_std"]) ** 2

    return nll


def _eval_worker(args):
    ci, pv, names, drugs = args
    return ci, _log_posterior_single(pv, names, drugs)


def run_mcmc(
    initial_ploidy, observed_burdens, observed_end_ploidy, observed_drugs,
    pk_params_to_fit, haploid_n, sample_name,
    fn_get_end_ploidy=None, fn_simulate_burden=None, fn_fill_gaps=None,
    caller_module_name="beam_search_flip_rate_wgd", verbose=True,
):
    # Setup
    base_names = [
        "r", "log_K", "log_k_kill", "log_k_clear",
        "log_k_cyto", "log_beta_dose", "log_mu_base", "log_mu_conf",
        "log_sigma_B",
    ]
    pk_param_names, pk_param_drugs = [], []
    pk_init_vals, pk_init_steps, pk_bounds = [], [], []

    obs_drug_set = {drug.lower() for _, _, drug, _ in observed_drugs}
    for drug, params in pk_params_to_fit.items():
        if drug.lower() not in obs_drug_set:
            continue
        for pname, cfg in params.items():
            pk_param_names.append(pname)
            pk_param_drugs.append(drug)
            pk_init_vals.append(math.log(cfg["init"]))
            pk_init_steps.append(cfg.get("step", 0.1))
            pk_bounds.append((math.log(cfg["init"]) - 3 * cfg["prior_log_std"],
                              math.log(cfg["init"]) + 3 * cfg["prior_log_std"]))

    n_base = len(base_names)
    n_dim  = n_base + len(pk_param_names)

    if verbose:
        print(f"\n{'='*60}\nMCMC Configuration\n{'='*60}")
        print(f"  Dimensions: {n_dim}  Chains: {N_CHAINS}  "
              f"Burn-in: {N_BURNIN}  Samples: {N_SAMPLES}  Thin: {THIN}")
        print(f"  Base params: {base_names}")
        print(f"  PK drugs fitted: {sorted(set(pk_param_drugs))}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if verbose:
        print(f"  Device: {device}")

    # Initialization.  Hard bounds (lo/hi) are enforced on init and on every
    # proposal below — see the PARAM_BOUNDS module docstring for why.
    init_vec = np.array([INIT_VALS[n] for n in base_names] + pk_init_vals, dtype=np.float64)
    step_vec = np.array([INIT_STEP[n] for n in base_names] + pk_init_steps, dtype=np.float64)
    lo_vec   = np.array([PARAM_BOUNDS[n][0] for n in base_names] + [b[0] for b in pk_bounds], dtype=np.float64)
    hi_vec   = np.array([PARAM_BOUNDS[n][1] for n in base_names] + [b[1] for b in pk_bounds], dtype=np.float64)

    rng = np.random.default_rng(42)
    chains_np = (np.tile(init_vec, (N_CHAINS, 1))
                 + rng.normal(0, step_vec * 0.1, size=(N_CHAINS, n_dim)))
    chains_np = np.clip(chains_np, lo_vec, hi_vec)

    chains = torch.tensor(chains_np, device=device, dtype=torch.float64)
    steps  = torch.tensor(step_vec,  device=device, dtype=torch.float64)
    lo_t   = torch.tensor(lo_vec,    device=device, dtype=torch.float64)
    hi_t   = torch.tensor(hi_vec,    device=device, dtype=torch.float64)

    worker_state = {
        "initial_ploidy": dict(initial_ploidy),
        "obs":            dict(observed_burdens),
        "end_ploidy":     np.array(observed_end_ploidy),
        "drugs":          list(observed_drugs),
        "pk_cfg":         dict(pk_params_to_fit),
        "sample_name":    sample_name,
    }

    max_from_env = int(os.environ.get("HARVEST_MAX_WORKERS", 0))
    n_workers = min(N_CHAINS, max_from_env or (os.cpu_count() or 4))

    pool = Pool(n_workers, initializer=_init_worker,
                initargs=(caller_module_name, worker_state))

    def _eval_all(chain_t):
        c_np = chain_t.cpu().numpy()
        args = [(ci, c_np[ci], pk_param_names, pk_param_drugs) for ci in range(N_CHAINS)]
        results = [None] * N_CHAINS
        for ci, nll in pool.map(_eval_worker, args):
            results[ci] = nll
        return torch.tensor(results, device=device, dtype=torch.float64)

    if verbose:
        print(f"\nEvaluating initial energies...")
    t0 = time()
    energies = _eval_all(chains)
    if verbose:
        print(f"  Done in {time()-t0:.1f}s  (best: {energies.min().item():.2f})")

    total_iters  = N_BURNIN + N_SAMPLES
    accept_count = torch.zeros(N_CHAINS, device=device, dtype=torch.float64)
    all_trace, post_samples = [], []
    best_energy = energies.min().item()
    best_params = chains[energies.argmin().item()].cpu().numpy().copy()

    if verbose:
        print(f"\nRunning MCMC ({total_iters} iters × {N_CHAINS} chains)...")

    for it in range(total_iters):
        is_burnin = it < N_BURNIN
        it_t0 = time()

        proposals     = chains + torch.randn_like(chains) * steps.unsqueeze(0)
        proposals     = torch.clamp(proposals, lo_t, hi_t)
        prop_energies = _eval_all(proposals)

        log_alpha  = energies - prop_energies
        accept_mask = torch.rand(N_CHAINS, device=device, dtype=torch.float64).log() < log_alpha
        chains   = torch.where(accept_mask.unsqueeze(1), proposals, chains)
        energies = torch.where(accept_mask, prop_energies, energies)
        accept_count += accept_mask.float()

        ci_best = energies.argmin().item()
        if energies[ci_best].item() < best_energy:
            best_energy = energies[ci_best].item()
            best_params = chains[ci_best].cpu().numpy().copy()

        if is_burnin and (it + 1) % ADAPT_INTERVAL == 0:
            mr = (accept_count / (it + 1)).mean().item()
            if mr < TARGET_ACCEPT - 0.05:
                steps *= 0.8
            elif mr > TARGET_ACCEPT + 0.05:
                steps *= 1.2
            if verbose and (it + 1) % 100 == 0:
                print(f"  {it+1:4d} burn-in | acc={mr:.3f} | best={best_energy:.2f} | {time()-it_t0:.1f}s")

        # Track MAP trace
        p_map = best_params.copy()
        pk_map_trace = {}
        for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
            pk_map_trace.setdefault(dr, {})[pn] = math.exp(p_map[n_base + i])
        all_trace.append({
            "iter": it, "burnin": is_burnin,
            "r":         p_map[0],
            "K":         math.exp(p_map[1]),
            "k_kill":    math.exp(p_map[2]),
            "k_clear":   math.exp(p_map[3]),
            "k_cyto":    math.exp(p_map[4]),
            "beta_dose": math.exp(p_map[5]),
            "mu_base":   math.exp(p_map[6]),
            "mu_conf":   math.exp(p_map[7]),
            "pk": pk_map_trace,
            "logP": -best_energy,
        })

        if not is_burnin and (it - N_BURNIN) % THIN == 0:
            c_all = chains.cpu().numpy()
            for ci in range(N_CHAINS):
                pv = c_all[ci]
                pkd = {}
                for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
                    pkd.setdefault(dr, {})[pn] = math.exp(pv[n_base + i])
                post_samples.append({
                    "r":          pv[0],
                    "K":          math.exp(pv[1]),
                    "k_kill":     math.exp(pv[2]),
                    "k_clear":    math.exp(pv[3]),
                    "k_cyto":     math.exp(pv[4]),
                    "beta_dose":  math.exp(pv[5]),
                    "mu_base":    math.exp(pv[6]),
                    "mu_conf":    math.exp(pv[7]),
                    "sigma_B":    math.exp(pv[8]),
                    "pk": pkd,
                    "energy": energies[ci].item(),
                })

        if verbose and not is_burnin and (it + 1) % 100 == 0:
            mr = (accept_count / (it + 1)).mean().item()
            print(f"  {it+1:4d} sample | acc={mr:.3f} | best={best_energy:.2f} | "
                  f"n={len(post_samples)} | {time()-it_t0:.1f}s")

    # MAP
    pk_map = {}
    for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_map.setdefault(dr, {})[pn] = math.exp(best_params[n_base + i])
    for drug, params in pk_params_to_fit.items():
        if drug not in pk_map:
            pk_map[drug] = {p: c["init"] for p, c in params.items()}

    map_params = {
        "r":          best_params[0],
        "K":          math.exp(best_params[1]),
        "k_kill":     math.exp(best_params[2]),
        "k_clear":    math.exp(best_params[3]),
        "k_cyto":     math.exp(best_params[4]),
        "beta_dose":  math.exp(best_params[5]),
        "mu_base":    math.exp(best_params[6]),
        "mu_conf":    math.exp(best_params[7]),
        "sigma_B":    math.exp(best_params[8]),
        "pk": pk_map,
        "energy": best_energy,
    }
    weights = np.ones(len(post_samples), dtype=float) / max(len(post_samples), 1)

    if verbose:
        print(f"\n{'='*60}\nMCMC Results\n{'='*60}")
        print(f"  MAP: r={map_params['r']:.4f}  K={map_params['K']:.3e}  "
              f"k_kill={map_params['k_kill']:.3e}  β={map_params['beta_dose']:.4f}")
        print(f"  mu_base={map_params['mu_base']:.4f}  mu_conf={map_params['mu_conf']:.4f}  "
              f"k_cyto={map_params['k_cyto']:.3e}  k_clear={map_params['k_clear']:.4f}")
        print(f"  Energy: {best_energy:.4f}  Samples: {len(post_samples)}")
        if post_samples:
            for k in ("r", "K", "k_kill", "beta_dose"):
                vals = [s[k] for s in post_samples]
                print(f"  E[{k}] = {np.mean(vals):.6g} ± {np.std(vals):.4g}")

    pool.close()
    pool.join()

    return map_params, post_samples, weights, all_trace


# ─────────────────────────────────────────────────────────────────────────────
# Joint (hierarchical) MCMC — all mice fitted simultaneously
# ─────────────────────────────────────────────────────────────────────────────
#
# Parameter vector layout (n_mice mice, n_pk PK params):
#
#   [0 : N_BIO_PARAMS]                          global log bio params (8)
#   [N_BIO_PARAMS]                               log_sigma_B  (global, index 8)
#   [9 : 9 + n_mice * N_BIO_PARAMS]             per-mouse epsilons
#                                                  mouse m → indices
#                                                  [9 + m*8 : 9 + (m+1)*8]
#   [9 + n_mice * N_BIO_PARAMS : ...]            global PK log-params
#
# Per-mouse effective log param:
#   log_param_m = log_param_global + epsilon_m
#
# Epsilon prior:
#   epsilon_m,k ~ N(0, EPSILON_PRIOR_STD²)   (iid, independent across mice/params)
# ─────────────────────────────────────────────────────────────────────────────

def _init_worker_joint(caller_module_name, pk_cfg, mice_states,
                       epsilon_prior_std=None):
    """
    Worker initializer for joint hierarchical fitting.

    Stores ALL mice's data in the worker-global _W dict.  The log-posterior
    function loops over mice sequentially inside each evaluation.
    """
    caller = importlib.import_module(caller_module_name)
    _W["caller"] = caller
    _W["mice"]   = mice_states   # list of per-mouse state dicts
    _W["pk_cfg"] = pk_cfg        # global PK config (same for all mice)
    _W["fn_combined"] = caller._simulate_burden_and_ploidy
    # Use the runtime value supplied by run_mcmc_joint; fall back to the
    # module-level default so the worker always has a valid prior width.
    _W["epsilon_prior_std"] = (epsilon_prior_std
                               if epsilon_prior_std is not None
                               else EPSILON_PRIOR_STD)


def _log_posterior_joint(param_vec, n_mice, pk_param_names, pk_param_drugs):
    """
    Joint log-posterior (returned as NLL, lower = better).

    Each mouse's likelihood is evaluated with its own effective parameters
    (global + epsilon_m), and the epsilon prior is added on top.
    """
    eps_start = N_BIO_PARAMS + 1          # index 9
    pk_start  = eps_start + n_mice * N_BIO_PARAMS

    # Global parameters.  Index 0 (r) is natural-space; indices 1-7 are
    # log-space (see BIO_PARAM_NAMES / PARAM_BOUNDS).
    global_bio = param_vec[:N_BIO_PARAMS]
    log_sB     = param_vec[N_BIO_PARAMS]
    sigma_B    = math.exp(log_sB)

    # Global PK state
    pk_state: dict = {}
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_state.setdefault(drug, {})[pname] = math.exp(param_vec[pk_start + i])

    caller    = _W["caller"]
    mice_list = _W["mice"]
    pk_cfg    = _W["pk_cfg"]

    total_nll = 0.0

    # ── Epsilon prior  N(0, σ²) for every (mouse, param) pair ─────────────────
    _eps_std = _W.get("epsilon_prior_std", EPSILON_PRIOR_STD)
    inv_var  = 1.0 / (_eps_std ** 2)
    for m in range(n_mice):
        o     = eps_start + m * N_BIO_PARAMS
        eps_m = param_vec[o : o + N_BIO_PARAMS]
        total_nll += 0.5 * float(np.dot(eps_m, eps_m)) * inv_var

    # ── Per-mouse burden likelihoods ───────────────────────────────────────────
    for m, mouse in enumerate(mice_list):
        o     = eps_start + m * N_BIO_PARAMS
        eps_m = param_vec[o : o + N_BIO_PARAMS]

        # Effective bio params for this mouse (global + epsilon, each in its
        # own fit-space: natural for r, log for the rest)
        eff = global_bio + eps_m

        (r, log_K, log_k_kill, log_k_clear,
         log_k_cyto, log_beta_dose, log_mu_base, log_mu_conf) = eff

        K         = math.exp(log_K)
        k_kill    = math.exp(log_k_kill)
        k_clear   = math.exp(log_k_clear)
        k_cyto    = math.exp(log_k_cyto)
        beta_dose = math.exp(log_beta_dose)
        mu_base   = math.exp(log_mu_base)
        mu_conf   = math.exp(log_mu_conf)

        obs       = mouse["obs"]
        obs_drugs = mouse["drugs"]

        # Point the caller module's globals at this mouse's data so that
        # _simulate_burden_and_ploidy (and the ODE helpers it calls) use
        # the correct INITIAL_PLOIDY and FIRST_TX_DAY.
        caller.INITIAL_PLOIDY              = mouse["initial_ploidy"]
        caller.FIRST_TX_DAY                = mouse["first_tx_day"]
        caller.OBSERVED_TUMOR_BURDENS      = obs
        caller.OBSERVED_DRUGS_ADMINISTERED = obs_drugs

        if len(obs) <= 1:
            continue   # no burden data → skip likelihood for this mouse

        timeline, _ = caller._simulate_burden_and_ploidy(
            obs_drugs, r, K, k_kill, k_clear, k_cyto,
            beta_dose, mu_base, mu_conf,
            pk_state=pk_state,
        )
        if timeline is None:
            return 1e12

        pred     = {round(t, 2): v for t, v in timeline}
        obs_days = sorted(d for d in obs if d > 0)
        for d in obs_days:
            v_obs     = obs[d]
            closest_t = min(pred.keys(), key=lambda t: abs(t - d))
            v_pred    = pred[closest_t]
            if v_pred <= 0 or v_obs <= 0:
                return 1e12
            residual   = (math.log(v_obs) - math.log(v_pred)) ** 2
            total_nll += (
                math.log(sigma_B)
                + 0.5 * math.log(2 * math.pi)
                + residual / (2 * sigma_B ** 2)
            )

    # ── Global PK priors ───────────────────────────────────────────────────────
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        cfg        = pk_cfg[drug][pname]
        total_nll += 0.5 * ((param_vec[pk_start + i] - cfg["prior_log_mean"])
                             / cfg["prior_log_std"]) ** 2

    return total_nll


def _eval_worker_joint(args):
    ci, pv, n_mice, names, drugs = args
    return ci, _log_posterior_joint(pv, n_mice, names, drugs)


def run_mcmc_joint(
    mice_data: list[dict],
    pk_params_to_fit: dict,
    sample_names: list[str],
    caller_module_name: str = "beam_search_flip_rate_wgd",
    verbose: bool = True,
    epsilon_prior_std: float = EPSILON_PRIOR_STD,
) -> dict:
    """
    Hierarchical MCMC fitting of all mice simultaneously.

    All mice share global biological parameters.  Each mouse additionally has
    per-parameter log-space epsilon offsets so that individual variation is
    captured while all data jointly inform the global parameters.

    Parameters
    ----------
    mice_data : list of dicts, one per mouse.  Required keys:
        initial_ploidy  – dict[int, float]
        obs             – dict[int, float]  (burden observations, cells)
        end_ploidy      – np.ndarray        (may be empty)
        drugs           – list of (start, end, drug, dose) tuples
        first_tx_day    – float
    pk_params_to_fit : global PK configuration (same structure as in
        beam_search_flip_rate_wgd._DK["PK_PARAMS_TO_FIT"])
    sample_names : list of harvest names (same length as mice_data, for logging)

    Returns
    -------
    dict with keys:
        global_map      – MAP global biological params + sigma_B + pk
        mice_eps_map    – list of per-mouse MAP epsilon dicts (natural-space
                            offset for r; log-space offsets for the rest)
        mice_map_params – list of per-mouse effective MAP params (natural units)
        post_samples    – list of posterior sample dicts; each has:
                            "global"      : global natural-scale params
                            "epsilons"    : list of per-mouse epsilon arrays
                            "mice_params" : list of per-mouse natural-scale params
                            "sigma_B", "pk", "energy"
        weights         – uniform weight array (n_samples,)
        all_trace       – iteration-by-iteration MAP trace.  Each entry has
                            the global params (for reference/logging) AND a
                            "mice" list (len n_mice) of each mouse's own
                            EFFECTIVE params (global + epsilon) at that
                            iteration — use "mice"[m] to animate mouse m's
                            fit, never the global-only values.
        sample_names    – passed through for downstream use
    """
    n_mice = len(mice_data)
    if n_mice == 0:
        raise ValueError("mice_data is empty")

    # ── PK params: union of all drugs seen across all mice ─────────────────────
    all_drug_set: set[str] = set()
    for mouse in mice_data:
        all_drug_set.update(drug.lower() for _, _, drug, _ in mouse["drugs"])

    pk_param_names: list[str] = []
    pk_param_drugs: list[str] = []
    pk_init_vals:   list[float] = []
    pk_init_steps:  list[float] = []
    pk_bounds:      list[tuple] = []

    for drug, params in pk_params_to_fit.items():
        if drug.lower() not in all_drug_set:
            continue
        for pname, cfg in params.items():
            pk_param_names.append(pname)
            pk_param_drugs.append(drug)
            pk_init_vals.append(math.log(cfg["init"]))
            pk_init_steps.append(cfg.get("step", 0.1))
            pk_bounds.append((
                math.log(cfg["init"]) - 3 * cfg["prior_log_std"],
                math.log(cfg["init"]) + 3 * cfg["prior_log_std"],
            ))

    n_pk  = len(pk_param_names)
    eps_start = N_BIO_PARAMS + 1          # index 9
    pk_start  = eps_start + n_mice * N_BIO_PARAMS
    n_dim     = pk_start + n_pk

    if verbose:
        print(f"\n{'='*60}")
        print(f"Joint Hierarchical MCMC Configuration")
        print(f"{'='*60}")
        print(f"  Mice          : {n_mice}")
        print(f"  Total dims    : {n_dim}  (9 global + {n_mice}×{N_BIO_PARAMS} ε + {n_pk} PK)")
        print(f"  Chains        : {N_CHAINS}  Burn-in: {N_BURNIN}  Samples: {N_SAMPLES}")
        print(f"  ε prior σ     : {epsilon_prior_std:.3f} "
              f"(log-space for K/k_kill/.../mu_conf, ≈{100*(math.exp(epsilon_prior_std)-1):.0f}% natural; "
              f"natural-space (additive) for r)")
        print(f"  Global params : {BIO_PARAM_NAMES}")
        for i, name in enumerate(sample_names):
            print(f"  Mouse {i:2d}      : {name}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if verbose:
        print(f"  Device        : {device}")

    # ── Build init / step / bounds vectors ────────────────────────────────────
    # Layout: [global bio (8)] [log_sigma_B (1)] [epsilons (n_mice×8)] [pk (n_pk)]
    bio_init  = [INIT_VALS[n]          for n in BIO_PARAM_NAMES]
    bio_step  = [INIT_STEP[n]          for n in BIO_PARAM_NAMES]
    bio_lo    = [PARAM_BOUNDS[n][0]    for n in BIO_PARAM_NAMES]
    bio_hi    = [PARAM_BOUNDS[n][1]    for n in BIO_PARAM_NAMES]

    eps_lo   = -3.0 * epsilon_prior_std
    eps_hi   = +3.0 * epsilon_prior_std
    # Proposal step scales with prior width so it stays well inside the
    # support regardless of how tight the prior is set.  A fixed step of
    # 0.05 would wildly overshoot when epsilon_prior_std < ~0.05 and the
    # MCMC chain would oscillate only at the ±3σ clipping boundaries.
    eps_step =  epsilon_prior_std / 4.0

    init_vec = np.array(
        bio_init
        + [INIT_VALS["log_sigma_B"]]
        + [0.0] * (n_mice * N_BIO_PARAMS)      # epsilons start at 0
        + pk_init_vals,
        dtype=np.float64,
    )
    step_vec = np.array(
        bio_step
        + [INIT_STEP["log_sigma_B"]]
        + [eps_step] * (n_mice * N_BIO_PARAMS)
        + pk_init_steps,
        dtype=np.float64,
    )
    lo_vec = np.array(
        bio_lo
        + [PARAM_BOUNDS["log_sigma_B"][0]]
        + [eps_lo] * (n_mice * N_BIO_PARAMS)
        + [b[0] for b in pk_bounds],
        dtype=np.float64,
    )
    hi_vec = np.array(
        bio_hi
        + [PARAM_BOUNDS["log_sigma_B"][1]]
        + [eps_hi] * (n_mice * N_BIO_PARAMS)
        + [b[1] for b in pk_bounds],
        dtype=np.float64,
    )

    rng = np.random.default_rng(42)
    chains_np = (np.tile(init_vec, (N_CHAINS, 1))
                 + rng.normal(0, step_vec * 0.1, size=(N_CHAINS, n_dim)))
    chains_np = np.clip(chains_np, lo_vec, hi_vec)

    chains = torch.tensor(chains_np, device=device, dtype=torch.float64)
    steps  = torch.tensor(step_vec,  device=device, dtype=torch.float64)
    lo_t   = torch.tensor(lo_vec,    device=device, dtype=torch.float64)
    hi_t   = torch.tensor(hi_vec,    device=device, dtype=torch.float64)

    # Serialise per-mouse state for worker processes
    mice_states = [
        {
            "initial_ploidy": dict(mouse["initial_ploidy"]),
            "obs":            dict(mouse["obs"]),
            "end_ploidy":     np.array(mouse.get("end_ploidy", [])),
            "drugs":          list(mouse["drugs"]),
            "first_tx_day":   float(mouse["first_tx_day"]),
        }
        for mouse in mice_data
    ]

    max_from_env = int(os.environ.get("HARVEST_MAX_WORKERS", 0))
    n_workers    = min(N_CHAINS, max_from_env or (os.cpu_count() or 4))

    pool = Pool(
        n_workers,
        initializer=_init_worker_joint,
        initargs=(caller_module_name, pk_params_to_fit, mice_states, epsilon_prior_std),
    )

    def _eval_all(chain_t):
        c_np = chain_t.cpu().numpy()
        args = [(ci, c_np[ci], n_mice, pk_param_names, pk_param_drugs)
                for ci in range(N_CHAINS)]
        results = [None] * N_CHAINS
        for ci, nll in pool.map(_eval_worker_joint, args):
            results[ci] = nll
        return torch.tensor(results, device=device, dtype=torch.float64)

    if verbose:
        print(f"\nEvaluating initial joint energies...")
    t0 = time()
    energies = _eval_all(chains)
    if verbose:
        print(f"  Done in {time()-t0:.1f}s  (best: {energies.min().item():.2f})")

    total_iters  = N_BURNIN + N_SAMPLES
    accept_count = torch.zeros(N_CHAINS, device=device, dtype=torch.float64)
    all_trace, post_samples = [], []
    best_energy = energies.min().item()
    best_params = chains[energies.argmin().item()].cpu().numpy().copy()

    if verbose:
        print(f"\nRunning joint MCMC ({total_iters} iters × {N_CHAINS} chains)...")

    for it in range(total_iters):
        is_burnin = it < N_BURNIN
        it_t0 = time()

        proposals     = chains + torch.randn_like(chains) * steps.unsqueeze(0)
        proposals     = torch.clamp(proposals, lo_t, hi_t)
        prop_energies = _eval_all(proposals)

        log_alpha   = energies - prop_energies
        accept_mask = torch.rand(N_CHAINS, device=device, dtype=torch.float64).log() < log_alpha
        chains   = torch.where(accept_mask.unsqueeze(1), proposals, chains)
        energies = torch.where(accept_mask, prop_energies, energies)
        accept_count += accept_mask.float()

        ci_best = energies.argmin().item()
        if energies[ci_best].item() < best_energy:
            best_energy = energies[ci_best].item()
            best_params = chains[ci_best].cpu().numpy().copy()

        if is_burnin and (it + 1) % ADAPT_INTERVAL == 0:
            mr = (accept_count / (it + 1)).mean().item()
            if mr < TARGET_ACCEPT - 0.05:
                steps *= 0.8
            elif mr > TARGET_ACCEPT + 0.05:
                steps *= 1.2
            if verbose and (it + 1) % 100 == 0:
                print(f"  {it+1:4d} burn-in | acc={mr:.3f} | best={best_energy:.2f} | {time()-it_t0:.1f}s")

        # MAP trace.  Stores the GLOBAL params for reference/logging, plus a
        # per-mouse "mice" list with each mouse's EFFECTIVE params (global +
        # that mouse's epsilon) at this iteration.  Downstream animation code
        # must plot each mouse's own effective trajectory, not the shared
        # global-only one — otherwise every mouse's fitting GIF shows the
        # population-average curve, which can look badly wrong for any mouse
        # whose epsilon deviates from the global value, and which washes out
        # per-mouse dosing/kill dynamics (see beam_search_flip_rate_wgd.py).
        p     = best_params
        pk_tr = {}
        for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
            pk_tr.setdefault(dr, {})[pn] = math.exp(p[pk_start + i])

        mice_tr = []
        for m in range(n_mice):
            o     = eps_start + m * N_BIO_PARAMS
            eff   = p[:N_BIO_PARAMS] + p[o : o + N_BIO_PARAMS]
            mice_tr.append({
                "r":         eff[0],
                "K":         math.exp(eff[1]),
                "k_kill":    math.exp(eff[2]),
                "k_clear":   math.exp(eff[3]),
                "k_cyto":    math.exp(eff[4]),
                "beta_dose": math.exp(eff[5]),
                "mu_base":   math.exp(eff[6]),
                "mu_conf":   math.exp(eff[7]),
            })

        all_trace.append({
            "iter": it, "burnin": is_burnin,
            "r":         p[0],
            "K":         math.exp(p[1]),
            "k_kill":    math.exp(p[2]),
            "k_clear":   math.exp(p[3]),
            "k_cyto":    math.exp(p[4]),
            "beta_dose": math.exp(p[5]),
            "mu_base":   math.exp(p[6]),
            "mu_conf":   math.exp(p[7]),
            "pk": pk_tr,
            "mice": mice_tr,
            "logP": -best_energy,
        })

        # Posterior samples (stored thinned, after burn-in)
        if not is_burnin and (it - N_BURNIN) % THIN == 0:
            c_all = chains.cpu().numpy()
            for ci in range(N_CHAINS):
                pv = c_all[ci]

                # Global natural-scale params
                g = {
                    "r":         pv[0],
                    "K":         math.exp(pv[1]),
                    "k_kill":    math.exp(pv[2]),
                    "k_clear":   math.exp(pv[3]),
                    "k_cyto":    math.exp(pv[4]),
                    "beta_dose": math.exp(pv[5]),
                    "mu_base":   math.exp(pv[6]),
                    "mu_conf":   math.exp(pv[7]),
                }
                pkd = {}
                for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
                    pkd.setdefault(dr, {})[pn] = math.exp(pv[pk_start + i])

                # Per-mouse epsilons and effective params
                eps_list:   list[np.ndarray] = []
                mice_plist: list[dict]        = []
                for m in range(n_mice):
                    o     = eps_start + m * N_BIO_PARAMS
                    eps_m = pv[o : o + N_BIO_PARAMS].copy()
                    eff   = pv[:N_BIO_PARAMS] + eps_m
                    eps_list.append(eps_m)
                    mice_plist.append({
                        "r":         eff[0],
                        "K":         math.exp(eff[1]),
                        "k_kill":    math.exp(eff[2]),
                        "k_clear":   math.exp(eff[3]),
                        "k_cyto":    math.exp(eff[4]),
                        "beta_dose": math.exp(eff[5]),
                        "mu_base":   math.exp(eff[6]),
                        "mu_conf":   math.exp(eff[7]),
                        "sigma_B":   math.exp(pv[N_BIO_PARAMS]),
                        "pk":        pkd,
                    })

                post_samples.append({
                    "global":      g,
                    "epsilons":    eps_list,
                    "mice_params": mice_plist,
                    "sigma_B":     math.exp(pv[N_BIO_PARAMS]),
                    "pk":          pkd,
                    "energy":      energies[ci].item(),
                })

        if verbose and not is_burnin and (it + 1) % 100 == 0:
            mr = (accept_count / (it + 1)).mean().item()
            print(f"  {it+1:4d} sample | acc={mr:.3f} | best={best_energy:.2f} | "
                  f"n={len(post_samples)} | {time()-it_t0:.1f}s")

    # ── Extract MAP results ────────────────────────────────────────────────────
    p = best_params

    pk_map: dict = {}
    for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_map.setdefault(dr, {})[pn] = math.exp(p[pk_start + i])
    for drug, params in pk_params_to_fit.items():
        if drug not in pk_map:
            pk_map[drug] = {pr: c["init"] for pr, c in params.items()}

    global_map = {
        "r":         p[0],
        "K":         math.exp(p[1]),
        "k_kill":    math.exp(p[2]),
        "k_clear":   math.exp(p[3]),
        "k_cyto":    math.exp(p[4]),
        "beta_dose": math.exp(p[5]),
        "mu_base":   math.exp(p[6]),
        "mu_conf":   math.exp(p[7]),
        "sigma_B":   math.exp(p[N_BIO_PARAMS]),
        "pk":        pk_map,
        "energy":    best_energy,
    }

    mice_eps_map:    list[dict] = []
    mice_map_params: list[dict] = []
    for m in range(n_mice):
        o     = eps_start + m * N_BIO_PARAMS
        eps_m = p[o : o + N_BIO_PARAMS]
        eff   = p[:N_BIO_PARAMS] + eps_m
        mice_eps_map.append({
            BIO_PARAM_LABELS[k]: float(eps_m[k])
            for k in range(N_BIO_PARAMS)
        })
        mice_map_params.append({
            "r":         eff[0],
            "K":         math.exp(eff[1]),
            "k_kill":    math.exp(eff[2]),
            "k_clear":   math.exp(eff[3]),
            "k_cyto":    math.exp(eff[4]),
            "beta_dose": math.exp(eff[5]),
            "mu_base":   math.exp(eff[6]),
            "mu_conf":   math.exp(eff[7]),
            "sigma_B":   math.exp(p[N_BIO_PARAMS]),
            "pk":        pk_map,
            "energy":    best_energy,
        })

    weights = np.ones(len(post_samples), dtype=float) / max(len(post_samples), 1)

    if verbose:
        print(f"\n{'='*60}\nJoint MCMC Results\n{'='*60}")
        print(f"  Global MAP : r={global_map['r']:.4f}  K={global_map['K']:.3e}  "
              f"k_kill={global_map['k_kill']:.3e}  β={global_map['beta_dose']:.4f}")
        print(f"               mu_base={global_map['mu_base']:.4f}  "
              f"mu_conf={global_map['mu_conf']:.4f}  "
              f"k_cyto={global_map['k_cyto']:.3e}  k_clear={global_map['k_clear']:.4f}")
        print(f"  Energy     : {best_energy:.4f}   Posterior samples: {len(post_samples)}")
        print(f"  Per-mouse MAP epsilons (natural-space for r, log-space for the rest):")
        for m, (name, eps) in enumerate(zip(sample_names, mice_eps_map)):
            eps_vals = "  ".join(f"{k}={v:+.3f}" for k, v in eps.items())
            print(f"    [{m:2d}] {name}: {eps_vals}")

    pool.close()
    pool.join()

    return {
        "global_map":      global_map,
        "mice_eps_map":    mice_eps_map,
        "mice_map_params": mice_map_params,
        "post_samples":    post_samples,
        "weights":         weights,
        "all_trace":       all_trace,
        "sample_names":    sample_names,
    }

# ─────────────────────────────────────────────────────────────────────────────
# Model selection: AIC / BIC against the (ploidy-derived) tumour-burden data
# ─────────────────────────────────────────────────────────────────────────────
#
# The MCMC objective (_log_posterior_single / _log_posterior_joint) returns a
# negative *log-posterior* — i.e. the burden log-likelihood PLUS the PK priors
# (and, in the joint model, the per-mouse epsilon prior).  AIC and BIC are
# defined purely in terms of the log-*likelihood* of the data, so the stored
# `energy` cannot be reused directly.  The helpers below recompute the pure
# burden log-likelihood at a fixed (MAP) parameter point, using the exact same
# log-normal noise model as the objective, then form
#
#       AIC = 2k − 2·logL
#       BIC = k·ln(n) − 2·logL
#
# where n is the number of burden datapoints the model is compared against and
# k is the number of free parameters.  The burden observations here are the
# tumour burdens derived from the ploidy-count timeline (see
# load_treatment_day_ploidy_from_oxygen / the harvest burden loader in
# beam_search_flip_rate_wgd.py), so this is the AIC/BIC of the model *relative
# to the ploidy burden datapoints*.


def _burden_loglik_at_map(caller, initial_ploidy, first_tx_day,
                          obs, obs_drugs, map_params, pk_state, sigma_B):
    """
    Pure burden log-likelihood (no priors) at a fixed parameter point.

    Uses the SAME per-observation log-normal noise model as the MCMC
    objective:

        log L = Σ_d  −[ log σ_B + 0.5·log(2π)
                        + (log v_obs − log v_pred)² / (2 σ_B²) ]

    The caller module globals (INITIAL_PLOIDY, FIRST_TX_DAY, …) are set so the
    simulation reproduces exactly what the sampler saw for this mouse.

    Returns
    -------
    (loglik, n_points)
        loglik   : float  — burden log-likelihood (nan if the simulation is
                            infeasible at these params, e.g. tumour goes
                            extinct / non-finite).
        n_points : int    — number of burden observations (days > 0) compared.
    """
    caller.INITIAL_PLOIDY = dict(initial_ploidy)
    if first_tx_day is not None:
        caller.FIRST_TX_DAY = float(first_tx_day)
    caller.OBSERVED_TUMOR_BURDENS      = dict(obs)
    caller.OBSERVED_DRUGS_ADMINISTERED = list(obs_drugs)

    obs_days = sorted(d for d in obs if d > 0)
    if len(obs) <= 1 or not obs_days:
        return 0.0, 0

    timeline, _ = caller._simulate_burden_and_ploidy(
        obs_drugs,
        map_params["r"], map_params["K"], map_params["k_kill"],
        map_params["k_clear"], map_params["k_cyto"], map_params["beta_dose"],
        map_params["mu_base"], map_params["mu_conf"],
        pk_state=pk_state,
    )
    if timeline is None:
        return float("nan"), len(obs_days)

    pred = {round(t, 2): v for t, v in timeline}
    loglik = 0.0
    n = 0
    for d in obs_days:
        v_obs     = obs[d]
        closest_t = min(pred.keys(), key=lambda t: abs(t - d))
        v_pred    = pred[closest_t]
        if v_pred <= 0 or v_obs <= 0:
            return float("nan"), len(obs_days)
        residual = (math.log(v_obs) - math.log(v_pred)) ** 2
        loglik  += -(math.log(sigma_B)
                     + 0.5 * math.log(2 * math.pi)
                     + residual / (2 * sigma_B ** 2))
        n += 1
    return loglik, n


def _aic_bic(loglik, k, n):
    """Standard AIC / BIC from a log-likelihood, #params k and #datapoints n."""
    aic = 2 * k - 2 * loglik
    bic = (k * math.log(n) - 2 * loglik) if n > 0 else float("nan")
    return aic, bic


def _count_fitted_pk_params(pk_params_to_fit, obs_drug_set):
    """
    Number of PK parameters the sampler actually fits — i.e. PK params whose
    drug appears in the observed schedule.  Mirrors the filtering in run_mcmc /
    run_mcmc_joint so the parameter count k is exact.
    """
    n_pk = 0
    for drug, params in pk_params_to_fit.items():
        if drug.lower() not in obs_drug_set:
            continue
        n_pk += len(params)
    return n_pk


def compute_aic_bic_single(caller, map_params, initial_ploidy, obs, obs_drugs,
                           pk_params_to_fit, first_tx_day=None):
    """
    AIC / BIC for a single-mouse fit against its ploidy-burden datapoints.

    Free parameters (k):
        8 biological params (r, K, k_kill, k_clear, k_cyto, beta_dose,
        mu_base, mu_conf) + sigma_B + fitted PK params.
        (k_tr is interpolated from ploidy, not fitted, so it is NOT counted.)

    Returns a dict with loglik, n_datapoints, k_params, aic, bic.
    """
    sigma_B  = map_params["sigma_B"]
    pk_state = map_params.get("pk")

    loglik, n = _burden_loglik_at_map(
        caller, initial_ploidy, first_tx_day, obs, obs_drugs,
        map_params, pk_state, sigma_B,
    )

    obs_drug_set = {drug.lower() for _, _, drug, _ in obs_drugs}
    n_pk = _count_fitted_pk_params(pk_params_to_fit, obs_drug_set)
    k = N_BIO_PARAMS + 1 + n_pk        # 8 bio + sigma_B + PK

    aic, bic = _aic_bic(loglik, k, n)
    return {
        "loglik":        loglik,
        "n_datapoints":  n,
        "k_params":      k,
        "n_bio_params":  N_BIO_PARAMS,
        "n_pk_params":   n_pk,
        "aic":           aic,
        "bic":           bic,
    }


def compute_aic_bic_joint(caller, global_map, mice_map_params, mice_data,
                          pk_params_to_fit, sample_names=None):
    """
    AIC / BIC for the joint hierarchical fit against ALL mice's ploidy-burden
    datapoints pooled together.

    Free parameters (k):
        8 global biological params + sigma_B
        + n_mice × 8 per-mouse epsilon offsets
        + fitted PK params (shared across mice).

    The log-likelihood is the sum of every mouse's burden log-likelihood,
    each evaluated with that mouse's EFFECTIVE MAP params (global + epsilon)
    and the shared global sigma_B / PK state.

    Returns a dict with the pooled totals plus a per-mouse breakdown.
    """
    sigma_B  = global_map["sigma_B"]
    pk_state = global_map.get("pk")

    total_loglik = 0.0
    total_n      = 0
    per_mouse    = []
    any_infeasible = False

    for m, (mouse, eff) in enumerate(zip(mice_data, mice_map_params)):
        ll, n = _burden_loglik_at_map(
            caller,
            mouse["initial_ploidy"], mouse.get("first_tx_day"),
            mouse["obs"], mouse["drugs"],
            eff, pk_state, sigma_B,
        )
        name = (sample_names[m] if sample_names and m < len(sample_names)
                else f"mouse_{m}")
        per_mouse.append({
            "harvest":       name,
            "loglik":        ll,
            "n_datapoints":  n,
        })
        if math.isnan(ll):
            any_infeasible = True
        else:
            total_loglik += ll
        total_n += n

    n_mice = len(mice_data)

    all_drug_set = set()
    for mouse in mice_data:
        all_drug_set.update(drug.lower() for _, _, drug, _ in mouse["drugs"])
    n_pk = _count_fitted_pk_params(pk_params_to_fit, all_drug_set)

    # k = 8 global bio + sigma_B + n_mice*8 epsilons + PK
    k = (N_BIO_PARAMS + 1) + n_mice * N_BIO_PARAMS + n_pk

    loglik = float("nan") if any_infeasible else total_loglik
    aic, bic = _aic_bic(loglik, k, total_n)

    return {
        "loglik":         loglik,
        "n_datapoints":   total_n,
        "n_mice":         n_mice,
        "k_params":       k,
        "n_global_bio":   N_BIO_PARAMS,
        "n_epsilon":      n_mice * N_BIO_PARAMS,
        "n_pk_params":    n_pk,
        "aic":            aic,
        "bic":            bic,
        "per_mouse":      per_mouse,
    }
