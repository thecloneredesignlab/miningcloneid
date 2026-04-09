from __future__ import annotations

import math
import os
from multiprocessing import Pool
from time import time

import numpy as np
import torch
import importlib

N_CHAINS       = 32
N_BURNIN       = 300
N_SAMPLES      = 700
THIN           = 2
ADAPT_INTERVAL = 50
TARGET_ACCEPT  = 0.234 # https://www.sciencedirect.com/science/article/pii/S0304414919306982#sec3

PARAM_BOUNDS = {
    "log_r_base":  (math.log(0.01),  math.log(2.0)),
    "log_k_cap":   (math.log(1e8),   math.log(1e12)),
    "log_beta":    (math.log(1e-4),  math.log(1.0)),
    "log_sigma_B": (math.log(0.01),  math.log(5.0)),
    "log_sigma_C": (math.log(0.1),   math.log(50.0)),
}

INIT_VALS = {
    "log_r_base":  math.log(0.28),
    "log_k_cap":   math.log(6e10),
    "log_beta":    math.log(0.05),
    "log_sigma_B": math.log(0.5),
    "log_sigma_C": math.log(2.0),
}

INIT_STEP = {
    "log_r_base":  0.05,
    "log_k_cap":   0.10,
    "log_beta":    0.10,
    "log_sigma_B": 0.05,
    "log_sigma_C": 0.05,
}

_W = {}


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

    _W["fn_burden"]  = caller._simulate_burden_timeline
    _W["fn_ploidy"]  = caller.get_observed_end_ploidy

# Objective function; lower is better
def _log_posterior_single(param_vec, pk_param_names, pk_param_drugs):
    log_r, log_k, log_b, log_sB, log_sC = param_vec[:5]
    r_base  = math.exp(log_r)
    k_cap   = math.exp(log_k)
    beta    = math.exp(log_b)
    sigma_B = math.exp(log_sB)
    sigma_C = math.exp(log_sC)

    pk_state = {}
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_state.setdefault(drug, {})[pname] = math.exp(param_vec[5 + i])

    obs        = _W["obs"]
    obs_drugs  = _W["drugs"]
    end_pl     = _W["end_ploidy"]
    pk_cfg     = _W["pk_cfg"]

    # L_B = (1/n) * Σ_(i,n) (log(sigma_B) + 0.5log(2pi) + (log|v_observed| - log |v_hat|) ** 2 / (2 * sigma_B) )
    nll = 0.0
    if len(obs) > 1:
        timeline = _W["fn_burden"](obs_drugs, r_base, k_cap,
                                   pk_state=pk_state, beta=beta)
        if timeline is None:
            return 1e12
        pred = {round(t, 2): v for t, v in timeline}
        obs_days = sorted(d for d in obs if d > 0)
        if obs_days:
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
            nll += burden_nll / len(obs_days)

    # L_C
    if end_pl.size > 0:
        pred_ploidy = _W["fn_ploidy"](r_base, k_cap, pk_state, beta=beta)
        if pred_ploidy is None:
            return 1e12
        total = sum(pred_ploidy.values())
        if total <= 0:
            return 1e12
        chr_arr = np.array(list(pred_ploidy.keys()), dtype=float)
        w_arr   = np.array(list(pred_ploidy.values()), dtype=float)
        w_arr  /= w_arr.sum()

        ploidy_nll = 0.0
        for z in end_pl:
            log_mix = -np.inf
            for N_val, w in zip(chr_arr, w_arr):
                if w <= 0: # population size < 0
                    continue
                lc = (math.log(w)
                      - 0.5 * math.log(2 * math.pi * sigma_C ** 2)
                      - (z - N_val) ** 2 / (2 * sigma_C ** 2))
                if log_mix == -np.inf:
                    log_mix = lc
                else:
                    log_mix = max(log_mix, lc) + math.log(
                        1 + math.exp(-abs(log_mix - lc)))
            ploidy_nll -= log_mix
        nll += ploidy_nll / len(end_pl)

    # PK priors
    for i, (pname, drug) in enumerate(zip(pk_param_names, pk_param_drugs)):
        cfg = pk_cfg[drug][pname]
        nll += 0.5 * ((param_vec[5 + i] - cfg["prior_log_mean"]) / cfg["prior_log_std"]) ** 2

    return nll


def _eval_worker(args):
    ci, pv, names, drugs = args
    return ci, _log_posterior_single(pv, names, drugs)


def run_mcmc(
    initial_ploidy, observed_burdens, observed_end_ploidy, observed_drugs,
    pk_params_to_fit, haploid_n, sample_name,
    fn_get_end_ploidy=None, fn_simulate_burden=None, fn_fill_gaps=None,
    caller_module_name="test", verbose=True,
):
    base_names = ["log_r_base", "log_k_cap", "log_beta", "log_sigma_B", "log_sigma_C"]
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
        print(f"  PK drugs fitted: {sorted(set(pk_param_drugs))}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if verbose:
        print(f"  Device: {device}")
    
    init_vec = np.array([INIT_VALS[n] for n in base_names] + pk_init_vals, dtype=np.float64)
    step_vec = np.array([INIT_STEP[n] for n in base_names] + pk_init_steps, dtype=np.float64)
    lo_vec = np.array([PARAM_BOUNDS[n][0] for n in base_names] + [b[0] for b in pk_bounds], dtype=np.float64)
    hi_vec = np.array([PARAM_BOUNDS[n][1] for n in base_names] + [b[1] for b in pk_bounds], dtype=np.float64)

    rng = np.random.default_rng(42)
    chains_np = np.tile(init_vec, (N_CHAINS, 1)) + rng.normal(0, step_vec * 0.1, size=(N_CHAINS, n_dim))
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
    n_workers = min(N_CHAINS, os.cpu_count() or 4)

    def _eval_all(chain_t):
        c_np = chain_t.cpu().numpy()
        args = [(ci, c_np[ci], pk_param_names, pk_param_drugs) for ci in range(N_CHAINS)]
        results = [None] * N_CHAINS
        with Pool(n_workers, initializer=_init_worker,
                  initargs=(caller_module_name, worker_state)) as pool:
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

        proposals = chains + torch.randn_like(chains) * steps.unsqueeze(0) # random step
        proposals = torch.clamp(proposals, lo_t, hi_t)
        prop_energies = _eval_all(proposals) # objective function run; lower is better

        log_alpha = energies - prop_energies # difference in objective; objective( X_(n-1) ) - objective( X_(n) ); positive log_alpha means better new parameter
        accept_mask = torch.rand(N_CHAINS, device=device, dtype=torch.float64).log() < log_alpha # log( (0,1) ) -> (-inf, 0) ; if difference in objective(theta) is greater than randomness, accept new proposal (flipped d.t. objective being negative)
        chains   = torch.where(accept_mask.unsqueeze(1), proposals, chains)
        energies = torch.where(accept_mask, prop_energies, energies)
        accept_count += accept_mask.float()

        ci_best = energies.argmin().item() # keep track of best parameters
        if energies[ci_best].item() < best_energy:
            best_energy = energies[ci_best].item()
            best_params = chains[ci_best].cpu().numpy().copy()

        if is_burnin and (it + 1) % ADAPT_INTERVAL == 0:
            mr = (accept_count / (it + 1)).mean().item()
            if mr < TARGET_ACCEPT - 0.05:
                steps *= 0.8 # if lots of acceptances, grow step size and explore
            elif mr > TARGET_ACCEPT + 0.05:
                steps *= 1.2 # if too few acceptances, decrease step size and exploit
            if verbose and (it + 1) % 100 == 0:
                print(f"  {it+1:4d} burn-in | acc={mr:.3f} | best={best_energy:.2f} | {time()-it_t0:.1f}s")

        # Animation
        p0 = chains[0].cpu().numpy()
        pk0 = {}
        for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
            pk0.setdefault(dr, {})[pn] = math.exp(p0[n_base + i])
        all_trace.append({"iter": it, "burnin": is_burnin,
                          "r_base": math.exp(p0[0]), "k_cap": math.exp(p0[1]),
                          "beta": math.exp(p0[2]), "pk": pk0,
                          "logP": -energies[0].item()})

        # After tuning, start sampling
        if not is_burnin and (it - N_BURNIN) % THIN == 0:
            c_all = chains.cpu().numpy()
            for ci in range(N_CHAINS):
                pv = c_all[ci]
                pkd = {}
                for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
                    pkd.setdefault(dr, {})[pn] = math.exp(pv[n_base + i])
                post_samples.append({
                    "r_base": math.exp(pv[0]), "k_cap": math.exp(pv[1]),
                    "beta": math.exp(pv[2]), "sigma_B": math.exp(pv[3]),
                    "sigma_C": math.exp(pv[4]), "pk": pkd,
                    "energy": energies[ci].item(),
                })

        if verbose and not is_burnin and (it + 1) % 100 == 0:
            mr = (accept_count / (it + 1)).mean().item()
            print(f"  {it+1:4d} sample | acc={mr:.3f} | best={best_energy:.2f} | "
                  f"n={len(post_samples)} | {time()-it_t0:.1f}s")

    # MAP <- best parameter
    pk_map = {}
    for i, (pn, dr) in enumerate(zip(pk_param_names, pk_param_drugs)):
        pk_map.setdefault(dr, {})[pn] = math.exp(best_params[n_base + i])
    for drug, params in pk_params_to_fit.items():
        if drug not in pk_map:
            pk_map[drug] = {p: c["init"] for p, c in params.items()}

    map_params = {
        "r_base": math.exp(best_params[0]), "k_cap": math.exp(best_params[1]),
        "beta": math.exp(best_params[2]), "sigma_B": math.exp(best_params[3]),
        "sigma_C": math.exp(best_params[4]), "pk": pk_map, "energy": best_energy,
    }
    weights = np.ones(len(post_samples), dtype=float) / max(len(post_samples), 1)

    if verbose:
        print(f"\n{'='*60}\nMCMC Results\n{'='*60}")
        print(f"  MAP: r={map_params['r_base']:.5f} K={map_params['k_cap']:.3e} "
              f"β={map_params['beta']:.5f} σ_B={map_params['sigma_B']:.4f} σ_C={map_params['sigma_C']:.4f}")
        print(f"  Energy: {best_energy:.4f}  Samples: {len(post_samples)}")
        if post_samples:
            for k in ("r_base", "k_cap", "beta"):
                vals = [s[k] for s in post_samples]
                print(f"  E[{k}] = {np.mean(vals):.6g} ± {np.std(vals):.4g}")

    return map_params, post_samples, weights, all_trace
