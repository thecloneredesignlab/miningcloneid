"""
regenerate_observed_fits.py
────────────────────────────
Backfill observed_fit_tumor_burden.png for runs that were completed before
that plot existed, WITHOUT re-running the MCMC fit.

It reuses each mouse's already-fitted parameters — the per-mouse
"effective_params" (global MAP × that mouse's fitted epsilons) and the shared
"pk_map" — stored in results_joint/joint_global_params.json, reloads each
mouse's observed drug schedule and burden data, and calls
beam_search_flip_rate_wgd.plot_observed_fit() to render the fitted model on
the drugs the mouse actually received, overlaid on the data points.

The output is written next to the other per-mouse plots, i.e.
    results_joint/<harvest_name>/observed_fit_tumor_burden.png
so that compile_fit_grid.py picks it up automatically.

In addition (unless --no-uniform is passed) it renders a SECOND plot per
mouse,
    results_joint/<harvest_name>/observed_fit_tumor_burden_uniform.png
which is identical to the first except that every mouse's y-axis is fixed to
the SAME (shared) bounds. Those shared bounds are computed in a first pass as
the min/max of every mouse's predicted curve and observed data points, padded
a little on the log scale. The per-mouse plot the first pass produces still
autoscales individually; only the *_uniform.png copies share bounds. Tiling
those copies (compile_fit_grid.py --uniform) yields a grid whose panels are
directly comparable in height.

Usage:
    python regenerate_observed_fits.py
    python regenerate_observed_fits.py --params results_joint/joint_global_params.json
    python regenerate_observed_fits.py --results-dir results_joint
    python regenerate_observed_fits.py --no-uniform          # skip shared-bounds copies
    python regenerate_observed_fits.py --uniform-pad 2.0     # more headroom on shared bounds
"""

import argparse
import json
import os

import beam_search_flip_rate_wgd as bs

UNIFORM_FILENAME = "observed_fit_tumor_burden_uniform.png"


def _simulated_and_observed_values(md, eff, pk_map):
    """
    Positive burden values (predicted curve + observed points) for one mouse.

    Used only to determine the shared y-axis range; mirrors exactly what
    plot_observed_fit() draws. Non-positive values are dropped because the
    plot uses a log y-axis. Requires the module globals used by
    _simulate_burden_timeline (INITIAL_PLOIDY, FIRST_TX_DAY) to already be set
    for this mouse.
    """
    vals = []
    tl = bs._simulate_burden_timeline(
        md["drugs"],
        eff["r"], eff["K"], eff["k_kill"], eff["k_clear"],
        eff["k_cyto"], eff["beta_dose"], eff["mu_base"], eff["mu_conf"],
        pk_state=pk_map,
    )
    if tl is not None:
        vals.extend(b for _, b in tl)
    vals.extend(md["obs"].values())
    return [v for v in vals if v is not None and v > 0]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=os.path.join("results_joint",
                                                     "joint_global_params.json"),
                    help="Path to joint_global_params.json (holds per-mouse "
                         "effective_params and the shared pk_map)")
    ap.add_argument("--results-dir", default="results_joint",
                    help="Directory containing per-harvest result subfolders")
    ap.add_argument("--no-uniform", dest="uniform", action="store_false", default=True,
                    help="Only render the individually-autoscaled "
                         "observed_fit_tumor_burden.png; skip the shared-y-bounds "
                         f"{UNIFORM_FILENAME} copies")
    ap.add_argument("--uniform-filename", default=UNIFORM_FILENAME,
                    help=f"Filename for the shared-y-bounds copy (default: {UNIFORM_FILENAME})")
    ap.add_argument("--uniform-pad", type=float, default=1.5,
                    help="Multiplicative padding applied to the shared y-range "
                         "(log scale): ymin is divided by and ymax multiplied by "
                         "this factor. Default 1.5; use 1.0 for a tight fit")
    args = ap.parse_args()

    with open(args.params) as fh:
        gp = json.load(fh)

    pk_map = gp.get("pk_map", {})
    mice_params = gp.get("mice", [])
    if not mice_params:
        raise RuntimeError(f"No per-mouse entries found in {args.params}")

    print(f"Loaded fitted parameters for {len(mice_params)} mice from {args.params}")

    def _set_globals_for(md):
        """Point the module globals _simulate_burden_timeline relies on at this mouse."""
        bs.INITIAL_PLOIDY              = md["initial_ploidy"]
        bs.FIRST_TX_DAY                = md["first_tx_day"]
        bs.OBSERVED_DRUGS_ADMINISTERED = md["drugs"]
        bs.OBSERVED_TUMOR_BURDENS      = md["obs"]

    # ── Pass 1: render each mouse's individually-autoscaled plot ─────────────
    # Cache (harvest, results_dir, eff, md) for the ready mice so the optional
    # uniform pass doesn't have to reload / re-validate them.
    n_ok = n_fail = 0
    ready: list[tuple[str, str, dict, dict]] = []
    global_min = None
    global_max = None

    for entry in mice_params:
        harvest = entry["harvest"]
        eff = entry.get("effective_params")
        if eff is None:
            print(f"  SKIP {harvest}: no effective_params in JSON")
            n_fail += 1
            continue

        results_dir = os.path.join(args.results_dir, harvest)
        if not os.path.isdir(results_dir):
            print(f"  SKIP {harvest}: results dir {results_dir} not found")
            n_fail += 1
            continue

        # Reload this mouse's observed data + schedule, and set the module
        # globals that _simulate_burden_timeline() relies on.
        loaded, _valid = bs.load_all_harvest_data_for_joint([harvest])
        if not loaded:
            print(f"  SKIP {harvest}: could not reload observed data")
            n_fail += 1
            continue
        md = loaded[0]

        _set_globals_for(md)

        out = bs.plot_observed_fit(
            results_dir, harvest,
            map_params=eff,
            pk_state=pk_map,
            observed_schedule=md["drugs"],
            observed_burdens=md["obs"],
        )
        if out is not None:
            n_ok += 1
            ready.append((harvest, results_dir, eff, md))
            if args.uniform:
                vals = _simulated_and_observed_values(md, eff, pk_map)
                if vals:
                    lo, hi = min(vals), max(vals)
                    global_min = lo if global_min is None else min(global_min, lo)
                    global_max = hi if global_max is None else max(global_max, hi)
        else:
            n_fail += 1

    print(f"\nAutoscaled pass done: {n_ok} regenerated, {n_fail} skipped/failed.")

    # ── Pass 2 (optional): re-render every ready mouse on shared y-bounds ────
    if args.uniform:
        if global_min is None or global_max is None or global_min >= global_max:
            print("Uniform pass skipped: could not determine a shared y-range "
                  "(no positive burden values across mice).")
        else:
            pad = max(args.uniform_pad, 1.0)
            shared_ylim = (global_min / pad, global_max * pad)
            print(f"\nUniform pass: shared y-bounds = "
                  f"[{shared_ylim[0]:.3e}, {shared_ylim[1]:.3e}] "
                  f"(pad x{pad}) across {len(ready)} mice")

            n_uni = 0
            for harvest, results_dir, eff, md in ready:
                _set_globals_for(md)
                out = bs.plot_observed_fit(
                    results_dir, harvest,
                    map_params=eff,
                    pk_state=pk_map,
                    observed_schedule=md["drugs"],
                    observed_burdens=md["obs"],
                    filename=args.uniform_filename,
                    ylim=shared_ylim,
                )
                if out is not None:
                    n_uni += 1
            print(f"Uniform pass done: {n_uni} shared-bounds plots "
                  f"written as {args.uniform_filename}.")


if __name__ == "__main__":
    main()