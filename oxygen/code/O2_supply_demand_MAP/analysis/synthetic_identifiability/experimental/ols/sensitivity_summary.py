"""Compare checked real-data fits and coarse profile probes under changed assumptions."""
from inference import THRESHOLD, likelihood_ratio
from practical import HERE, MECHANISMS, read, write


def summarize():
    baseline = read(HERE / "main_constrained_best.tsv")[0]
    report = []
    for case in ("noise", "passage"):
        fit = read(HERE / f"{case}_constrained_best.tsv")[0]
        primary_loss = float(fit["loss"])
        profiles = read(HERE / f"{case}_profiles.tsv")
        retuned = HERE / f"{case}_profiles_retuned.tsv"
        if retuned.exists():
            profiles += read(retuned)
        valid = [r for r in profiles if r["feasible"] == "True"]
        scores = [(likelihood_ratio(float(r["loss"]), primary_loss, 114), r) for r in valid]
        reference_improvement = min(lr for lr, _ in scores)
        for p in MECHANISMS:
            x = float(fit[f"x_{p}"])
            shift = x - float(baseline[f"x_{p}"])
            outside = [(lr, r) for lr, r in scores if abs(float(r[f"x_{p}"]) - x) > .10000001]
            tested = {r["fixed_normalized"] for _, r in outside if r["profile_parameter"] == p}
            projected = {r[f"x_{p}"] for lr, r in outside if lr <= THRESHOLD and r["profile_parameter"] != p}
            levels = [lr for lr,r in outside if r["profile_parameter"] == p or lr <= THRESHOLD]
            lowest = min(levels, default=float("inf"))
            status = "reference_unresolved" if reference_improvement < -.10 else (
                "outside_precision_supported" if lowest <= THRESHOLD
                else "no_outside_support_at_tested_points")
            report.append(dict(case=case, parameter=p, baseline_estimate=baseline[p], variant_estimate=fit[p],
                               normalized_point_estimate_shift=shift, shift_larger_than_tolerance=abs(shift)>.1,
                               tested_points_beyond_tolerance=len(tested),
                               accepted_fit_projections_beyond_tolerance=len(projected),
                               minimum_lr_beyond_tolerance=lowest, outside_precision_supported=lowest<=THRESHOLD,
                               reference_KKT_norm=fit["max_projected_gradient"],
                               best_profile_lr_against_reference=reference_improvement, sensitivity_status=status))
    write(HERE / "sensitivity_summary.tsv", report)


if __name__ == "__main__":
    summarize()
