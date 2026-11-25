"""Collect the 14 J curves and project calibration cost from the completed pilot."""
from collections import Counter, defaultdict
from math import ceil
from pathlib import Path
import sys

from scan import PILOT, Run, independent_score, read, write
import numpy as np


def collect(out):
    profile = read(out/"profile_grid.tsv")
    assert len({(r["parameter"],r["value"]) for r in profile})==len(profile)
    support = {float(r["value"]):r for r in read(out/"lam_max_supply_envelope.tsv")
               if r["violates_live_supply_envelope"]=="True"}
    for r in profile:
        envelope_failed = r["parameter"]=="lam_max" and float(r["value"]) in support
        r["feasibility_evidence"] = "Live supply envelope rules out point" if envelope_failed else ""
    write(out/"profile_grid.tsv",profile)
    groups = defaultdict(list)
    for row in profile:
        groups[row["parameter"]].append(row)
    assert len(groups)==14
    base,=read(out/"baseline.tsv")
    assert base["status"]=="verified"
    assert all(r["baseline_loss"]==base["selected_loss"] for r in profile)
    initial = read(out/"baseline_vector.tsv")[0]
    counts = Counter(row["status"] for row in profile)
    audit = Run(out)
    model = audit.model
    p = audit.conf.unmap(audit.baseline_x())[0]
    prediction = independent_score(model,p)
    branches = []
    for index in np.flatnonzero(model.death_branch):
        slot = model.classifier_slots[index]
        count = model.death_log_index[index]
        assert slot==model.indices[count] and model.count_branch[count]
        observed = model.dead[index]/model.eligible[index]
        fitted = prediction["fraction"][index]
        branches.append(dict(event=model.slot_ids[slot],observed_count=np.exp(model.y[count]),
            fitted_visible_count=np.exp(prediction["predicted_logs"][count]),death_denominator=model.eligible[index],
            observed_death_fraction=observed,fitted_death_fraction=fitted,death_error_percentage_points=100*(fitted-observed)))
    assert len(branches)==4
    write(out/"a7k_baseline_residuals.tsv",branches)
    a7k_errors=[r["death_error_percentage_points"] for r in branches]
    summary=[]
    for name,rows in groups.items():
        states=Counter(r["status"] for r in rows)
        ruled_out=sum(bool(r["feasibility_evidence"]) for r in rows)
        check=[r for r in rows if r["status"]=="verified"]
        ref=[r for r in rows if r["rise"]]
        assert len(rows)>=10 and check
        def at_boundary(x):
            point=min(rows,key=lambda r:abs(float(r["normalized"])-x))
            assert abs(float(point["normalized"])-x)<1e-10
            return point["status"],point["rise"]
        low,high=at_boundary(0),at_boundary(1)
        summary.append(dict(parameter=name,estimate=initial[name],
                estimate_normalized=initial["x_"+name],
                estimated_box_face=next((s.split(":")[1] for s in base["mechanistic_box_faces"].split()
                                         if s.split(":")[0]==name),""),
                grid_points=len(rows),verified=states["verified"],provisional=states["provisional"],
                supply_infeasible=ruled_out,feasibility_unresolved=states["feasibility_unresolved"]-ruled_out,
                lower_status=low[0],lower_rise=low[1],
                upper_status=high[0],upper_rise=high[1],
                lowest_grid_rise=min(float(r["rise"]) for r in ref),
                highest_verified_rise=max(float(r["rise"]) for r in check)))
    write(out/"mechanism_summary.tsv",summary)
    zeros=[r for r in profile if r["parameter"] in PILOT[:2] and r["exact_zero"]=="True"]
    assert len(zeros)==2
    write(out/"exact_zero_summary.tsv",zeros)
    paired = [r for n in PILOT for r in read(out/"bootstrap_pilot"/n/"pairs.tsv")]
    assert len(paired)==12 and all(r["conditional_death_replay"]=="True" for r in paired)
    assert max(float(r["likelihood_scale_agreement"]) for r in paired)<2e-7
    write(out/"bootstrap_pilot_pairs.tsv",paired)
    valid=sum(r["pair_status"]=="verified" for r in paired)
    earlier = [r for n in PILOT for r in read(out/"bootstrap_pilot"/n/"pairs_initial.tsv")]
    preliminary=sum(r["pair_status"]=="verified" for r in earlier) if len(earlier)==12 else None
    time=np.array([float(r["total_cpu_seconds"]) for r in paired])
    p50,p75,p90=(float(np.quantile(time,q)) for q in (.50,.75,.90))
    planning = 1.25*p75
    scenarios = []
    for label,anchors,n in (("Four pilot anchors",4,199),("Thirty candidate values",30,999),
                            ("Sixty candidate values",60,999),
                            ("Extra per difficult candidate",1,1000)):
        scenarios.append(dict(scope=label,candidate_values=anchors,replicates_each=n,paired_refits=anchors*n,
                 empirical_median_core_hours=anchors*n*p50/3600.,
                 empirical_p90_core_hours=anchors*n*p90/3600.,
                 planning_cpu_seconds_per_pair=planning,planning_core_hours=ceil(anchors*n*planning/3600.),
                 status="Requires approval"))
    write(out/"calibration_budget.tsv",scenarios)
    observed=[r for name in groups for r in read(out/"profiles"/name/"attempts.tsv")]
    baseline=read(out/"baseline"/"attempts.tsv")
    def phase(name,rows):
        return dict(phase=name,optimizer_attempts=sum(r["origin"]=="fit" for r in rows),
                    replayed_attempts=sum(r["origin"]=="archived Step 3" for r in rows),
                    cpu_seconds=sum(float(r["cpu_seconds"]) for r in rows))
    runtime=[phase("Observed baseline",baseline),phase("Observed profile grid and pilot",observed),
             dict(phase="Twelve bootstrap pilot pairs",optimizer_attempts=sum(int(r["optimization_attempts"]) for r in paired),
                  replayed_attempts=0,cpu_seconds=float(np.sum(time)))]
    write(out/"runtime_summary.tsv",runtime)
    obs_optimizers=sum(r["optimizer_attempts"] for r in runtime[:2])
    fidelity=max(float(r["scaled_error"]) for r in read(out/"gradient_checks.tsv"))
    rows=["| Mechanism | J estimate | Fitted face | Verified | Provisional | Supply infeasible | Feasibility unresolved |",
          "| --- | ---: | --- | ---: | ---: | ---: | ---: |"]
    for r in summary:
        rows.append(f"| `{r['parameter']}` | {float(r['estimate']):.7g} | {r['estimated_box_face'] or 'Interior'} | "
                    f"{r['verified']} | {r['provisional']} | {r['supply_infeasible']} | {r['feasibility_unresolved']} |")
    zero_lines=[]
    for r in zeros:
        change=f"{float(r['rise']):.6f}" if r["rise"] else "Unavailable"
        zero_lines.append(f"- `{r['parameter']} = 0` had objective rise {change} and status {r['status']}. "
                          f"Inactive carrier `{r['inactive_at_endpoint'] or ('k_o_mis' if r['parameter']=='p_misseg' else 'buffer_n_exp')}`.")
    references=[]
    for p in PILOT:
        v=[r for r in paired if r["parameter"]==p]
        hits=sum(r["pair_status"]=="verified" for r in v)
        seconds=np.array([float(r["total_cpu_seconds"]) for r in v])
        references.append(f"| `{p}` | {float(v[0]['value']):.7g} | {hits}/3 | {np.median(seconds):.1f} | "
                          f"{np.min(seconds):.1f}–{np.max(seconds):.1f} |")
    hypothesis="The one reserve pass used the same pseudo datasets." if preliminary is not None else ""
    delta=f"Before this reserve pass, {preliminary}/12 pairs verified. " if preliminary is not None else ""
    report=f"""# Three-readout profiles for 14 in vitro mechanism parameters

## Observed-data result

The strict donor-inheritance J reference is **{float(base['selected_loss']):.9f}** with
{base['agreeing_starts']} independently initialized qualifying starts. There are
**{len(profile)} fixed-value points across 14 mechanisms**. Of these, **{counts['verified']}**
are verified, **{counts['provisional']}** retain a replayed feasible candidate,
**{len(support)}** violate the donor-supply envelope, and
**{counts['feasibility_unresolved']-len(support)}** have unresolved feasibility.
Each nominal box endpoint was attempted. The eight values ruled out by the
envelope belong to `lam_max`.

For the supply check, each living class can contribute at most two daughters per
division, the division rate is bounded by `lam_max`, and survival and crowding
weights are at most one. The upper population multiplier per Euler step is
`1 + 0.05 * lam_max`. This bound also applies when negative updates are
clipped because `0.05 * lam_max <= 1`. Multiplying it across the
donor's observed passage and comparing against the grouped withdrawal demand
provides a necessary feasibility check. The most restrictive checked low-rate
donor requires about 0.46225/day even under this growth ceiling.

[Donor supply envelope](lam_max_supply_envelope.tsv)

The J objective uses corrected visible counts, measured death classification and
chromosome distributions for the verified events. Donor propagation, 118 conditioned
seed totals, two fixed ancestry proxies and the 110 supply constraints are from the
V4 design. The visible count is distinct from the latent living component. The 20
flow records and two A19 histograms remain held under their specimen gates.

The scorer binds 225 counts, 195 classification events and six 133-bin chromosome
samples totaling 104 cells. The count residual used in each death term maps to the
same event. The maximum checked relative-scaled gradient discrepancy was {fidelity:.3g}.
Each passing profile row has two distinct numerically qualified starts near the
lowest replayed feasible score. The archived stationarity and score agreement gates
were applied and each qualified endpoint received a feasible-direction audit.
A profile row remains conditional on this resolution and finite set of starts.
{obs_optimizers} new observed-data optimizer attempts were recorded.

{chr(10).join(rows)}

## Exact mechanism-null fits

{chr(10).join(zero_lines)}

These are rises on the fitted objective scale. The interval threshold is pending
calibration. The original fitting bounds and all selected endpoint nuisance
coordinates are available in the receipts. Root Gaussian, A7K composition, and
bound-release profile sensitivity are reserved for assessment at the candidate
interval endpoints.

[All curve values](profile_grid.tsv) · [Mechanism summary](mechanism_summary.tsv) ·
[Baseline vector](baseline_vector.tsv) · [Both exact zero fits](exact_zero_summary.tsv) ·
[Four-parameter release gate](pilot_release.tsv) · [Source and file hashes](source_receipt.tsv)

![All 14 J profile objective traces](joint_profile_objectives.png)

## Paired bootstrap engineering test

Three pseudo datasets at each of four observed-data restricted J endpoints were
generated with the registered sample schedule, count covariance, count-dependent
death classifier, and specimen-specific chromosome sample sizes. The paired
restricted and unrestricted fits reused the same data and exact donor constraints.
The independent full joint log likelihood difference and the objective difference
agreed on all 12 simulated datasets at both fixed and varying parameters. Exact-null
carriers remained inactive where the target was fixed to zero.

**{valid}/12 pairs passed both arm and nesting gates.** {delta}{hypothesis}
An unresolved pair retains both fitted candidates and does not enter a nominal
critical value. These deliberately small samples test the numerical pipeline.

| Generating restriction | Fixed value | Verified pairs | Median CPU seconds per pair | Range |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(references)}

The observed median was **{p50:.1f} CPU seconds per pair**, with 75th percentile
**{p75:.1f}**, 90th percentile **{p90:.1f}** and maximum **{max(time):.1f}**. This includes
retained unresolved pairs and recorded reserve fits. CPU is summed for paired
work and can differ from elapsed time on the local multi-process run.
The planning allowance of {planning:.1f} CPU seconds per pair uses the observed
75th percentile plus 25%. Twelve pairs give a preliminary cost model. Solver work
will also vary by anchor, optimizer quality and dataset tail.

[Paired receipts and seeds](bootstrap_pilot_pairs.tsv) · [Measured cost](runtime_summary.tsv) ·
[Budget projections](calibration_budget.tsv) · [A7K baseline residuals](a7k_baseline_residuals.tsv)

## Batch calibration gate

The budget scenarios include 199 paired refits at each of the four engineering
anchors. The unresolved rate in this pilot calls for an optimizer reliability
check on the saved seeds before a calibration batch. An expansion to 999 replicates per candidate would require candidate-specific profile crossing brackets and its own cost accounting.
Quantiles and tail classifications must bound unresolved fits, retain all attempted
replicates and account for Monte Carlo error. A numerical failure fraction that
prevents classification at the 5% tail requires optimizer recovery or an unresolved
interval endpoint. A7K fit discrepancy and nuisance alternatives remain part of the
calibration sensitivity assessment. At the current J baseline, the four A7K death
errors in fitted-minus-observed percentage points range from {min(a7k_errors):.2f}
to {max(a7k_errors):.2f}. Any eventual nominal support sets concern
individual parameters conditional on the observation law and initialization
design. This run supplies the curves and benchmark for that next decision.

All fits in this release were local. Batch calibration and HPC submission are
awaiting approval.
"""
    (out/"REPORT.md").write_text(report)
    print(f"{len(profile)} grid rows: {dict(counts)}. Bootstrap pairs verified: {valid}/12.")
    print(f"Median pair CPU {p50:.2f} s, planning rate {planning:.2f} CPU seconds per pair.")


if __name__=="__main__":
    collect(Path(sys.argv[1]))
