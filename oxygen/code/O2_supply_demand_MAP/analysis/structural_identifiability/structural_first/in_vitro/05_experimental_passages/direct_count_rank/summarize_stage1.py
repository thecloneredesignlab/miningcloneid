"""Verify the saved certificates and prepare the ten-parameter evidence table."""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from flint import arb_mat, ctx

from certify_counts import COLUMNS, TARGETS, HERE, COUNTS, INHERITANCE, decode, independent_minor_check, read, write


def summary():
    ctx.prec=512
    expected=read(HERE/"point_selection.tsv")
    results={r["case"]:r for r in read(HERE/"certificates.tsv")}
    assert len(results)==len(expected)==19 and set(results)=={r["case"] for r in expected}
    for row in read(HERE/"source_provenance.tsv"):
        assert sha256(Path(row["path"]).read_bytes()).hexdigest()==row["sha256"],row["path"]
    measured=read(HERE/"observation_audit.tsv")[0]
    assert [int(measured[q]) for q in ("count_records","fixed_day0_records","post_seed_records","ols_outputs",
                                      "donor_constraints","rows_with_exact_audited_weights")]==[335,114,221,114,110,335]
    parts=defaultdict(lambda: defaultdict(dict))
    for row in read(HERE/"certified_minors.tsv"):
        parts[row["case"]][(int(row["row_no"]),row["record_id"])][row["column"]]=decode(row)
    assert set(parts)==set(results)
    verify=[]
    for case, record in results.items():
        mat=parts[case]
        chosen=sorted(mat)
        assert len(chosen)==16 and len({r[1] for r in chosen})==16
        assert all(set(mat[key])==set(COLUMNS) for key in chosen)
        enclosed=arb_mat([[mat[key][n] for n in COLUMNS] for key in chosen])
        good,bound=independent_minor_check(enclosed)
        assert good and bound<1 and int(record["rank_certificate"])==16
        assert int(record["n_log_count_rows"])==221
        assert int(record["n_mechanism_columns"])==14 and int(record["n_visibility_columns"])==2
        assert float(record["native_count_max_rel_error"])<1e-12
        assert float(record["native_gradient_max_scaled_error"])<1e-7
        feasible=int(record["violated_donor_constraints"])==0
        if case=="archived_real_fit":
            assert not feasible and float(record["strict_supply_minimum_cells_lower"])<0
            assert record["scope_of_rank"]=="analytic_map_extension_at_infeasible_point"
        else:
            assert feasible and float(record["strict_supply_minimum_cells_lower"])>0
            assert record["scope_of_rank"] in ("local_at_feasible_interior","restricted_local_on_box_face")
        assert float(record["minimum_visible_count"])>1
        assert float(record["minimum_logistic_gap"])>0
        assert float(record["minimum_pruning_gap"])>0
        verify.append(dict(case=case,rank=16,physical_feasibility=feasible,
                           exact_rational_replay_residual_upper=str(bound)))
    write(HERE/"certificate_validation.tsv",verify)
    physical=[r for r in results.values() if int(r["violated_donor_constraints"])==0]
    assert len(physical)==18
    by_point=defaultdict(dict)
    for row in read(HERE/"exact_points.tsv"):
        by_point[row["case"]][row["parameter"]]=Fraction(row["exact_value"])
    assert len({tuple(v[q] for q in COLUMNS) for v in by_point.values()})==len(expected)
    face=by_point["mu_hp_lower_face"]
    above=by_point["mu_hp_inward"]
    move=above["mu_hp"]-face["mu_hp"]
    assert face["mu_hp"]==Fraction(1,200) and move>0
    assert all(face[n]==above[n] for n in COLUMNS if n!="mu_hp")
    assert results["mu_hp_lower_face"]["declared_box_faces"]=="mu_hp:lower"
    assert results["mu_hp_inward"]["declared_box_faces"]==""
    assert results["mu_hp_lower_face"]["operator_branch_masks"]==results["mu_hp_inward"]["operator_branch_masks"]
    default_ols=next(r for r in read(INHERITANCE/"rank_summary.tsv")
                     if (r["observation"],r["dp_root"],r["visibility"],r["transfer"])==("ols","hist","shared","exact"))
    assert default_ols["certified_rank_lower"]=="16"
    rankref=results["original_interior"]
    assert [int(rankref[q]) for q in ("old_counts_interval_overlaps","old_count_J_interval_overlaps",
                                     "old_ols_J_interval_overlaps")]==[335,221*16,114*16]
    practical={r["parameter"]:r for r in read(COUNTS/"parameter_comparison.tsv")}
    evidence=[]
    groups=defaultdict(list)
    for row in physical:
        if row["scope"]=="nominal_profile_support":
            groups[row["profile_parameter"]].append(row["case"])
    assert sum(map(len,groups.values()))==14
    for p in TARGETS:
        points=sorted(groups[p],key=lambda case:by_point[case][p])
        assert points
        scopes=["Certified lower-bound face and inward direction"] if p=="mu_hp" else []
        if any(results[q]["declared_box_faces"] for q in points):
            scopes.append("Includes feasible points on parameter bounds")
        q=practical[p]
        evidence.append(dict(parameter=p,prior_ols_rank="16/16",direct_count_original_rank="16/16",
                             direct_count_near_real_fit_rank="16/16",named_profile_points=len(points),
                             profile_values_natural=",".join(str(float(by_point[k][p])) for k in points),
                             profile_case_ids=",".join(points),direct_count_profile_ranks=",".join("16/16" for _ in points),
                             boundary_scope="、".join(scopes),all_feasible_points_local="18/18",
                             actual_precision_under_primary_noise=q["counts_pass"],
                             global_status="open_for_this_target"))
    write(HERE/"ten_parameter_evidence.tsv",evidence)
    fit=results["archived_real_fit"]
    near=results["near_fit_interior"]
    face_out=results["mu_hp_lower_face"]
    inward=results["mu_hp_inward"]
    documentation=[
        "# Local structural identifiability check using 221 in vitro counts","",
        "**Conclusion.** At each of the 18 parameter points that strictly satisfy donor supply constraints, the 221 × 16 Jacobian has full column rank. "
        "The 16 jointly varying unknowns are 14 mechanism parameters and 2 shared dead-cell visibility fractions. "
        "The ten parameters of interest are locally structurally identifiable in the feasible neighborhood of each point.","",
        "## Observation map","",
        f"We rechecked 116 passage segments, 114 fixed seeding records, 221 post-seeding counts, and 110 donor supply constraints. "
        f"We recomputed all 335 OLS weights as exact fractions from the original calendar days. Fixed seeding contributes a constant term. "
        f"The original OLS slope is that constant plus a fixed linear combination of the post-seeding log counts. The largest absolute error when reconstructing the measured slopes was "
        f"{float(measured['recorded_ols_max_abs_error']):.3g}/day. "
        "At the original witness, the intervals for all 335 counts, 3,536 count Jacobian entries, and 1,824 OLS Jacobian entries overlap the earlier intervals. "
        "The full column rank of the original OLS map therefore carries over to the individual-count map under the same conditions. We also certified the count Jacobian directly.","",
        "## Near the experimental fit and at the boundary","",
        f"The saved numerical fit falls short of one donor constraint by {-float(fit['strict_supply_minimum_cells_lower']):.6g} predicted live cells. "
        "The analytic extension at that point has full column rank. The result for the feasible domain uses a nearby interior point that satisfies all supply constraints. "
        f"The interior point is at most {float(near['max_normalized_distance_from_fit']):.3g} normalized fitting units from the saved fit. "
        f"The minimum supply margin is {float(near['strict_supply_minimum_cells_lower']):.6g} predicted live cells. "
        f"The count likelihood-ratio statistic relative to the saved fit is {float(near['count_profile_lr']):.3g}.","",
        f"Setting `mu_hp = 0.005` exactly gives a minimum supply margin of "
        f"{float(face_out['strict_supply_minimum_cells_lower']):.6g} predicted live cells. "
        f"Increasing only `mu_hp` by {float(move):.3g} preserves donor feasibility, "
        f"with a minimum margin of {float(inward['strict_supply_minimum_cells_lower']):.6g}. Both points use the same probability-pruning branches. "
        "At the lower-bound point, strict margins for supply, logarithms, and model branches ensure that the analytic map extends to a neighborhood. Full column rank gives uniqueness within that neighborhood restricted to the feasible domain.","",
        "## Evidence for the ten parameters","",
        "| Parameter | Original certificate witness | Near the experimental fit | Support points from this parameter's profile |",
        "|---|---:|---:|---:|",
    ]
    for e in evidence:
        documentation.append(f"| `{e['parameter']}` | 16/16 | 16/16 | {e['named_profile_points']}, all 16/16 |")
    documentation.extend(["",
        "Profile points came from the nominal 95% support set of the earlier count fit. Donor supply was checked at each point with interval arithmetic. "
        "For each minor, an independent exact rational midpoint inverse was used to recompute the residual norm. "
        "The certification uses karyotype inheritance, both dead-cell pools, the original Hill exponent, the Euler calendar, and the original probability pruning. "
        "The evidence covers the listed points and their feasible neighborhoods. Global uniqueness over the specified parameter domain remains open. "
        "The 10% tolerance and 90% recovery rate are criteria for practical precision.","",
        "[Evidence for the ten parameters](ten_parameter_evidence.tsv) · [Certificates by point](certificates.tsv) · "
        "[Certified minors](certified_minors.tsv) · [Independent recomputation](certificate_validation.tsv) · "
        "[Exact parameter points](exact_points.tsv) · [Observation map audit](observation_audit.tsv)","",
    ])
    (HERE/"REPORT.md").write_text("\n".join(documentation),encoding="utf-8")
    print("Verified",len(physical),"feasible certificates,",len(verify),"minor replays, and",len(evidence),"parameter rows")


if __name__=="__main__":
    summary()
