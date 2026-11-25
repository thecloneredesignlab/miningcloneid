"""Verify saved certificates and compile the registered parameter verdicts."""

import argparse
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
from time import process_time

from flint import ctx

import joint_rank as rank
from algebra import decode, independent_minor_check
from paths import FROZEN_DESIGN, PROFILES, output_directory

RUNS = {
    "free_mixtures": ("free_384", "free_512_second"),
    "fixed_proxy_inheritance": ("inherit_512", "inherit_384_second"),
    "gaussian_root_inheritance": ("gaussian_384",),
    "a7k_tv_inheritance": ("tv_384",),
    "elapsed_inheritance": ("elapsed_inherit_384",),
    "elapsed_free_mixtures": ("elapsed_free_384",),
    "free_mixtures_free_viable_seed": ("free_seeds_384",),
}


def verify_archive(out):
    ctx.prec = 512
    checks = []
    pathways = {folder for folders in RUNS.values() for folder in folders}
    pathways.update(("pilot_free", "pilot_inherit_validated", "free_512_same_root"))
    for folder in sorted(pathways):
        for file in sorted((out/folder).glob("*_minors.tsv")):
            certificates = defaultdict(list)
            for row in rank.read(file):
                certificates[row["certificate"]].append(row)
            for name, rows in certificates.items():
                started = process_time()
                matrix,bound = rank.verify_minor_rows(rows)
                verified = bound < 1
                if matrix.nrows() <= 40:
                    algorithm = "exact_rational_midpoint_inverse"
                    independent,error = independent_minor_check(matrix)
                else:
                    algorithm = "outward_determinant"
                    determinant = matrix.det()
                    independent = determinant.is_finite() and not determinant.contains(0)
                    error = determinant
                checks.append(dict(run=folder,mask=file.stem.removesuffix("_minors"),certificate=name,order=matrix.nrows(),
                                   replayed_bound=float(bound),alternate_method=algorithm,alternate_result=str(error),
                                   alternate_pass=independent,verified=verified,
                                   cpu_seconds=process_time()-started))
                assert verified,(folder,file,name,bound)
        print(f"{folder}: certificate replay complete",flush=True)
    rank.write(out/"certificate_verification.tsv",checks)


def precision_checks(out):
    comparisons = (
        ("pilot_inherit_validated","inherit_512"),
        ("pilot_free","free_512_same_root"),
    )
    receipts = []
    for first,second in comparisons:
        def by_coordinate(folder):
            rows = rank.read(out/folder/"CDK_minors.tsv")
            return {(r["certificate"],r["output"],r["coordinate"]):decode(r) for r in rows}
        p1,p2 = by_coordinate(first),by_coordinate(second)
        checked = p1.keys() & p2.keys()
        counts = defaultdict(int)
        failures = []
        for key in checked:
            counts[key[0]] += 1
            if not p1[key].overlaps(p2[key]):
                failures.append(key)
        assert not failures, failures[:3]
        assert counts["global_after_root_elimination"] > 0
        if first=="pilot_free":
            assert counts["SUM-159_NLS_2N_dp_seed"] == 132*132
        for certificate,n in sorted(counts.items()):
            receipts.append(dict(original=first,comparison=second,certificate=certificate,overlapping_entries=n,
                            precision_bits="384 vs 512",passed=True))
    rank.write(out/"precision_overlap.tsv",receipts)


def case_conclusion(row, candidate):
    rows = int(candidate["rows"])
    root_rank = int(candidate["root_column_rank_lower"])
    retained = int(candidate["effective_global_rank_lower"])
    global_columns = int(candidate["global_columns"])
    n_local = int(row["deterministic_signal_coordinates"])-global_columns
    assert n_local>=root_rank and rows==int(row["total_output_rows"])
    unresolved_roots = bool(candidate["root_blocks_without_full_pivot"])
    root_rank_upper = n_local if unresolved_roots else root_rank
    eff_rows = int(candidate["effective_row_count"])
    adjusted_upper = rows if unresolved_roots else root_rank_upper+min(global_columns,eff_rows)
    rank_lower = root_rank+retained
    rank_upper = min(rows,adjusted_upper,n_local+global_columns)
    nuisance_lower = root_rank + max(0,retained-14)
    nuisance_upper = min(rows,root_rank_upper+global_columns-14,
                         rows if unresolved_roots else root_rank_upper+eff_rows)
    assert rank_lower<=rank_upper and nuisance_lower<=nuisance_upper
    if retained==global_columns:
        status = "locally_identifiable"
        proof = "Full column rank after exact elimination of local initial-state coordinates in the selected rows"
    elif root_rank==rows:
        status = "locally_nonidentifiable"
        proof = "Initial-state coordinates span every output row, so the implicit function theorem gives exact local compensation"
    else:
        status = "unresolved"
        proof = "Certified lower ranks leave an untested target-specific rank or invariance claim"
    return dict(output_rows=rows,total_columns=n_local+global_columns,
                nuisance_rank_lower=nuisance_lower,nuisance_rank_upper=nuisance_upper,
                augmented_rank_lower=rank_lower,augmented_rank_upper=rank_upper,
                effective_global_rank=retained,shared_columns=global_columns,
                all_mechanisms=status,proof=proof,unresolved_roots=candidate["root_blocks_without_full_pivot"])


def collect_verdicts(out, design_dir):
    definitions = rank.read(design_dir/"case_register.tsv")
    summary, verdicts, sensitivity = [], [], []
    for row in definitions:
        group = row["case"].removesuffix("_"+row["mask"])
        case_records = []
        for folder in RUNS[group]:
            data = next(r for r in rank.read(out/folder/"summary.tsv") if r["mask"]==row["mask"])
            point = {r["parameter"]:Fraction(r["value"]) for r in rank.read(out/folder/"witness.tsv")}
            domains = [r for r in rank.read(design_dir/"parameter_roles_by_case.tsv")
                       if r["case"]==row["case"] and r["role"]=="mechanism"]
            for role in domains:
                low,high = map(Fraction,role["domain"][1:-1].split(", "))
                assert low < point[role["coordinate"]] < high,(folder,role)
            result = dict(case=row["case"],mask=row["mask"],run=folder,witness=data["witness"],
                          precision_bits=data["precision"],**case_conclusion(row,data))
            case_records.append(result)
        outcomes = sorted({r["all_mechanisms"] for r in case_records})
        result = case_records[0]
        for r in case_records:
            summary.append(r)
        witness_status = "; ".join(f"{r['witness']} {r['all_mechanisms']}" for r in case_records)
        for target in rank.MECHANISMS:
            global_status = ("nonunique_on_domain_containing_witness"
                             if result["all_mechanisms"]=="locally_nonidentifiable" else "unresolved_global")
            verdicts.append(dict(case=row["case"],parameter=target,local_status=result["all_mechanisms"],
                                 witness_status=witness_status,global_status=global_status,
                                 evidence_run=result["run"],proof_rule=result["proof"]))
        if len(outcomes)>1:
            sensitivity.append(dict(case=row["case"],outcomes="|".join(outcomes),note="Witness conclusions differ"))
    rank.write(out/"case_results.tsv",summary)
    rank.write(out/"mechanism_verdicts.tsv",verdicts)
    if sensitivity:
        rank.write(out/"witness_disagreements.tsv",sensitivity)
    by_id = {r["case"]:r for r in summary if r["witness"]=="original"}
    comparisons = []
    for row in rank.read(design_dir/"contrast_register.tsv"):
        first,second = by_id[row["smaller"]],by_id[row["larger"]]
        comparisons.append(dict(smaller=first["case"],larger=second["case"],stream_added=row["changed_readout"],
                                status_smaller=first["all_mechanisms"],status_larger=second["all_mechanisms"],
                                proof_contrast=first["proof"]+" | "+second["proof"]))
    rank.write(out/"readout_contributions.tsv",comparisons)
    return by_id


def combine_with_practical(out, structural):
    mask_by_readout = {"death_fraction":"CK", "karyotype":"CD", "corrected_count":"DK"}
    table = []
    for r in rank.read(PROFILES/"paired_readout_comparisons.tsv"):
        if r["fixed_value"]!="0.0" or r["comparison_verified"]!="True":
            continue
        mask = mask_by_readout[r["assay"]]
        with_case = structural["fixed_proxy_inheritance_CDK"]
        without = structural[f"fixed_proxy_inheritance_{mask}"]
        table.append(dict(parameter=r["parameter"],readout_added=r["assay"],fixed_value=0,
                          structural_scope="Local positive-interior mechanistic signal under fixed proxy inheritance",
                          structural_with=with_case["all_mechanisms"],structural_without=without["all_mechanisms"],
                          practical_with_contract=r["with_contract"],practical_without_contract=r["without_contract"],
                          practical_rise_with=r["with_delta"],practical_rise_without=r["without_delta"],
                          practical_difference=r["delta_rise_difference"]))
    assert len(table)==6
    rank.write(out/"structural_and_practical_comparison.tsv",table)
    combinations = []
    for key in ("free_mixtures_C","free_mixtures_CD","free_mixtures_CK","free_mixtures_DK","free_mixtures_CDK",
                "fixed_proxy_inheritance_CDK"):
        case = structural[key]
        if case["all_mechanisms"]=="locally_identifiable":
            status = "locally_identifiable_as_function_of_identified_mechanisms"
        elif case["all_mechanisms"]=="locally_nonidentifiable":
            status = "locally_nonidentifiable_for_nonconstant_mechanistic_functions"
        else:
            status = "unresolved"
        for function in ("Effective per-chromosome missegregation at specified oxygen and parental ploidy",
                         "Post-missegregation survival at specified parental ploidy and nonzero error count"):
            combinations.append(dict(case=key,function=function,status=status,proof=case["proof"]))
    rank.write(out/"combination_verdicts.tsv",combinations)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=output_directory,required=True)
    parser.add_argument("--design",type=Path,default=FROZEN_DESIGN)
    args = parser.parse_args()
    verify_archive(args.out)
    precision_checks(args.out)
    combine_with_practical(args.out, collect_verdicts(args.out, args.design))
