"""Assemble the report from verified certificates and seal the delivery."""

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from importlib.metadata import version
from os.path import relpath
from pathlib import Path
import platform
import re
import subprocess
import sys

import joint_rank as rank
from collect_results import RUNS
from paths import ARCHIVED_RUN, FINITE_STATE, FIT, FROZEN_DESIGN, ROOT, SOURCE, SYNTHETIC, output_directory


def portable(path):
    path = path.resolve()
    return str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)


def verdict(row):
    table = {"locally_identifiable":"14/14 locally identifiable",
             "locally_nonidentifiable":"0/14 locally identifiable", "unresolved":"Unresolved"}
    return table[row["all_mechanisms"]]


def interval(row, name):
    a,b = row[name+"_lower"],row[name+"_upper"]
    return a if a==b else f"{a} to {b}"


def report(out, design_dir):
    all_results = rank.read(out/"case_results.tsv")
    base = {r["case"]:r for r in all_results if r["witness"]=="original"}
    second = {r["case"]:r for r in all_results if r["witness"]=="second"}
    checks = rank.read(out/"certificate_verification.tsv")
    assert len(base)==19 and len(second)==10 and all(r["verified"]=="True" for r in checks)
    alt_count = sum(r["alternate_pass"]=="True" for r in checks)
    verified_order = Counter((r["run"],r["mask"]) for r in checks)
    directions = rank.read(out/"free_512_second/count_compensation.tsv")
    assert all(r["all_enclosures_contain_zero"]=="True" for r in directions if r["output_kind"]=="C")
    roles = rank.read(out/"mechanism_verdicts.tsv")
    assert len(roles)==266 and len({(r["case"],r["parameter"]) for r in roles})==266
    labels = {
        "free_mixtures":"Free mixtures at 14 restart roots",
        "fixed_proxy_inheritance":"Strict inheritance with two fixed root proxies",
        "gaussian_root_inheritance":"Inheritance with four free Gaussian root coordinates",
        "a7k_tv_inheritance":"Inheritance with four free A7K tilt coordinates",
        "elapsed_inheritance":"Inheritance with A7K elapsed collection time",
        "elapsed_free_mixtures":"Free mixtures with A7K elapsed collection time",
        "free_mixtures_free_viable_seed":"Free mixtures and free viable seeding fractions",
    }
    names = {"C":"Corrected visible counts", "CD":"Counts and measured death", "CK":"Counts and karyotype",
             "DK":"Measured death and karyotype", "CDK":"Counts, measured death and karyotype"}
    lines = [
        "# Joint experimental in vitro structural identifiability",
        "",
        "**Steps 4–7 completion report | 22 November 2026**",
        "",
        "## Primary result",
        "",
        f"With the three verified readouts measured on their recorded schedule, the primary free-mixture result is "
        f"**{verdict(base['free_mixtures_CDK'])}**. The co-primary donor-inheritance result is "
        f"**{verdict(base['fixed_proxy_inheritance_CDK'])}**. The five-way comparisons below use the same specimen design "
        "and freeze the initial-state and observer assumptions within each row comparison.",
        "",
        "| Scored observations | Free starting mixtures | Strict donor inheritance |",
        "| --- | --- | --- |",
    ]
    for mask in ("CDK","C","CD","CK","DK"):
        a,b = base[f"free_mixtures_{mask}"],base[f"fixed_proxy_inheritance_{mask}"]
        lines.append(f"| {names[mask]} | {verdict(a)} | {verdict(b)} |")
    if all(base[f"free_mixtures_{mask}"]["all_mechanisms"]=="locally_nonidentifiable" for mask in ("C","CD")):
        lines.extend(["","With the 14 free mixtures, the count and count-plus-death maps each have "
                      "enough initial-state variation to reproduce every small change in their outputs. "
                      "Their nuisance-only Jacobians span all output directions."])
    if all(base[f"free_mixtures_{mask}"]["all_mechanisms"]=="locally_identifiable" for mask in ("CK","DK","CDK")):
        lines.extend(["","Adding the measured chromosome distributions breaks that compensation. "
                      "In the free-mixture cases, pairing these distributions with either count or death achieves local identification "
                      "of all 14 mechanisms. The death-plus-karyotype contract fixes shared count-residual loading at zero. "
                      "The joint all-three case remains the primary result."])
    if all(base[f"fixed_proxy_inheritance_{mask}"]["all_mechanisms"]=="locally_identifiable" for mask in names):
        lines.extend(["","Under inherited composition the recorded counts already identify all 14 mechanisms locally. "
                      "Measuring death and chromosomes can therefore affect fit and precision even when they add no new "
                      "mechanism directions to a locally full-rank count design."])
    lines.extend([
        "",
        "The structural results are local at the certified witnesses. The two primary contracts were evaluated at two distinct witnesses. "
        "The registered witnesses are feasible parameter vectors for the rank proof. The data overlay uses the archived fitted vector. "
        "Full injectivity across the positive parameter domain remains open. Wherever nuisance variation yields exact local compensation, "
        "global injectivity on a domain containing that neighborhood is excluded.",
        "",
        "[Full case table](case_results.tsv) · [All 14 parameters by case](mechanism_verdicts.tsv) · "
        "[Readout comparisons](readout_contributions.tsv)",
        "",
        "## What was observed",
        "",
        "The joint design uses 225 post-seeding corrected visible counts, 195 complete measured dead fractions, and six empirical "
        "autosome distributions with 104 cells. Two 4N A17 and four A7K events have all three measurements. Other events retain their "
        "available readouts. The histories condition on 118 recorded seeding totals, leave two auxiliary upstream seed totals fixed, "
        "and enforce the same 110 original donor supply times. Oxygen is prescribed. The main A7K model uses calendar day 4 "
        "and a paired sensitivity uses the recorded elapsed times.",
        "",
        "The 792 chromosome coordinates in the rank analysis are components of the model probabilities at six sampling events. "
        "The 104 measured cells determine the empirical histograms and their finite-sample precision.",
        "",
        "The count observer includes visible dead objects and estimates the live component from the model. "
        "The archived count and death exports lack a shared field and mask crosswalk. The joint observation law conditions the "
        "death classification on the measured log count and uses the original eligible denominator. The proof for the law is "
        "given alongside the deterministic signal certificate. The observation coverage uses six identified chromosome specimens. "
        "The 20 flow specimens and two A19 karyotypes await the collection and processing evidence listed in the sample ledger.",
        "",
        "[Joint experimental trajectories and histograms](joint_trajectories.pdf) · "
        "[Event-level coverage](../20261122_in_vitro_joint_structural_steps01_03_01/coverage_by_passage.tsv) · "
        "[Sample and assay gaps](../20261122_in_vitro_joint_structural_steps01_03_01/open_evidence.tsv)",
        "",
        "## Ranks and scope of the parameter evidence",
        "",
        "| Initial-state assumption | Observation mask | Outputs × unknowns | Nuisance rank | Joint rank | Mechanism verdict |",
        "| --- | --- | ---: | ---: | ---: | --- |",
    ])
    for group,label in labels.items():
        for mask in ("CDK","C","CD","CK","DK"):
            r = base.get(f"{group}_{mask}")
            if r is None:
                continue
            lines.append(f"| {label} | {mask} | {r['output_rows']} × {r['total_columns']} | "
                         f"{interval(r,'nuisance_rank')} | {interval(r,'augmented_rank')} | {verdict(r)} |")
    if all(r["all_enclosures_contain_zero"]=="False" for r in directions if r["output_kind"]!="C"):
        lines.extend([
            "",
            "An explicit count-preserving direction was reconstructed for each of inducible missegregation and ploidy-dependent "
            "survival in the longest 2N block. Both directions change predicted death and chromosome outputs. The full "
            "count-plus-death nuisance certificate establishes that other local adjustments can preserve those two readouts "
            "together. The full joint certificate accounts for all local adjustments to the initial mixtures.",
            "",
            "[Verified count-preserving directions](free_512_second/count_compensation.tsv)",
        ])
    lines.extend([
        "",
        "C denotes corrected visible counts, D measured death and K normalized autosome frequencies. "
        "The ranks are for the deterministic signal coordinates. The same mechanism verdicts extend to the specified count and "
        "conditional death observation law where applicable. The 5% A7K row tests the four tilt freedoms at the zero-tilt interior "
        "witness under the per-branch limit. The free viable-seeding row uses the declared 99% live inoculum witness and "
        "splits the remaining visible total equally between the two pools. These are local assumption sensitivities.",
        "",
        "The targets are " + ", ".join(f"`{k}`" for k in rank.MECHANISMS) + ". "
        "The four root mean and width coordinates are included in the Gaussian sensitivity. The original growth scale and "
        "karyotype smoothing width remain registered for their historical observers. They are absent in this deterministic "
        "experimental signal. Count scale, AR(1) correlation, shared count loading and beta-binomial dispersion are law parameters.",
        "",
        "The joint positive-interior results support the effective missegregation and survival functions wherever their "
        "mechanism coordinates are identified. At exact `p_misseg=0` the carrier `k_o_mis` disappears. At exact "
        "`buffer_beta=0` the exponent `buffer_n_exp` disappears. These are exact submodel invariances.",
        "",
        "[Parameter-role contract](../20261122_in_vitro_joint_structural_steps01_03_01/parameter_roles_by_case.tsv) · "
        "[Combination results](combination_verdicts.tsv) · [Certificate method and observation-law argument](method.md)",
        "",
        "## What the measurements contribute in the fitted data",
        "",
        "The plotted data overlay uses the validated archived joint fit. Each exact-zero comparison below uses its own fitted "
        "unrestricted reference. A larger objective rise means that the corresponding mechanism null is harder to compensate in that fit.",
        "",
        "| Mechanism at zero | All three, conditional death model | Without measured death | Without karyotype |",
        "| --- | ---: | ---: | ---: |",
    ])
    practical = {(r["parameter"],r["readout_added"]):r for r in rank.read(out/"structural_and_practical_comparison.tsv")}
    for target in ("p_misseg","buffer_beta"):
        d = practical[target,"death_fraction"]
        k = practical[target,"karyotype"]
        lines.append(f"| `{target}` | {float(d['practical_rise_with']):.2f} | "
                     f"{float(d['practical_rise_without']):.2f} | {float(k['practical_rise_without']):.2f} |")
    lines.extend([
        "",
        "The measured dead fractions strengthen both exact-zero contrasts in the existing fit. The nuclear karyotype substantially "
        "strengthens the survival null contrast. For inducible missegregation, the sign of the karyotype contribution changes "
        "across the checked parameter values. The profiles use the observed data and specified error model. The rank certificates "
        "use equality of ideal outputs and cover all 14 mechanisms. The fitted A7K predictions overestimate all four visible counts "
        "and dead fractions. The dead fraction errors are 3.70 to 4.64 percentage points in the archived primary fit. "
        "The provisional 5% A7K starting-mixture profile at zero survival buffering retains its previous optimization status.",
        "",
        "[Practical comparisons and exact-zero contrasts](structural_and_practical_comparison.tsv) · "
        "[Joint fit residuals](../20261122_experimental_joint_fit_step02_01/REPORT.md) · "
        "[Limited profile review](../20261122_experimental_joint_identifiability_limited_review_01/REPORT.md)",
        "",
        "## Verification and completion",
        "",
        f"The archived count controls reproduce rank 221 in the free-mixture nuisance maps and rank 16 in the inherited count map. "
        f"The joint run independently checks output assembly and derivatives using complex-step propagation. "
        f"The 384-bit and 512-bit minor entries overlap in the registered precision checks. "
        f"All {len(checks)} saved minors passed an independent serialized-data inverse-bound replay. "
        f"The alternate rational or determinant check passed on {alt_count} minors. "
        f"The primary joint cases have {verified_order['free_384','CDK']} saved minor certificates for the first free-mixture witness "
        f"and {verified_order['inherit_512','CDK']} for the first inherited witness.",
        "",
        "The report covers all 19 registered cases and 266 parameter-case entries. Final structural interpretation is scoped to "
        "the verified assay definitions, actual masked schedule, feasible witness neighborhoods and mathematical model equations. "
        "The analysis used local single-thread propagation and certificate checks. The joint prediction plots use the previous "
        "fit. No new fitting or HPC work entered this release.",
        "",
        "[Control replays](control_replays.tsv) · [Precision overlap](precision_overlap.tsv) · "
        "[Minor verification](certificate_verification.tsv) · [Resource receipt](resources.tsv) · "
        "[Reproduction instructions](reproduce.md) · [Claim-to-evidence register](verification_receipt.tsv)",
    ])
    if len(checks)!=alt_count:
        lines.insert(-4,"Alternate-method failures are recorded separately from successful outward inverse-bound checks.")
    report_text = "\n".join(lines)+"\n"
    link_targets = {
        "../20261122_in_vitro_joint_structural_steps01_03_01/": design_dir,
        "../20261122_experimental_joint_fit_step02_01/": FIT,
        "../20261122_experimental_joint_identifiability_limited_review_01/":
            SYNTHETIC / "20261122_experimental_joint_identifiability_limited_review_01",
    }
    for prefix, target in link_targets.items():
        report_text = report_text.replace(f"({prefix}", f"({relpath(target, out)}/")
    report_text = report_text.replace("(method.md)", f"({relpath(ARCHIVED_RUN / 'method.md', out)})")
    report_text = report_text.replace("(reproduce.md)", f"({relpath(SOURCE / 'README.md', out)})")
    (out/"REPORT.md").write_text(report_text)


def resources(out):
    folders = {f for names in RUNS.values() for f in names}
    folders.update(("pilot_inherit_validated","pilot_free","free_512_same_root"))
    if (out/"pilot_inherit/summary.tsv").is_file():
        folders.add("pilot_inherit")
    rows = []
    for folder in sorted(folders):
        row = rank.read(out/folder/"summary.tsv")[-1]
        rows.append(dict(run=folder,precision_bits=row["precision"],cpu_seconds=row["total_cpu_seconds"],
                         wall_seconds=row["wall_seconds"],peak_rss_bytes=row.get("peak_rss_bytes",""),
                         stage="engineering_or_precision_pilot" if folder.startswith("pilot") or folder.endswith("same_root") else "registered_cases"))
    archived_log = (out/"archived_free_control.log").read_text()
    match = re.search(r"^real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys\s+([\d.]+)", archived_log, re.M)
    if match is None:
        match = re.search(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", archived_log)
    assert match, "Archived free count control is missing its timing receipt"
    rows.append(dict(run="archived_free_control",precision_bits=768,
                     cpu_seconds=float(match[2])+float(match[3]),wall_seconds=match[1],peak_rss_bytes="",stage="archived_control"))
    rank.write(out/"resources.tsv",rows)


def seal(out, design_dir):
    evidence = [
        design_dir/"SHA256SUMS", design_dir/"event_mask.tsv", design_dir/"passage_graph.tsv", design_dir/"case_register.tsv",
        rank.FIT/"joint_kernel.cpp", rank.FIT/"readout_predictions.tsv", rank.FIT/"karyotype_distributions.tsv",
        rank.FIT/"a7k_event_contract.tsv", rank.ROOT_FILE, rank.ALIGN/"cell_alignment.tsv",
        Path(rank.inherit.__file__), Path(rank.inherit.cert.__file__),
        Path(rank.inherit.cert.model.__file__), FINITE_STATE/"algebra.py",
        rank.STRUCT/"20261118_culture_experimental_counts_local_01/certified_minors.tsv",
    ]
    evidence.extend(path for path in SOURCE.iterdir() if path.is_file() and path.suffix in (".py", ".sh", ".R", ".md"))
    rank.write(out/"source_receipts.tsv",(dict(path=portable(p),sha256=sha256(p.read_bytes()).hexdigest()) for p in evidence))
    commit = subprocess.check_output(["git","-C",str(ROOT),"rev-parse","HEAD"], text=True).strip()
    rank.write(out/"environment.tsv",(dict(field=k,value=v) for k,v in
                                       (("python",sys.version.split()[0]),("numpy",version("numpy")),("scipy",version("scipy")),
                                        ("python_flint",version("python-flint")),("platform",platform.platform()),("git_head",commit),
                                        ("numerical_threads",1))))
    records = (
        ("4", "Original count controls", "control_replays.tsv", "221 free-mixture nuisance pivots and 16 inherited columns"),
        ("4", "Joint interface and tangent validation", "pilot_inherit_validated/derivative_validation.tsv",
         "Independent observer and derivative checks at the inherited witness"),
        ("4", "Largest free-mixture tangent block", "free_512_same_root/derivative_validation.tsv",
         "Independent checks of 132 composition directions and all shared columns"),
        ("4", "Repeated outward arithmetic", "precision_overlap.tsv", "384-bit and 512-bit matching joint-certificate entries"),
        ("5", "Verified case-specific ranks", "certificate_verification.tsv", "Every serialized minor replayed"),
        ("5", "Named local verdicts and global scope", "mechanism_verdicts.tsv", "Fourteen mechanisms for each of 19 cases"),
        ("5", "Effective mechanism functions", "combination_verdicts.tsv", "Consequences of uniqueness or exact compensation"),
        ("6", "Effects of adding measurements", "readout_contributions.tsv", "Matches the frozen case pairs"),
        ("6", "Joint empirical readouts", "joint_trajectories.pdf", "Recorded calendar and full autosome distributions"),
        ("6", "Practical comparison", "structural_and_practical_comparison.tsv", "Uses validated archived exact-zero fits"),
        ("7", "Observation coverage and remaining evidence", relpath(design_dir/"open_evidence.tsv", out),
         "Verified subset and specimen release evidence"),
    )
    assert all((out/f).exists() for _,_,f,_ in records)
    rank.write(out/"verification_receipt.tsv",(dict(step=s,claim=c,evidence=f,acceptance=a) for s,c,f,a in records))
    event_rows = rank.read(design_dir/"event_mask.tsv")
    dead_rows = [r for r in event_rows if r["modality"]=="death_fraction"]
    count_events = {r["event"] for r in event_rows if r["modality"]=="corrected_count"}
    series = defaultdict(set)
    for r in event_rows:
        if r["modality"]=="corrected_count":
            series[r["passage"]].add(float(r["model_day"]))
    multi = sum(len(days)>1 for days in series.values())
    eligible_min = min(int(r["denominator"]) for r in dead_rows)
    specimens = [r for r in event_rows if r["modality"]=="karyotype"]
    chromosome_min = min(int(r["sample_size"]) for r in specimens)
    assert len(dead_rows)==195 and all(r["event"] in count_events for r in dead_rows)
    assert multi>0 and eligible_min>1 and chromosome_min>0
    rank.write(out/"observation_law_checks.tsv",[
        dict(check="classification_events_with_count",observed=len(dead_rows),pass_condition="Every complete classification has same-event count"),
        dict(check="minimum_classifier_denominator",observed=eligible_min,pass_condition="Eligible denominator exceeds one"),
        dict(check="count_series_with_distinct_days",observed=multi,pass_condition="At least one within-series positive day gap"),
        dict(check="minimum_chromosome_sample_size",observed=chromosome_min,pass_condition="Every karyotype has at least one cell"),
    ])
    files = sorted(f for f in out.rglob("*") if f.is_file() and f != out/"SHA256SUMS")
    for f in files:
        if f.suffix in (".md",".tsv"):
            assert not re.search(r"[\u3400-\u9fff]",f.read_text()),f
        assert f.suffix not in (".json",".jsonl",".js",".java")
        if f.suffix==".md":
            for link in re.findall(r"\[[^\]]+\]\(([^)]+)\)",f.read_text()):
                assert (f.parent/link).exists(),(f,link)
    (out/"SHA256SUMS").write_text("".join(f"{sha256(f.read_bytes()).hexdigest()}  {f.relative_to(out)}\n" for f in files))
    print(f"Sealed {len(files)} local files",flush=True)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=output_directory,required=True)
    parser.add_argument("--design",type=Path,default=FROZEN_DESIGN)
    args = parser.parse_args()
    resources(args.out)
    report(args.out, args.design.resolve())
    seal(args.out, args.design.resolve())
