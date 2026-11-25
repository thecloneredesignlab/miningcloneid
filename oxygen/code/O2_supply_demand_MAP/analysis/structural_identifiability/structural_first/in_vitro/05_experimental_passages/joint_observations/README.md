# Joint experimental culture observations

This directory is the working source for Steps 1–7 of the in vitro joint structural identifiability analysis. It covers the verified visible counts, measured death fractions and chromosome distributions at the recorded sampling events. The two primary initial-state assumptions allow free restart mixtures and strict donor inheritance.

| Step | Entry | Evidence written |
| --- | --- | --- |
| 1–3 | `freeze_design.py` | Specimen decisions, passage graph, case register, frozen parameter roles and archived endpoint replay |
| 4 | `validate_joint.py` | Archived count controls, observer fixtures and independent complex-step derivative checks |
| 5 | `joint_rank.py` | Feasible witnesses, interval minors, effective derivatives and case ranks |
| 5–6 | `collect_results.py`, `plot_joint.R` | Verified certificates, parameter verdicts, readout contrasts and joint trajectory PDF |
| 7 | `finalize_release.py` | Report, evidence and resource receipts, and output checksums |

## Inputs and runtime

Use Python 3.12 with NumPy, SciPy and python-flint 0.8.0, plus base R. The runner uses one numerical thread. Set `CULTURE_PYTHON` to the interpreter path when these dependencies are installed in a separate environment. Steps 1–3 replay the existing macOS joint library and verify the original source locators and hashes. Those locators include historical local archives outside this checkout.

The archived data and certificates are under `oxygen/results`. The original Steps 1–3 evidence is in `synthetic_identifiability/20261122_in_vitro_joint_structural_steps01_03_01`. The original Steps 4–7 evidence, method and source snapshots are in `synthetic_identifiability/20261122_in_vitro_joint_structural_steps04_07_01`. The deliberate archived dependencies include the joint fit and observation freeze, the earlier experimental count controls, and the existing practical profile comparisons. The interval state operator and its independent counterpart come from the existing in vitro source tree.

## Reproduce all seven steps

From the repository root choose a fresh output directory. Set the script and output paths before running the entry point.

```sh
joint=oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/in_vitro/05_experimental_passages/joint_observations
out=oxygen/results/synthetic_identifiability/joint_structural_replay
bash "$joint/run_cases.sh" all "$out"
(cd "$out" && shasum -a 256 -c SHA256SUMS)
```

The runner writes and seals `OUTPUT_DIR/design/` first, executes the old controls and new structural cases, then collects and seals the report. Each case can resume from its `summary.tsv` in the same unsealed run directory. `case_design.sha256` verifies that resumption and collection use the recorded case design. The original result packages contain the source snapshots used for those releases. Use a separate output directory for any replay or regeneration.

## Run by stage

```sh
bash "$joint/run_cases.sh" prepare "$out"
bash "$joint/run_cases.sh" certify "$out"
bash "$joint/run_cases.sh" collect "$out"
```

`certify` and `collect` use `OUTPUT_DIR/design/` when present and otherwise read the existing frozen Steps 1–3 design. Pass an explicit design directory as the third argument to either stage, or to `all` to reuse an existing freeze and start at Step 4. The input design directory must contain its recorded `SHA256SUMS` for final source receipts. `collect` runs from completed case folders and produces new reports in the specified unsealed output directory.

To validate one case against the published study without running all 19, use a separate output directory and the existing frozen design.

```sh
"${CULTURE_PYTHON:-python3}" "$joint/joint_rank.py" \
  --case inherit --bits 384 --masks CDK --validate \
  --out "$out/inherited_check"
```

The optional `--design` argument selects a regenerated design. The standalone `plot_joint.R` takes the output directory, design directory and archived fit directory in that order. Case table and minor values are reproducible for the pinned inputs and settings. Runtime receipts and plot metadata can vary.
