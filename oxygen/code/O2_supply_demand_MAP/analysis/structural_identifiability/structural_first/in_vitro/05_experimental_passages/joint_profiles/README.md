# Three-readout culture profiles

Run from the repository root with the archived macOS arm64 Step 3 library
and Python 3.12.7, NumPy 2.3.1 and SciPy 1.16.0. Use a fresh output directory
under `oxygen/results/synthetic_identifiability`. All scientific processes
use one numerical thread each. The observed-data dispatch can run several
such local processes.

```sh
code=oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/in_vitro/05_experimental_passages/joint_profiles
out=oxygen/results/synthetic_identifiability/joint_profiles_replay
python -B "$code/scan.py" baseline "$out"
for parameter in p_misseg buffer_beta k_o_mis n_O; do
  python -B "$code/scan.py" pilot "$out" "$parameter"
  python -B "$code/bootstrap_pilot.py" "$out" "$parameter"
done
python -B "$code/run_profiles.py" "$out"
python -B "$code/run_profiles.py" "$out" --refine
python -B "$code/supply_envelope.py" "$out"
python -B "$code/summarize.py" "$out"
Rscript --vanilla "$code/plot_profiles.R" "$out"
```

The bootstrap entry point generates exactly three datasets per target. Its
optional `--reserve` revisits those seeds and allows one further local fit
per unresolved arm. Preserve the initial pair receipt separately before
using it. This code implements the engineering pilot. Batch calibration
requires a separate worked release decision and budget.

The profile grid uses the frozen J ranges and layer definitions. A profile
arm compares checked scores within J. Exact-zero carriers are removed when
their parent mechanism is explicitly fixed to zero. When an optimum lands
at a zero-parent face, the report records the inactive carrier and probes
release of that parent at multiple carrier values. All raw fits, actual
starts, boundary flags, feasibility gaps and lower candidate scores remain
in the tables. A failed feasible-start search does not certify infeasibility.

V4 and the source archives are read-only inputs here. Their source hashes and
the numerical reproduction receipt are saved in the new outputs. The plot
shows discrete objective rises and optimization evidence. Formal statistical
support strips require the separate calibration stage.

## Review the twelve frozen bootstrap pairs

```sh
review=oxygen/results/synthetic_identifiability/bootstrap_frozen_pair_replay
python -B "$code/recover_bootstrap_pairs.py" "$review"
```

This entry point reads the sealed pilot arrays and verifies their byte hashes
and event order before it starts any optimizer. It replays the old endpoints
and recovers selected arms using centered, diagonally scaled SLSQP with the
registered bounds and donor constraints. The final score and KKT checks use
the original coordinates. Each dispatched SLSQP invocation is charged against
the 31-invocation allocation, with a single subsequent reserve pass and a
2-core-hour CPU ceiling. Use a fresh result directory. The output files include
fit starts and endpoints, all embed-and-replay nesting checks, a CPU ledger
and the numerical assessment for a possible future calibration batch.

## Calibrate J at four recorded anchors

```sh
batch=oxygen/results/synthetic_identifiability/bootstrap_J_4anchors_B199_replay
python -B "$code/calibrate_bootstrap.py" test
python -B "$code/calibrate_bootstrap.py" freeze "$batch"
python -B "$code/calibrate_bootstrap.py" run "$batch"
```

Use a fresh directory. The first command checks seed uniqueness and the
failure-inclusive tail and order-statistic rules. The freeze command verifies
the sealed twelve-pair review and replays the four conditional generators and
the joint observed-data fit. It seals the 796 seed streams before any draw.
The run command manages up to four local single-thread workers. Each fits
its anchor's first 25 pairs and records all dispatches, endpoints and
qualification receipts. The remaining replicates launch only after the
first-wave numerical, provenance and CPU forecast checks. Each sample is
written once and used by both arms. Each unresolved arm may use one reserve
start. An anchor stops at ten unresolved pairs. Worker CPU is measured at
process exit and the initial and subsequent worker CPU allocations share
one 28-core-hour ceiling.

The run output includes the signed pair contrasts, seed outcomes including
unrun replicates, conservative add-one tail bounds, exact binomial Monte
Carlo bounds and critical-value rank bounds. The four results concern the
specified generating anchors under strict donor inheritance. Inspect the
recorded process exit status and partial outputs if a worker is interrupted.
The controller protects existing dispatch receipts from a second launch.

## First endpoint release

Use an empty path under the synthetic identifiability results directory for a replay.
The first command creates the 14-parameter register and recovers failed arms on the
12 saved exact-zero datasets. The second requires the completed recovery receipts,
freezes the four new generators and 796 seeds, and checks each anchor after its
first 25 pairs. A release to 199 requires at least 23 verified pairs, complete
provenance, and a forecast within that anchor's seven-core-hour allowance.
Four single-threaded processes run locally, including the coordinating worker.
The last command independently verifies the dispositions and the resource ledger,
writes the English report and seals the result files.

```sh
out=oxygen/results/synthetic_identifiability/joint_support_first_release_replay
python -Bu "$code/endpoint_zero_review.py" "$out"
python -Bu "$code/endpoint_bootstrap_release.py" run "$out"
python -Bu "$code/seal_endpoint_release.py" "$out"
```

The first release has a 36-core-hour budget. When an anchor fails its gate, the
unrun seeds stay in the tail bounds. The two exact-zero hypotheses use the prior
199 draws with the frozen-data recovery applied. The release contains no
additional 999-draw samples. An interrupted release retains its dispatch receipts;
use a fresh directory for an independent replay.

### Second endpoint calibration release

Use the same environment and run these commands in order from the repository
root. The first command checks the sealed sources and reviews four paired fits
on the original observations. The next command verifies the four observed
restricted fits and freezes 199 seeds for each accepted anchor.

```sh
out=oxygen/results/synthetic_identifiability/joint_support_second_release_replay
python -Bu "$code/second_endpoint_prepare.py" recover "$out"
python -Bu "$code/second_endpoint_prepare.py" prepare "$out" 1 2 3 4
python -Bu "$code/second_endpoint_workers.py" primary "$out"
python -Bu "$code/second_endpoint_workers.py" reserve "$out"
python -Bu "$code/seal_second_endpoints.py" "$out"
```

Each 25-pair wave needs at least 23 verified pairs and a projected cost within
seven core-hours to continue. The first four anchors test the lower midpoints
for `p_misseg` and `lam_max`, followed by the registered `buffer_beta` and
`buffer_smax` values. A complete parent run and sufficient projected capacity
can release at most two further anchors. The reservation is 25.46 core-hours,
including 22.96 for paired workers, 1.50 for preparation and 1.00 for
coordination and verification. The coordinator shares the four-process limit.
The final command recomputes the failure-inclusive decisions, writes both
endpoint directions for all 14 parameters and seals the files.

### V6 joint support intervals

Run from the repository root with the same Python environment. The local
controller waits for the second endpoint release, including its registered
reserve decisions, and verifies its full seal before starting new fits.

```sh
python -Bu "$code/v6_endpoints.py" test
python -Bu "$code/test_v6_endpoints.py"
python -Bu "$code/v6_endpoints.py" run
```

The controller owns `oxygen/results/synthetic_identifiability/20261124_joint_support_v6_01`.
Its lock protects an existing run. The progress table and controller log
record the current stage. Each batch freezes up to four anchors before paired
refitting, with one numerical thread per worker. Each anchor gets at most
seven paired calibration core-hours. Preparation and verification are recorded
separately. The first 25 draws must pass the original quality and cost gates
to continue to 199. Unresolved and unrun draws remain in the tail bounds.

Both directions of all 14 parameters use the strict donor-inheritance J
contract. A retained–excluded bracket must span at most 0.02 on the registered
coordinate. The register distinguishes sampled crossings, supported domain
contacts and sides still open after the midpoint and flanking checks. Each
sealed batch refreshes the English V6 report alongside V5. The final audit
records both endpoint directions, natural-unit bounds, point decisions,
resource use and the HTML digest. Inspect the last batch's receipts and
process state after an interruption before reconciling a partially written
draw. Each frozen draw belongs to precisely one anchor.

The first V6 controller writes log files alongside the slot directories.
Its original collection pattern also matches those logs. For the registered
release, the following continuation waits for the active owner, verifies the
first batch's saved slots, and collects the slots listed in the batch's
frozen seed manifest. The original worker, per-slot seals and CPU accounts
are reused. The handoff log is `collection_controller.log`.

```sh
python -Bu "$code/v6_batch_controller.py"
```
