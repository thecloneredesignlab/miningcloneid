"""Verify and seal the first endpoint release from saved receipts only."""

import argparse
from collections import Counter
from pathlib import Path
from scan import cpu, read, write, digest
from endpoint_zero_review import ANCHORS, PREVIOUS, SOURCE
import calibrate_bootstrap as batch
from scipy.stats import binom


def check_bounds(row, values):
    n=int(row["frozen"])
    t=float(row["observed_statistic"])
    k=sum(v>=t for v in values)
    f=n-len(values)
    assert n==199 and len(values)==int(row["verified"]) and k==int(row["verified_exceedances"])
    assert f==int(row["unresolved"])+int(row["unrun"])
    lo,hi=(float(row[col]) for col in ("tail_lower_add_one","tail_upper_add_one"))
    lmc,hmc=(float(row[col]) for col in ("monte_carlo_95_lower","monte_carlo_95_upper"))
    assert (lo,hi)==((k+1)/200,(k+f+1)/200)
    if k: assert abs(binom.sf(k-1,n,lmc)-.025)<2e-10
    else: assert lmc==0.
    if k+f<n: assert abs(binom.cdf(k+f,n,hmc)-.025)<2e-10
    else: assert hmc==1.
    decision=("Excluded" if max(hi,hmc)<.05 else
              "Retained" if min(lo,lmc)>.05 else "Inconclusive")
    assert decision==row["decision"]
    return values


def checked_release(out):
    assert (out/"resources.tsv").exists()
    assert not (out/"SHA256SUMS").exists()
    batch.check_seal(out/"FREEZE.sha256")
    for manifest in [SOURCE/"SHA256SUMS",PREVIOUS/"SHA256SUMS"]:
        batch.check_seal(manifest)
    original=read(PREVIOUS/"pair_results.tsv")
    old={r["seed"]:r for r in original}
    zero_rows=read(out/"zero_seed_dispositions.tsv")
    patches={r["seed"]:r for r in read(out/"zero_pair_reviews.tsv")}
    assert len(patches)==12 and len(zero_rows)==398
    for r in zero_rows:
        prev=old[r["seed"]]
        if r["seed"] not in patches:
            assert prev["status"]==r["status"]=="verified" and prev["statistic"]==r["statistic"]
            continue
        latest=patches[r["seed"]]
        assert latest["prior_status"]==prev["status"]=="unresolved"
        assert latest["status"]==r["status"] and latest["statistic"]==r["statistic"]
        pair=Path(r["parameter"])/f"replicate_{int(r['replicate']):03d}"
        assert digest(out/"zero_review"/pair/"observations.tsv")==digest(PREVIOUS/"pairs"/pair/"observations.tsv")
        assert prev["observations_sha256"]==latest["observations_sha256"]
    zeros=read(out/"zero_calibration.tsv")
    fresh=read(out/"fresh_calibration.tsv")
    seeds=read(out/"fresh_seed_results.tsv")
    assert len(seeds)==len({r["seed"] for r in seeds})==796 and not (set(old)&{r["seed"] for r in seeds})
    assert len({(r["parameter"],r["fixed_value"]) for r in zeros+fresh})==6
    gates=[]
    attempts=0
    reserve=0
    for r in zeros:
        check_bounds(r,[float(v["statistic"]) for v in zero_rows if v["parameter"]==r["parameter"] and v["status"]=="verified"])
    for r in fresh:
        name=r["parameter"]
        assert float(r["fixed_value"])==ANCHORS[name]
        g,=read(out/f"gate_{name}.tsv")
        close,=read(out/f"close_{name}.tsv")
        subset=[s for s in seeds if s["parameter"]==name]
        records=batch.read_pairs(out,name)
        lookup={v["seed"]:v for v in records}
        assert len(records)==int(r["processed"])==int(close["processed"])
        for seed in subset:
            raw=lookup.get(seed["seed"])
            if not raw:
                assert seed["status"]=="Not run"
                continue
            receipt=batch.pair_receipt(out,raw)
            assert raw["observations_sha256"]==seed["observations_sha256"]
            assert receipt["status"]==seed["status"] and raw["statistic"]==seed["statistic"]
            attempts+=int(raw["optimizer_invocations"])
            reserve+=int(raw["reserve_invocations"])
        first=[p for p in records if int(p["replicate"])<=25]
        pass25=sum(p["status"]=="verified" for p in first)
        assert pass25==int(g["first_wave_verified"])
        admitted=(len(first)==25 and pass25>=23 and
                  float(g["projected_anchor_cpu_seconds"])<=7*3600 and g["integrity_gate"]=="True")
        assert (g["release_remaining"]=="True")==admitted
        if not admitted: assert len(records)<=25
        if len(records)<199 and admitted:
            assert int(r["unresolved"])>=10 or "CPU" in close["finish_reason"]
        gates.append(dict(**g,final_processed=len(records),finish_reason=close["finish_reason"]))
        check_bounds(r,[float(v["statistic"]) for v in records if v["status"]=="verified"])
    write(out/"gate_summary.tsv",gates)
    write(out/"calibration_summary.tsv",zeros+fresh)
    status=[]
    legacy=[r for r in read(PREVIOUS/"calibration.tsv") if r["parameter"] in {"n_O","k_o_mis"}]
    for r in read(out/"endpoint_register.tsv"):
        name=r["parameter"]
        tested=[f'{v["fixed_value"]}: {v["decision"]}' for v in zeros+fresh if v["parameter"]==name]
        inherited=[f'{v["fixed_value"]}: retained in prior calibration' for v in legacy if v["parameter"]==name]
        status.append(dict(**r,calibrated_points=" | ".join(tested+inherited) or "No calibrated anchor",
                           interval_release="Endpoint refinement deferred; listed classifications apply at the tested values"))
    write(out/"endpoint_status.tsv",status)
    source=read(out/"resources.tsv")[0]
    write(out/"final_audit_code.tsv",[dict(code=Path(__file__).name,sha256=digest(Path(__file__)))])
    check_cpu=cpu()
    # A margin covers the short final checksum/verification pass after the ledger is written.
    seal_allowance=30
    measured=float(source["measured_cpu_seconds"])+check_cpu
    inspection_allowance=120
    conservative=float(source["preflight_incident_allowance_seconds"])+seal_allowance+inspection_allowance
    charged=measured+conservative
    assert charged<=36*3600
    write(out/"final_resources.tsv",[dict(
        charged_core_hours=charged/3600, measured_core_hours_before_seal=measured/3600,
        zero_recovery_cpu_seconds=source["zero_recovery_and_register_cpu_seconds"],
        fresh_worker_cpu_seconds=source["fresh_worker_cpu_seconds"],
        release_initialization_coordinator_cpu_seconds=source["initialization_collection_and_coordinator_cpu_seconds"],
        final_audit_cpu_seconds_before_seal=check_cpu,preflight_allowance_seconds=source["preflight_incident_allowance_seconds"],
        seal_allowance_seconds=seal_allowance, status_monitor_and_inspection_allowance_seconds=inspection_allowance, remaining_approved_core_hours=36-charged/3600,
        total_approved_core_hours=36,fresh_optimizer_dispatches=attempts,fresh_reserve_dispatches=reserve,
        zero_review_dispatches=sum(int(r["recovery_fits"]) for r in patches.values()))])
    brief=["# First conditional support calibration release","",
           "The endpoint register covers all 14 mechanisms under strict donor inheritance. This release tests two exact zeros and four prespecified fixed values. It uses corrected visible counts, conditional death classifications, and six verified full karyotype histograms.","",
           "| Parameter | Fixed value | Observed rise | Verified / unresolved / unrun | Exceedances | Add-one bounds | 95% Monte Carlo envelope | Decision |",
           "|---|---:|---:|---:|---:|---|---|---|"]
    for r in zeros+fresh:
        brief.append(f'| {r["parameter"]} | {float(r["fixed_value"]):.8g} | {float(r["observed_statistic"]):.6g} | {r["verified"]} / {r["unresolved"]} / {r["unrun"]} | {r["verified_exceedances"]} | {float(r["tail_lower_add_one"]):.3f}–{float(r["tail_upper_add_one"]):.3f} | {float(r["monte_carlo_95_lower"]):.5f}–{float(r["monte_carlo_95_upper"]):.5f} | {r["decision"]} |')
    brief+=["","## Frozen-data recovery","",f'The original two zero anchors contributed 398 frozen pairs. Recovery revisited all 12 unresolved pairs. {sum(v["status"]=="verified" for v in patches.values())} met the original agreement and nesting gates after recovery. The old draws, fits and audit are preserved. The inactive carrier is k_o_mis at zero inducible missegregation and buffer_n_exp at zero ploidy-dependent survival.',"",
            "## First-wave release decisions","","| Parameter | Passed / attempted | Projected core-hours at the gate | Disposition |","|---|---:|---:|---|"]
    for g in gates:
        brief.append(f'| {g["parameter"]} | {g["first_wave_verified"]} / {g["first_wave_pairs"]} | {float(g["projected_anchor_cpu_seconds"])/3600:.3f} | {g["finish_reason"]} |')
    brief+=["","Both arms used the same frozen draw. Two distinct qualifying starts within 0.02 of the feasible floor were required in each arm, with independently replayed joint likelihood and nesting at the selected endpoints and floors. Unresolved and unrun draws enter the upper tail bound as exceedances. Both probability bounds and their Monte Carlo envelope must be strictly on the same side of 0.05 for classification. A retained candidate and an excluded candidate can guide the next search. These six tests do not set all 14 interval endpoints.","",
            f'Local compute used {measured/3600:.5f} measured core-hours before sealing. The budget is charged {charged/3600:.5f} core-hours including conservative startup, status inspection and seal allowances. {36-charged/3600:.5f} core-hours remain in the approved first release. New anchors dispatched {attempts} optimizations including {reserve} reserve starts.',"",
            "[14-parameter status](endpoint_status.tsv) · [Six calibration results](calibration_summary.tsv) · [Zero recovery](zero_pair_reviews.tsv) · [First-wave gates](gate_summary.tsv) · [Fresh seeds](fresh_seed_results.tsv) · [Resource ledger](final_resources.tsv) · [Frozen inputs](FREEZE.sha256)"]
    (out/"REPORT.md").write_text("\n".join(brief)+"\n")
    count=batch.seal_results(out)
    assert cpu()-check_cpu<=seal_allowance
    print(f"Verified and sealed {count} release files. Measured final seal CPU {cpu()-check_cpu:.3f} s.")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    checked_release(args.output)
