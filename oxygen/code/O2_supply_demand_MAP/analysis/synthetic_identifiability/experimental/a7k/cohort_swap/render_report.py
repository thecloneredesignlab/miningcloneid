"""Verify delivered comparisons and write the analysis report."""
from pathlib import Path
import csv

ROOT=Path(__file__).resolve().parent
MODES={"group_050":"5%","group_025":"2.5%"}
COMBINATIONS={"HH":"2N upward, 4N upward","LH":"2N downward, 4N upward",
              "HL":"2N upward, 4N downward","LL":"2N downward, 4N downward"}
CONTRASTS={"2N_given_4N_high":"Switch 2N from upward to downward with 4N upward",
           "2N_given_4N_low":"Switch 2N from upward to downward with 4N downward",
           "4N_given_2N_high":"Switch 4N from upward to downward with 2N upward",
           "4N_given_2N_low":"Switch 4N from upward to downward with 2N downward",
           "both":"Switch both cohorts from upward to downward",
           "interaction":"Interaction"}


def read(name):
    with (ROOT/f"{name}.tsv").open() as file:
        return list(csv.DictReader(file,delimiter="\t"))


def format_interval(a,b,scale=1.):
    lo,hi=(scale*float(x) for x in (a,b))
    if scale==100.:
        return f"[{lo:+.1f}, {hi:+.1f}]"
    return f"[{lo:+.4f}, {hi:+.4f}]"


def verified(mode):
    all_rows=[]
    replay=[]
    seed=[]
    for phase,n in (("pilot",80),("expanded",320)):
        rows=read(f"fits_{mode}_{phase}")
        assert len(rows)==n
        assert all(r["target_sha256"]==r["reverse_target_sha256"] for r in rows)
        assert all(tuple(int(r[f"n_{k}"]) for k in ("2N_O1","2N_O2","4N_O1","4N_O2"))
                   ==(19,15,12,18) for r in rows)
        counts=(ROOT/f"slurm_{mode}_{phase}.tsv").read_text().splitlines()
        assert len(counts)==n//8 and all(r.split("|")[1:]==["COMPLETED","0:0"] for r in counts)
        all_rows.extend(rows)
        replay.extend(float(r["dev_replay_error"]) for r in rows)
        seed.extend(float(r["max_seed_error"]) for r in rows)
    assert len({(r["case_id"],r["replicate"],r["branch"]) for r in all_rows})==400
    for j in ("1","2"):
        for g in ("2N","4N"):
            ss=sorted(int(r["replicate"]) for r in all_rows if r["case_id"]==j and r["branch"]==g)
            assert ss==list(range(1,101))
    check=read(f"crossfit_checks_{mode}")
    assert len(check)==800
    flagged=read(f"crossfit_rechecks_{mode}") if (ROOT/f"crossfit_rechecks_{mode}.tsv").exists() else []
    assert len(flagged)==sum(float(r["gain_from_cross_evaluation"])>.025 for r in check)
    assert all(r["validated"]=="True" for r in flagged)
    return dict(validated=sum(r["validated"]=="True" for r in all_rows),
                replay=max(replay),seed=max(seed),rechecked=len(flagged),
                changed=sum(r["changed"]=="True" for r in flagged))


def build():
    results={m:read(f"factorial_metrics_{m}") for m in MODES}
    contrasts={m:read(f"factorial_contrasts_{m}") for m in MODES}
    qa={m:verified(m) for m in MODES}
    s="""# Swapping both A7K branches by cohort and recovery of `p_misseg`

The primary analysis uses a 5% A7K perturbation limit, with a check at 2.5%. Each simulation uses the original paired seeds and one of the two original truths. Counts, A17 karyotypes, cells sampled per branch, and recorded seeding counts are fixed. Fits retain the original within-passage calendar-day AR(1) error model and shared A7K perturbations within each ploidy cohort, and re-estimate the original free parameters.

For a given seed, swapping the two 2N branches from an upward background or swapping the two 4N branches from a downward background produces identical mixed data (2N downward, 4N upward). The converse combination is symmetric. At each perturbation size, the two original directions are the controls and two new mixed combinations are fitted. Each combination and truth group starts with 20 datasets and expands to 100 after numerical checks. Recovery means being within 0.10 of the truth on the normalized fitting scale. Both groups in a given combination must reach at least 90/100.

## Conclusion

Switching both 4N karyotypes from upward to downward is the main contributor to reduced underestimation and increased recovery. At both 5% and 2.5%, when 4N is shifted downward, both truth groups pass the precision threshold under either 2N background. Downward shifts in 4N bring the mean errors closer to zero in both truth groups and increase recovery in both 2N backgrounds. In truth group 2, the increases are 7 and 9 percentage points at 5%, and 5 and 3 percentage points at 2.5%.

At 5% in truth group 2 with 2N shifted upward, switching 4N downward changes the mean error from −0.0445 to −0.0284 and reduces its standard deviation from 0.0535 to 0.0464. Also switching 2N downward brings the mean error closer to zero, increases the standard deviation back to 0.0516, and lowers recovery from 93% to 90%. At 2.5%, the corresponding recovery decreases from 92% to 90%. The 2N shift can therefore affect the amount of underestimation and the number of tail outcomes outside tolerance. The improvement in recovery from 4N persists in both 2N backgrounds. The attribution applies to the simulated karyotypes, truths, and noise settings in this analysis.

## Precision criteria and error distributions

| Perturbation limit | Simulated karyotype | Truth group 1, recoveries／100 | Truth group 2, recoveries／100 | Both groups pass | Group 1 mean error [SD] | Group 2 mean error [SD] |
|---|---|---:|---:|---|---:|---:|\n"""
    for mode,label in MODES.items():
        rr=results[mode]
        assert len(rr)==8
        for cell,name in COMBINATIONS.items():
            a,b=sorted((r for r in rr if r["configuration"]==cell),key=lambda r:r["case_id"])
            rate=[int(r["recovered"]) for r in (a,b)]
            err=[f"{float(r['signed_mean']):+.4f} [{float(r['sd']):.4f}]" for r in (a,b)]
            s+=f"| {label} | {name} | {rate[0]} | {rate[1]} | {'Pass' if min(rate)>=90 else 'Fail'} | {err[0]} | {err[1]} |\n"
    s+="""

## Paired comparisons in the two backgrounds

Each comparison switches only the two branches in the indicated ploidy cohort from an upward to a downward shift. A positive change in mean error indicates higher estimates. Paired recovery changes are in percentage points. The interaction measures how the effect of the same swap differs under upward and downward backgrounds in the other ploidy cohort. The 95% intervals use resampling of the original paired seeds within each group.

| Perturbation limit | Truth group | Paired comparison | Change in mean error [95% interval] | Change in recovery (percentage points) [95% interval] | Entered tolerance／left tolerance |
|---|---:|---|---:|---:|---:|\n"""
    for mode,label in MODES.items():
        assert len(contrasts[mode])==12
        for r in contrasts[mode]:
            name=CONTRASTS[r["contrast"]]
            pair=f"{r['improved']}／{r['reversed']}" if r["improved"] else "Difference-in-differences"
            err=f"{float(r['signed_delta']):+.4f} {format_interval(r['signed_lo'],r['signed_hi'])}"
            recovery=f"{100*float(r['rate_delta']):+.0f} {format_interval(r['rate_lo'],r['rate_hi'],100.)}"
            s+=f"| {label} | {r['case_id']} | {name} | {err} | {recovery} | {pair} |\n"
    s+="""\n[5% metrics by combination](factorial_metrics_group_050.tsv) · [5% paired replicates](paired_runs_group_050.tsv) · [5% swap contrasts](factorial_contrasts_group_050.tsv)\n\n[2.5% metrics by combination](factorial_metrics_group_025.tsv) · [2.5% paired replicates](paired_runs_group_025.tsv) · [2.5% swap contrasts](factorial_contrasts_group_025.tsv)\n\n## Numerical verification\n\n"""
    for m,label in MODES.items():
        a=qa[m]
        s+=(f"- At {label}, {a['validated']} of the 400 new fits passed numerical verification. "
            f"The largest deviance replay error was {a['replay']:.2e}, and the maximum numerical error in recorded seeding counts was {a['seed']:.2e} cells. "
            f"Cross-checking estimates from the four combinations against the other observations from the same seed required {a['rechecked']} refits and changed {a['changed']} verdicts."
            f" [Cross-checks](crossfit_checks_{m}.tsv)\n")
    s+="\nThe jobs in both phases completed at each perturbation size with exit code 0:0. [Job records](slurm_jobs.tsv)\n"
    (ROOT/"REPORT.md").write_text(s,encoding="utf-8")
    print("Verified two amplitudes and all four paired configurations.",flush=True)


if __name__=="__main__":
    build()
