"""Audit the staged fit receipts and assemble the comparison tables."""
import csv
from pathlib import Path
import re

OUT=Path(__file__).resolve().parent
MODES={"group_050":"5%","group_025":"2.5%"}
BACKGROUNDS={"up":"Upward shift","down":"Downward shift"}
CELLS={"HH":"O1 and O2 both shifted upward","LH":"Only O1 shifted downward","HL":"Only O2 shifted downward","LL":"O1 and O2 both shifted downward"}
CONTRASTS={
    "O1_given_O2_high":"Swap O1 with O2 shifted upward",
    "O1_given_O2_low":"Swap O1 with O2 shifted downward",
    "O2_given_O1_high":"Swap O2 with O1 shifted upward",
    "O2_given_O1_low":"Swap O2 with O1 shifted downward",
    "O1_minus_O2":"O1 alone shifted downward versus O2 alone",
    "both":"Switch both from upward to downward",
    "interaction":"Paired difference-in-differences",
}


def rows(name):
    with (OUT/f"{name}.tsv").open() as handle:
        return list(csv.DictReader(handle,delimiter="\t"))


def validate(variant,mode,backgrounds):
    claimed=[]
    for phase,length in (("pilot",20),("expanded",80)):
        s=f"{variant}_{mode}_{phase}"
        fs=rows(f"fits_{s}")
        cells=("LH","HL") if variant=="shared" else tuple(CELLS)
        n=length*2*len(backgrounds)*len(cells)
        assert len(fs)==n
        ids={(r["background"],r["case_id"],r["cell"],int(r["replicate"])) for r in fs}
        assert len(ids)==n and ids=={(b,j,q,i) for b in backgrounds for j in ("1","2")
                                      for q in cells for i in (range(1,21) if length==20 else range(21,101))}
        assert all(r["target_sha256"]==r["reciprocal_sha256"] for r in fs)
        assert all((r["cells_4N_O1"],r["cells_4N_O2"])==("12","18") for r in fs)
        assert all(float(r["deviance_replay_error"])<1e-7 and float(r["max_seed_error"])<1e-7 for r in fs)
        receipt=(OUT/f"slurm_{s}.tsv").read_text().splitlines()
        assert len(receipt)==n//8 and all(r.split("|")[1:]==["COMPLETED","0:0"] for r in receipt)
        claimed.extend(fs)
    nval=sum(r["validated"]=="True" for r in claimed)
    return len(claimed),nval,max(float(r["deviance_replay_error"]) for r in claimed),max(
        float(r["max_seed_error"]) for r in claimed)


def render():
    gates=rows("independent_fit_gate")
    studies=[]
    for model in ("shared","split4N"):
        for mode in MODES:
            bs=tuple(BACKGROUNDS) if model=="shared" else tuple(r["background"] for r in gates
                                                           if r["mode"]==mode and r["run_split"]=="True")
            if bs:
                studies.append((model,mode,bs,validate(model,mode,bs)))
    text="""# Paired recovery after swapping individual A7K 4N O1 and O2 branches

The comparison proceeds in the order 5%, then 2.5%. Within each of two fixed 2N backgrounds, the simulated karyotype of 4N O1 or O2 is switched from an upward to a downward shift. For each seed, the same combination constructed through the reverse swap was verified to match exactly. Counts, A17 karyotypes, and branch sample sizes retain their original values, with 12 and 18 cells observed for 4N O1 and O2, respectively. Seeding counts remain fixed in each passage. Fits retain the within-passage calendar-day AR(1) model and re-estimate the original free parameters.

The primary analysis fits four new combinations at each perturbation size. Each 2N background has an O1-only and an O2-only swap. The original paired results supply controls in which both 4N branches shift upward or both shift downward. Each new combination and truth group expands from 20 to 100 datasets. Recovery means being within 0.10 of the truth on the normalized fitting scale. A combination passes if both truth groups recover in at least 90/100 datasets. When the two 4N perturbations are estimated separately, 2N still shares perturbation coefficients and all four combinations for that background are refitted under the same model.

## Recovery by combination

| Fitted model | Perturbation limit | 2N background | 4N karyotype | Recoveries in truth groups 1／2 | Both groups pass | Mean error in groups 1／2 |
|---|---|---|---|---:|---|---:|\n"""
    detail=[]
    for model,mode,bs,_ in studies:
        style="Shared 4N perturbation" if model=="shared" else "Independent 4N perturbations"
        a=rows(f"metrics_{model}_{mode}")
        assert len(a)==8*len(bs)
        paired=rows(f"paired_{model}_{mode}")
        assert len(paired)==800*len(bs)
        for m in a:
            bucket=[r for r in paired if all(r[key]==m[key] for key in ("background","case_id","cell"))]
            assert sorted(int(r["replicate"]) for r in bucket)==list(range(1,101))
            assert sum(r["recovered"]=="True" for r in bucket)==int(m["recovered"])
            assert all((r["recovered"]=="True")==(
                r["validated"]=="True" and abs(float(r["normalized_error"]))<=.10) for r in bucket)
        for b in bs:
            for q,label in CELLS.items():
                x,z=sorted((r for r in a if r["background"]==b and r["cell"]==q),key=lambda r:r["case_id"])
                rec=[int(r["recovered"]) for r in (x,z)]
                e="／".join(f"{float(r['bias']):+.4f}" for r in (x,z))
                text+=f"| {style} | {MODES[mode]} | {BACKGROUNDS[b]} | {label} | {rec[0]}／{rec[1]} | {'Pass' if min(rec)>=90 else 'Fail'} | {e} |\n"
        t=rows(f"contrasts_{model}_{mode}")
        assert len(t)==14*len(bs)
        detail.append(f"### {style}，{MODES[mode]}\n\n")
        detail.append("| 2N background | Truth group | Comparison | Change in mean error [95% interval] | Change in recovery (percentage points) [95% interval] | Entered／left tolerance |\n|---|---:|---|---:|---:|---:|\n")
        for r in t:
            err=f"{float(r['signed_delta']):+.4f} [{float(r['signed_lo']):+.4f}, {float(r['signed_hi']):+.4f}]"
            hit=f"{100*float(r['rate_delta']):+.0f} [{100*float(r['rate_lo']):+.1f}, {100*float(r['rate_hi']):+.1f}]"
            n=f"{r['improved']}／{r['reversed']}" if r["improved"] else "Difference-in-differences"
            detail.append(f"| {BACKGROUNDS[r['background']]} | {r['case_id']} | {CONTRASTS[r['contrast']]} | {err} | {hit} | {n} |\n")
        detail.append(f"\n[Metrics by combination](metrics_{model}_{mode}.tsv) · [Paired replicates](paired_{model}_{mode}.tsv) · [Paired changes](contrasts_{model}_{mode}.tsv)\n\n")
    text+="\n## Effects within the same seed\n\nA positive change in mean error indicates higher estimates. Intervals use resampling of the original paired seeds within each truth group. The difference-in-differences tests whether the effect of swapping one branch depends on the other branch's direction.\n\n"
    text+="".join(detail)
    text+="## Numerical verification\n\n"
    for model,mode,bs,s in studies:
        style="Shared 4N perturbation" if model=="shared" else "Independent 4N perturbations"
        n,good,replay,seed=s
        text+=(f"- {style} at {MODES[mode]} had {n} fits, of which {good} passed the numerical checks. "
               f"The maximum deviance replay error was {replay:.2e}, and the maximum seeding error was {seed:.2e} cells.\n")
    text+="\n[Criteria for independent perturbation checks](independent_fit_gate.tsv) · [Job records](slurm_jobs.tsv)\n"
    for mode,pct in MODES.items():
        name=f"design_replay_{mode}.tsv"
        if (OUT/name).exists():
            text+=f"\n[Original control replay and model agreement at {pct}]({name})\n"
    claims=OUT/"interpretation.md"
    if claims.exists():
        header,rest=text.split("\n\n",1)
        text=header+"\n\n"+claims.read_text().strip()+"\n\n"+rest
    assert all((OUT/p).is_file() for p in re.findall(r"\]\(([^)]+)\)",text))
    (OUT/"REPORT.md").write_text(text,encoding="utf-8")
    print("Verified",sum(s[3][0] for s in studies),"new fits,",len(studies),"model and amplitude stages.")


if __name__=="__main__":
    render()
