"""Keep the strongest verified fit per dataset and put profiles on one reference."""
from collections import defaultdict
import math

from inference import THRESHOLD, likelihood_ratio
from practical import HERE, MECHANISMS, read, write
from recovery import pilot_gate


def combine_recovery(root=HERE):
    groups=defaultdict(list)
    for name in ("recovery_fits.tsv","recovery_repairs.tsv","pilot_bank_a.tsv","pilot_bank_b.tsv",
                 "pilot_near_branch_confirmation.tsv","pilot_near_branch_confirmation_late.tsv",
                 "expanded_refits.tsv","pilot_near_branch_confirmation_expansion.tsv",
                 "pilot_near_branch_confirmation_expansion_b.tsv",
                 "strict_recovery_repair_a.tsv","strict_recovery_repair_b.tsv"):
        path=root/name
        if path.exists():
            for r in read(path):
                groups[r["scenario"],int(r["dataset"])].append((name,r))
    merged=[]
    for (scenario,dataset), candidates in sorted(groups.items()):
        feasible=[(name,r) for name,r in candidates if float(r["min_aggregate_donor_margin_fraction"])>=-1e-8]
        source,winner=min(feasible or candidates,key=lambda x:float(x[1]["loss"]))
        conf=[(f,r) for f,r in feasible if r["qualified"]=="True"
              and likelihood_ratio(float(r["loss"]),float(winner["loss"]),114)<=.10]
        if conf:
            source,winner=min(conf,key=lambda x:float(x[1]["loss"]))
        merged.append(dict(winner,selected_run=source,all_batches=len(candidates)))
    write(root/"qualified_recovery.tsv",merged)
    for case in sorted({r["scenario"] for r in merged}):
        r=[z for z in merged if z["scenario"]==case and int(z["dataset"])<=20]
        print(case,len(r),sum(z["qualified"]=="True" for z in r))
    if sum(int(r["dataset"])<=20 for r in merged)==40:
        print("Expansion passes",pilot_gate(merged))
    return merged


def profile_envelope(root=HERE):
    reference=read(root/"main_constrained_best.tsv")[0]
    filed=defaultdict(list)
    for name in ("main_profiles.tsv","main_profiles_fast.tsv","main_profiles_refit.tsv",
                 "main_profiles_edge_refit.tsv","main_profiles_extension.tsv",
                 "main_profiles_final_crossings.tsv"):
        p=root/name
        if p.exists():
            for row in read(p):
                filed[row["profile_parameter"],float(row["fixed_normalized"])].append((name,row))
    if not filed:
        return []
    output=[]
    names=set()
    for (param,x), suggestions in sorted(filed.items()):
        names.add(param)
        valid=[(f,r) for f,r in suggestions if r["feasible"]=="True"]
        label,best=min(valid or suggestions,key=lambda q:float(q[1]["loss"]))
        choose=dict(best,source_profile_run=label)
        choose["lr"]=(likelihood_ratio(float(best["loss"]),float(reference["loss"]),114)
                      if best["feasible"]=="True" else math.inf)
        output.append(choose)
    prototype=output[0]
    reference_loss=float(reference["loss"])
    supports=[(r["source_profile_run"],r) for r in output if r["feasible"]=="True"]
    for filename in ("main_constrained_fits.tsv","main_constrained_refinement.tsv","profile_reference_refits.tsv"):
        path=root/filename
        if path.exists():
            supports.extend((filename,r) for r in read(path)
                            if float(r["min_aggregate_donor_margin_fraction"])>=-1e-8)
    by_point={(r["profile_parameter"],float(r["fixed_normalized"])):r for r in output}
    for source,candidate in supports:
        lr=likelihood_ratio(float(candidate["loss"]),reference_loss,114)
        if lr>THRESHOLD:
            continue
        for key in names:
            point=float(candidate[f"x_{key}"])
            old=by_point.get((key,point))
            if old is not None and float(old["lr"])<=lr:
                continue
            r={f:candidate[f] if f in candidate else prototype[f] for f in prototype}
            r.update(profile_parameter=key,fixed_normalized=point,lr=lr,feasible=True,
                     source_profile_run=f"feasible projection from {source}")
            by_point[key,point]=r
    output=list(by_point.values())
    for key in sorted(names):
        point=reference[f"x_{key}"]
        r={f:reference[f] if f in reference else prototype[f] for f in prototype}
        r.update(profile_parameter=key,fixed_normalized=point,lr=0.,feasible=True,
                 source_profile_run="strict final experiment fit")
        output=[old for old in output if (old["profile_parameter"],float(old["fixed_normalized"]))!=(key,float(point))]
        output.append(r)
    print("Profile parameters",len(names))
    write(root/"main_profiles_envelope.tsv",sorted(output,key=lambda r:(r["profile_parameter"],float(r["fixed_normalized"]))))
    return output


if __name__=="__main__":
    profile_envelope()
    combine_recovery()
