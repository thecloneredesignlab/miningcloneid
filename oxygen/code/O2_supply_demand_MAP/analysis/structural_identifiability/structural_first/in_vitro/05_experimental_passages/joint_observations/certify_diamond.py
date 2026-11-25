"""Certify the strict donor-inheritance signal at the original Panel A diamond."""
import argparse
import os
from collections import Counter
from time import process_time
import resource

for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"
import joint_rank as rank
from validate_joint import compare, forward, physical
from algebra import decode, encode, independent_minor_check
from paths import ARCHIVED_RUN, FROZEN_DESIGN, ROOT, output_directory
from flint import arb, arb_mat, ctx


def joined(key):
    return "|".join(map(str,key))


def parameters(out):
    point=rank.read(out/"fitted_point.tsv")
    raw={"k_clear":"1/10000"}
    alias={"eta_H_2N":"eta_H","eta_M_2N":"eta_M"}
    for r in point:
        if r["role"] in ("Mechanism","Signal observer"):
            raw[alias.get(r["parameter"],r["parameter"])]=r["exact_binary_value"]
    assert set(raw)==set(rank.COMMON_COLS)|{"k_clear"}
    return raw


def material(row):
    return row['output'],row['coordinate']


def assemble(out,bits):
    ctx.prec=bits
    folder=out/f"arb_{bits}"
    folder.mkdir(exist_ok=True)
    raw=parameters(out)
    native={r["output"]:float(r["native_diamond"]) for r in rank.read(out/"signal_replay.tsv")}
    donors=rank.read(out/"native_supply.tsv")
    reference, independent_margins=forward("inherit",raw)
    ref={joined(k):physical(k,v) for k,v in reference.items()}
    assert ref.keys()==native.keys()
    errors={k: abs(v/native[k]-1) if k.startswith("C|") else abs(v-native[k]) for k,v in ref.items()}
    assert max(errors.values())<1e-10
    abs_by_key={(node,tick):cells for node,tick,cells in independent_margins}
    intervals={}
    blocks,checks,pred,used=rank.propagate("inherit","fitted_diamond",parameter_values=raw,supply_receipts=intervals)
    margin_rows=[]
    for row in donors:
        key=(row["donor"],int(row["tick"]))
        absolute=intervals[key]
        demand=int(row["demand_cells"])
        relative=absolute/demand
        expected=float(row["native_relative"])
        assert relative>0
        assert abs(float(relative.mid())-expected)<1e-10*max(1,abs(expected))
        assert abs(abs_by_key[key]/demand-expected)<1e-10*max(1,abs(expected))
        margin_rows.append(dict(**row,available_cells=str(demand+absolute),independent_float_cells=abs_by_key[key],
            absolute_margin_cells=str(absolute),relative_margin=str(relative),**encode(absolute)))
    assert len(intervals)==len(abs_by_key)==len(margin_rows)==110
    assert min(margin_rows,key=lambda r:float(decode(r).mid())/int(r["demand_cells"]))["index"]=="49"
    rank.write(folder/"donor_supply_audit.tsv",margin_rows)
    rank.write(folder/"branch_checks.tsv",checks)
    assert len(checks)==120
    assert {joined(k) for k in pred}==native.keys()
    replay=[]
    for key,v in pred.items():
        k=joined(key)
        a=float(v)
        err=abs(a/native[k]-1) if key[0]=="C" else abs(a-native[k])
        assert err<1e-10
        replay.append(dict(output=k,native=native[k],interval_mid=a,interval_to_native_error=err,
                           independent_to_native_error=errors[k]))
    rank.write(folder/"signal_replay.tsv",replay)
    # Analytic branch derivatives and an independently coded complex-step observer.
    if bits==384:
        check_rows=compare("inherit",raw,None,blocks,pred,FROZEN_DESIGN)
        rank.write(folder/"derivative_validation.tsv",check_rows)
        assert all(r["passed"] for r in check_rows), [r for r in check_rows if not r["passed"]]
    ci,rows,keys,_,local,unresolved,_=rank.rank_case(blocks,"CDK","inherit")
    assert Counter(k[0] for k in keys)=={"C":225,"D":195,"K":792}
    assert tuple(ci)==tuple(range(18)) and local==0 and not unresolved
    rank.write(folder/"derivatives.tsv",(dict(output=joined(key),coordinate=n,**encode(v))
               for key,row in zip(keys,rows) for n,v in zip(rank.COMMON_COLS,row)))
    print("Validated supply, branches and signal at",bits,"bits.",flush=True)


def certify(out):
    ctx.prec=768
    branches=rank.read(out/"operator_branches.tsv")
    assert branches and all(r["branch_match"]=="True" for r in branches)
    derivative_checks=rank.read(out/"arb_384"/"derivative_validation.tsv")
    assert derivative_checks and all(r["passed"]=="True" for r in derivative_checks)
    precision_rows={}
    for bits in (384,512):
        precision_rows[bits]={material(r):decode(r) for r in rank.read(out/f"arb_{bits}"/"derivatives.tsv")}
    lo,hi=precision_rows.values()
    assert len(lo)==1212*18 and lo.keys()==hi.keys()
    assert all(a.overlaps(hi[k]) for k,a in lo.items())
    supply={}
    for bits in (384,512):
        supply[bits]=rank.read(out/f"arb_{bits}"/"donor_supply_audit.tsv")
    for a,b in zip(supply[384],supply[512]):
        assert (a["index"],a["donor"],a["tick"])==(b["index"],b["donor"],b["tick"])
        assert decode(a).overlaps(decode(b)) and decode(a)>0 and decode(b)>0
    old=rank.read(ARCHIVED_RUN/"inherit_512"/"CDK_minors.tsv")
    keys=list(dict.fromkeys(r["output"] for r in old))
    assert len(keys)==18
    evidence=[]
    for pick in ("Archived row selection","Fresh row selection"):
        if pick.startswith("Fresh"):
            allkeys=list(dict.fromkeys(k for k,n in hi))
            full=[[hi[k,n] for n in rank.COMMON_COLS] for k in allkeys]
            keys=[allkeys[i] for i in rank.qr_take(full,18)]
        success=True
        for bits,data in precision_rows.items():
            mat=arb_mat([[data[k,n] for n in rank.COMMON_COLS] for k in keys])
            _, bound=rank.precondition(mat)
            evidence.append(dict(selection=pick,bits=bits,residual_norm_upper=str(bound),
                                 rank18_certified=bool(bound<1)))
            if not bound<1:
                success=False
                break
        if success:
            break
    rank.write(out/"rank_certificate.tsv",evidence)
    if not success:
        raise ArithmeticError("Neither 18-row selection certified the fitted point; examine the null directions.")
    rank.write(out/"minor_rows.tsv",[dict(row=i,output=k) for i,k in enumerate(keys)])
    rank.write(out/"minor_coordinates.tsv",[dict(column=i,coordinate=n) for i,n in enumerate(rank.COMMON_COLS)])
    for bits,data in precision_rows.items():
        mat=arb_mat([[data[k,n] for n in rank.COMMON_COLS] for k in keys])
        rank.write(out/f"minor_{bits}.tsv",rank.matrix_rows("Diamond",mat,[(k,) for k in keys],rank.COMMON_COLS))
    rank.write(out/"precision_overlap.tsv",[dict(derivative_entries=len(lo),overlaps=len(lo),supply_enclosures=110,
        row_selection=pick, certified_rows=18, first_bits=384, second_bits=512)])
    print("18-column certificate",evidence,flush=True)


def verify(out):
    # Start a separate process and rebuild the matrices only from saved dyadic bounds.
    ctx.prec=768
    keys=[r["output"] for r in rank.read(out/"minor_rows.tsv")]
    names=[r["coordinate"] for r in rank.read(out/"minor_coordinates.tsv")]
    assert names==list(rank.COMMON_COLS) and len(set(keys))==18
    regular=arb_mat([[1,2],[3,4]])
    duplicate=arb_mat([[1,2],[1,2]])
    assert independent_minor_check(regular)[0] and not independent_minor_check(duplicate)[0]
    rank.write(out/"verifier_controls.tsv",[
        dict(control="Nonsingular two-column control",expected="Pass",observed="Pass"),
        dict(control="Duplicate-row two-column control",expected="Reject",observed="Reject")])
    results=[]
    for bits in (384,512):
        rows=rank.read(out/f"minor_{bits}.tsv")
        indexed={(r["output"],r["coordinate"]):decode(r) for r in rows}
        full={material(r):decode(r) for r in rank.read(out/f"arb_{bits}"/"derivatives.tsv")}
        assert len(rows)==len(indexed)==324
        for i,k in enumerate(keys):
            for j,n in enumerate(names):
                r=rows[i*18+j]
                assert (r["row"],r["column"],r["output"],r["coordinate"])==(str(i),str(j),k,n)
                assert indexed[k,n].contains(full[k,n])
        M=arb_mat([[indexed[k,n] for n in names] for k in keys])
        ok,residual=independent_minor_check(M)
        determinant=M.det()
        assert ok and not determinant.contains(0)
        results.append(dict(bits=bits,order=18,entries=324,
             rational_inverse_residual_upper=str(residual),determinant_enclosure=str(determinant),passed=ok))
    rank.write(out/"independent_verification.tsv",results)
    rank.write(out/"mechanism_verdicts.tsv",[dict(parameter=n,verdict="Locally identifiable at the fitted point",
        assumption="Strict donor inheritance", evidence="Verified 18-column signal minor") for n in rank.MECHANISMS])
    print("Independent exact-preconditioner and determinant rechecks passed.",flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=("384","512","certify","verify"))
    p.add_argument("out",type=output_directory)
    a=p.parse_args()
    elapsed=sum(float(r["cpu_seconds"]) for f in a.out.glob("*_resources.tsv") for r in rank.read(f))
    available=int(7200-elapsed-15)
    assert available>0 and (a.out/"freeze_resources.tsv").is_file()
    resource.setrlimit(resource.RLIMIT_CPU,(available,available))
    start=process_time()
    try:
        if a.stage.isdigit(): assemble(a.out,int(a.stage))
        elif a.stage=="certify": certify(a.out)
        else: verify(a.out)
    finally:
        target=a.out/f"{a.stage}_resources.tsv"
        previous=rank.read(target) if target.exists() else []
        rank.write(target,previous+[dict(cpu_seconds=process_time()-start)])
