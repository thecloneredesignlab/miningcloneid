"""Certify the actual post-inoculation count Jacobian on specified model branches."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import multiprocessing
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from flint import arb, arb_mat, ctx

from prepare_points import (COUNTS, HERE, INHERITANCE, MECHANISMS, TARGETS, VISIBILITY,
                            Configuration, model, natural_array, read, write)
sys.path.insert(0,str(INHERITANCE))
import inherit as intervals
import rank_cases
from algebra import decode, encode, exact, independent_minor_check

COLUMNS=(*MECHANISMS,*VISIBILITY)
K=intervals.cert.CAPACITY
REFERENCE=INHERITANCE/"hist"


def audit_map(k):
    with TemporaryDirectory() as work:
        design,demand,_=intervals.build_design(Path(work))
        root=intervals.root_histograms(Path(work))
    ids={r["record_id"]:r for r in k.slots}
    forward={r["passage_id"]:r for r in k.graph}
    assert set(forward)==set(design)
    assert len(k.graph)==116 and k.n_constraints==110
    assert sum(map(len,demand.values()))==110
    assert set(np.flatnonzero(k.capture))==set(k.indices)
    count_ids=[]
    oats={}
    for p in design.values():
        candidate=forward[p.id]
        assert p.cohort==candidate["cohort"] and Fraction(p.oxygen)==Fraction(candidate["oxygen_pct"])
        assert (p.parent or "")==candidate["parent_id"]
        assert str(p.donor_tick)==candidate["donor_step"] if p.parent else candidate["donor_step"] in ("","None")
        assert p.seed_count==int(candidate["fixed_live_cells"])
        if not p.samples:
            continue
        times=[Fraction(t)*intervals.DT for t,_,_ in p.samples]
        mean=sum(times)/len(times)
        denom=sum((t-mean)**2 for t in times)
        initial=None
        for day,(tick,key,w) in zip(times,p.samples):
            slot=ids[key]
            assert slot["passage_id"]==p.id and int(slot["euler_step"])==tick
            assert day==int(slot["calendar_day"])
            expected=(day-mean)/denom
            assert expected==w==Fraction(int(slot["ols_weight_numerator"]),int(slot["ols_weight_denominator"]))
            assert w==Fraction(k.weights[k.slots.index(slot)]).limit_denominator(10**8)
            if tick:
                count_ids.append(key)
            else:
                initial=w
        assert initial is not None and p.samples[0][0]==0
        oats[p.id]=initial
    count_source=read(COUNTS/"count_inputs.tsv")
    native_ids=[r["record_id"] for r in k.slots]
    assert native_ids==[r["record_id"] for r in count_source]
    assert len(ids)==335 and len(set(count_ids))==len(count_ids)==221 and len(oats)==114
    assert set(count_ids)=={native_ids[i] for i in k.indices}
    assert all((r["fitted_count"]=="True")==bool(cap) for r,cap in zip(count_source,k.capture))
    histories=read(REFERENCE/"dp_proxy_histograms.tsv")
    for cohort in ("2N","4N"):
        assert sum(int(r["cells"]) for r in histories if r["cohort"]==cohort)==20
        ref={int(r["chromosome_count"]):int(r["cells"]) for r in histories if r["cohort"]==cohort}
        assert all((v-exact(Fraction(ref[n],20))).contains(0)
                   for n,v in zip(intervals.cert.model.GRID,root[cohort].entries()))
    assert {p.id for p in design.values() if p.parent is None}==set(intervals.ROOT_IDS)
    weights=np.zeros((k.ols_n,k.n))
    weights[k.group[k.capture],np.arange(k.n)]=k.weights[k.capture]
    offsets=np.bincount(k.group[~k.capture],k.weights[~k.capture]*k.observed_logs[~k.capture],minlength=k.ols_n)
    measured=np.asarray([float(r["corrected_count"]) for r in count_source])
    reproduced=offsets+weights@np.log(measured[k.capture])
    error=max(abs(reproduced-k.ols_y))
    assert error<1e-12
    return design,demand,root,weights,offsets,dict(passages=len(design),donor_constraints=110,
               count_records=len(ids),fixed_day0_records=114,post_seed_records=len(count_ids),
               ols_outputs=len(oats),rows_with_exact_audited_weights=len(ids),
               recorded_ols_max_abs_error=error)


def verified_operator(p,oxygens):
    previous=intervals.cert.model.operator
    metadata=[]
    def record_operator(parameters, oxygen):
        result=previous(parameters,oxygen)
        gap=min(r["minimum_margin"] for r in result["masks"])
        assert gap>0
        metadata.append((result["mask_sha256"],gap))
        return result
    # The shared operator already checks the inequalities. Capture its margins.
    intervals.cert.model.operator=record_operator
    try:
        operators={o:intervals.operator(p,o) for o in oxygens}
    finally:
        intervals.cert.model.operator=previous
    return operators,metadata


def exact_forward(p,eta,design,demand,roots):
    operators,prune=verified_operator(p,sorted({p.oxygen for p in design.values()},key=Fraction))
    donor={}
    all_counts={}
    log_jac={}
    slope={}
    slope_jac={}
    supply=[]
    crowds=[]
    for seg in design.values():
        if seg.parent is None:
            pi=roots[seg.cohort]
            j=arb_mat(intervals.N,len(COLUMNS))
        else:
            x,j,L,_,_=donor[seg.parent][seg.donor_tick]
            pi,j=intervals.normalized(x,j,L)
        known=exact(Fraction(seg.seed_count,K))
        samples={t for t,_,_ in seg.samples}
        requested=demand.get(seg.id,{})
        ticks=samples | set(requested)
        values,crowd=intervals.propagate(pi*known,j*known,ticks,operators[seg.oxygen],
                                          1-exact(intervals.DT)*p["k_clear"],*eta)
        crowds.append(crowd)
        donor[seg.id]={t:values[t] for t in requested}
        for t,qty in requested.items():
            L=values[t][2]
            surplus=L*K-qty
            supply.append((seg.id,t,surplus,surplus/qty))
        if not seg.samples:
            continue
        s=arb(0)
        g=[arb(0) for _ in COLUMNS]
        for t,rid,w in seg.samples:
            _,_,_,visible,derivatives=values[t]
            key=(seg.cohort,rid)
            all_counts[rid]=visible*K
            s+=exact(w)*(visible*K).log()
            if t==0:
                assert (visible-known).contains(0) and all(v.contains(0) for v in derivatives)
                continue
            log_jac[key]=[v/visible for v in derivatives]
            for i,v in enumerate(log_jac[key]):
                g[i]+=exact(w)*v
        slope[seg.id]=s
        slope_jac[seg.cohort,seg.id]=g
    assert len(all_counts)==335 and len(log_jac)==221 and len(supply)==110
    assert len(slope)==len(slope_jac)==114
    assert all(s[2]>0 or s[2]<0 for s in supply)
    return all_counts,log_jac,slope,slope_jac,supply,crowds,prune


def overlap_reference(all_counts,jac,slopes):
    overlap_counts=overlap_log_j=overlap_slope_j=0
    for cohort in ("2N","4N"):
        name_map={n:n for n in MECHANISMS}
        name_map.update({f"eta_H_{cohort}":VISIBILITY[0],f"eta_M_{cohort}":VISIBILITY[1]})
        old_counts={r["record_id"]:arb(r["expected_visible_witness_cells"])
                    for r in read(REFERENCE/f"{cohort}_counts_replay.tsv")}
        for key,value in old_counts.items():
            assert all_counts[key].overlaps(value),(cohort,key)
            overlap_counts+=1
        for r in read(REFERENCE/f"{cohort}_counts_derivatives.tsv"):
            if r["column"] not in name_map:
                continue
            key=r["output_id"]
            log_derivative=decode(r)/(all_counts[key]/K)
            actual=jac[cohort,key][COLUMNS.index(name_map[r["column"]])]
            assert actual.overlaps(log_derivative),(key,r["column"])
            overlap_log_j+=1
        for r in read(REFERENCE/f"{cohort}_ols_derivatives.tsv"):
            if r["column"] not in name_map:
                continue
            key=r["output_id"]
            actual=slopes[cohort,key][COLUMNS.index(name_map[r["column"]])]
            assert actual.overlaps(decode(r)),(key,r["column"])
            overlap_slope_j+=1
    assert (overlap_counts,overlap_log_j,overlap_slope_j)==(335,221*16,114*16)
    return overlap_counts,overlap_log_j,overlap_slope_j


def certificate(matrix):
    # A nonzero column scaling preserves rank. Row selection uses floating
    # midpoints; the complete selected submatrix is certified with balls.
    powers=[]
    for j in range(len(COLUMNS)):
        m=max(abs(row[j].mid()) for row in matrix)
        mantissa,exponent=m.man_exp()
        powers.append(0 if m.is_zero() else exponent+abs(int(mantissa)).bit_length())
    scaled=[[v/arb((1,q)) for v,q in zip(row,powers)] for row in matrix]
    rows=rank_cases.row_selector(scaled)
    assert len(rows)==len(COLUMNS)
    minor=arb_mat([scaled[i] for i in rows])
    bound=intervals.cert.midpoint_preconditioned_bound(minor)
    good,independent=independent_minor_check(minor)
    assert bound<1 and good and independent<1
    return rows,powers,minor,bound,independent


def native_checks(k,p,logs,chosen,jac):
    c=Configuration(k)
    numeric=np.array([[float(a.mid()) for a in row] for row in jac])
    errors=[]
    # Direct selected-row adjoints and dense contrasts exercise the full
    # native count-likelihood chain, including the shared visibility columns.
    _,baseline,_=k.constrained(p,np.ascontiguousarray(logs),c,True)
    contrasts=np.random.default_rng(401548).normal(size=(3,k.n))
    for w in [*np.eye(k.n)[chosen],*contrasts]:
        _,gradient,_=k.constrained(p,np.ascontiguousarray(logs-w),c,True)
        assert np.all(np.isfinite(gradient))
        compare=c.gradient(gradient-baseline,np.ones(len(COLUMNS)))
        expected=w@numeric
        denominators=np.maximum(1e-10,abs(w)@abs(numeric))
        errors.append(max(abs(compare-expected)/denominators))
    err=float(max(errors))
    assert err<1e-7,err
    return err


def worker(spec):
    detail,values,precision=spec
    ctx.prec=precision
    k=model()
    c=Configuration(k)
    assert tuple(c.names)==COLUMNS
    p={n:exact(values[n]) for n in MECHANISMS}
    p["k_clear"]=exact(Fraction(1,10000))
    eta=tuple(exact(values[n]) for n in VISIBILITY)
    design,demand,roots,W,offset,design_check=audit_map(k)
    count,jac,slope,slope_jac,supply,crowd,pruning=exact_forward(p,eta,design,demand,roots)
    profile_slope=native_counts=diagnostics=None
    pp=natural_array(c,values)
    profile_slope,native_counts,diagnostics,_,_=k.predict(pp)
    logs=np.log(native_counts[k.capture])
    assert np.all(native_counts>1)
    native_slopes=offset+W@logs
    certified_counts=np.array([float(count[z["record_id"]].mid()) for z in k.slots])
    forward_error=max(abs(native_counts/certified_counts-1))
    slope_error=max(abs(native_slopes-np.array([float(slope[q].mid()) for q in k.ids])))
    native_ols_error=max(abs(native_slopes-profile_slope[:k.ols_n]))
    assert forward_error<1e-12 and slope_error<1e-12 and native_ols_error<1e-12
    bad=[s for s in supply if s[2]<0]
    if detail["case"]=="archived_real_fit":
        # Validate the status independently of the optimizer's feasibility tolerance.
        assert bool(bad)==(diagnostics[1]<0)
    else:
        assert not bad,(detail["case"],bad)
    strict_supply=min(supply,key=lambda s:float(s[2].lower()))
    limiting=min(supply,key=lambda s:float(s[3].lower()))
    assert abs(float(limiting[3].mid())-diagnostics[1])<1e-12
    row_ids=[(k.slots[i]["cohort"],k.slots[i]["record_id"]) for i in k.indices]
    matrix=[jac[key] for key in row_ids]
    at_reference=detail["case"]=="original_interior"
    comparison=overlap_reference(count,jac,slope_jac) if at_reference else ("","","")
    rows,powers,minor,bound,independent=certificate(matrix)
    gradient_error=native_checks(k,pp,logs,rows,matrix)
    metadata=dict(detail,n_log_count_rows=221,n_mechanism_columns=14,n_visibility_columns=2,
                  rank_certificate=len(rows),scope_of_rank=("analytic_map_extension_at_infeasible_point" if bad
                       else "restricted_local_on_box_face" if detail["declared_box_faces"]
                       else "local_at_feasible_interior"),
                  working_precision_bits=precision,arb_inverse_residual_upper=str(bound),
                  exact_rational_inverse_residual_upper=str(independent),
                  strict_supply_minimum_cells=str(strict_supply[2]),
                  strict_supply_minimum_cells_lower=float(strict_supply[2].lower()),
                  violated_donor_constraints=len(bad),
                  minimum_supply_fraction=float(limiting[3].mid()),
                  smallest_surplus_passage=strict_supply[0],smallest_surplus_euler_step=strict_supply[1],
                  minimum_logistic_gap=float(min(crowd)),minimum_pruning_gap=float(min(v for _,v in pruning)),
                  minimum_visible_count=min(certified_counts),operator_branch_masks=",".join(m for m,_ in pruning),
                  native_count_max_rel_error=float(forward_error),native_ols_max_abs_error=float(slope_error),
                  weight_map_max_abs_error=float(native_ols_error),native_gradient_max_scaled_error=gradient_error,
                  old_counts_interval_overlaps=comparison[0],old_count_J_interval_overlaps=comparison[1],
                  old_ols_J_interval_overlaps=comparison[2])
    entries=[dict(case=detail["case"],row_no=a+1,record_id=row_ids[i][1],column=name,
                 scale_exponent=powers[j],**encode(minor[a,j]))
             for a,i in enumerate(rows) for j,name in enumerate(COLUMNS)]
    return metadata,entries,design_check


def run(workers,precision,resume):
    candidates=read(HERE/"point_selection.tsv")
    selected=defaultdict(dict)
    for r in read(HERE/"exact_points.tsv"):
        selected[r["case"]][r["parameter"]]=Fraction(r["exact_value"])
    args=[(r,selected[r["case"]],precision) for r in candidates]
    assert all(set(x[1])==set(COLUMNS) for x in args)
    summary_path=HERE/"certificates.tsv"
    minors_path=HERE/"certified_minors.tsv"
    records=read(summary_path) if resume and summary_path.exists() else []
    mins=read(minors_path) if records else []
    for r in records:
        assert sum(s["case"]==r["case"] for s in mins)==256
        assert int(r["rank_certificate"])==len(COLUMNS)
    completed={r["case"] for r in records}
    failures=[]
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        futures={pool.submit(worker,x):x[0]["case"] for x in args if x[0]["case"] not in completed}
        for f in as_completed(futures):
            try:
                r,mat,check=f.result()
            except Exception as e:
                failures.append((futures[f],e))
                print(f"Verification failed at {futures[f]}: {e!r}",flush=True)
                continue
            records.append(r)
            mins.extend(mat)
            write(summary_path,records)
            write(minors_path,mins)
            write(HERE/"observation_audit.tsv",[check])
            print(r["case"],r["rank_certificate"],r["scope_of_rank"],
                  r["strict_supply_minimum_cells_lower"],flush=True)
    if failures:
        raise RuntimeError(f"Verification failed at {[name for name,_ in failures]}") from failures[0][1]


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--workers",default=7,type=int)
    parser.add_argument("--bits",default=384,type=int)
    parser.add_argument("--resume",action="store_true")
    a=parser.parse_args()
    run(a.workers,a.bits,a.resume)
