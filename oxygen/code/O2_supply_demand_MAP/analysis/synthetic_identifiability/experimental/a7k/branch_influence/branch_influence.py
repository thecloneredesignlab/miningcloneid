"""Likelihood profiles after retaining selected A7K karyotype observations."""
import argparse
from concurrent.futures import as_completed
from functools import lru_cache
from pathlib import Path
import sys

import numpy as np
from scipy.special import xlogy

OUT = Path(__file__).resolve().parent
ARCH = OUT.with_name("20261119_A7K_inheritance_sensitivity_01")
sys.path.insert(0, str(ARCH))
from tilt_model import (BASE, EPS, PREV, THRESHOLD, Config, TiltModel, case, checked_fit,
                        pool, read, receipt, stable, take, write)

MODES = ("strict", "group_050")
NAMES = ("2N_O1", "2N_O2", "4N_O1", "4N_O2")


@lru_cache
def masked(mode, mask):
    k = TiltModel(mode)
    c = Config(k)
    assert tuple(c.kary_labels) == NAMES
    keep = np.r_[True, True, [(mask & (1 << j)) != 0 for j in range(4)]].astype(bool)
    k.kary = k.kary[keep].copy()
    k.kary_rows = k.kary_rows[keep].copy()
    k.kary_size = k.kary.size
    fixed = {}
    for j, label in enumerate(("2N", "4N")):
        if k.cap and not mask & (3 << (2 * j)):
            fixed.update({c.names.index(f"tilt_{q}_{label}"): .5 for q in ("mean", "shape")})
    return k, c, k.target(), fixed


def tag(mode, mask):
    return f"{mode}_{mask:02d}"


def solution(mode, mask, start, origin, beta=None, polish=False):
    k, c, y, fixed = masked(mode, mask)
    fixed = dict(fixed)
    if beta is not None:
        fixed[c.names.index("buffer_beta")] = beta
    fit = checked_fit(k, c, y, start, fixed=fixed,
                      limit=150 if polish else 85, divisor=17. if polish else 38.)
    return receipt(k, c, y, fit, mode=mode, mask=mask,
                   coordinate=-1. if beta is None else beta, origin=origin)


def archived_rows(mode):
    f = PREV/"profile_points_stock_a7k.tsv" if mode == "strict" else ARCH/f"profile_points_{mode}.tsv"
    raw = [r for r in read(f) if r["parameter"] == "buffer_beta" and
           float(r.get("all_constraint_margin", r["donor_margin"])) >= 0.]
    rows = {}
    for r in sorted(raw, key=lambda r: float(r["loss"])):
        rows.setdefault(float(r["coordinate"]), r)
    b = read(PREV/"best_stock_baseline.tsv")[0]
    a = read(ARCH/f"best_{mode}.tsv")[0]
    return b, a, rows


def fit_mask(mode, mask):
    k, c, y, _ = masked(mode, mask)
    base, full, old = archived_rows(mode)
    q = min(old.values(), key=lambda r: abs(float(r["coordinate"]) - .60))
    qhi = min(old.values(), key=lambda r: abs(float(r["coordinate"]) - .99))
    bank = [base, full, q, qhi]
    rs = [solution(mode, mask, c.anchor(row), f"free_{i}") for i, row in enumerate(bank)]
    best = take(rs)
    if not stable(best):
        rs.append(solution(mode, mask, c.anchor(best), "free_polish", polish=True))
        best = take(rs)
    ends = [solution(mode, mask, c.anchor(row), f"beta10_{i}", 1.)
            for i, row in enumerate([base, qhi, best])]
    at10 = take(ends)
    if not stable(at10):
        ends.append(solution(mode, mask, c.anchor(at10), "beta10_polish", 1., True))
        at10 = take(ends)
    if float(at10["loss"]) < float(best["loss"]) - 1e-5:
        rs.append(solution(mode, mask, c.anchor(at10), "free_from_beta10", polish=True))
        best = take(rs)
    prefix = tag(mode, mask)
    write(OUT/f"fit_points_{prefix}.tsv", rs+ends)
    write(OUT/f"fit_ref_{prefix}.tsv", [best, at10])
    return mode, mask, best["buffer_beta"], at10["loss"]-best["loss"], best["stationarity"], at10["stationarity"]


def profile(mode, mask):
    k, c, y, _ = masked(mode, mask)
    base, full, old = archived_rows(mode)
    prefix = tag(mode, mask)
    references = read(OUT/f"fit_ref_{prefix}.tsv")
    best, at10 = references
    j = c.names.index("buffer_beta")
    points = {1.: at10}
    raw = [dict(at10)]

    def evaluate(v, warm, alternative):
        near = min(old.values(), key=lambda r: abs(float(r["coordinate"])-v))
        starts = [warm]
        for anchor in ([near, base if v >= .6 else full] if alternative else []):
            if min(np.linalg.norm(np.delete(c.anchor(z)-c.anchor(anchor),j)) for z in starts) > .035:
                starts.append(anchor)
        rows = [solution(mode, mask, c.anchor(r), f"profile_{i}", v)
                for i, r in enumerate(starts)]
        chosen = take(rows)
        if not stable(chosen) and abs(float(chosen["loss"])-float(best["loss"])-THRESHOLD) < 12.:
            rows.append(solution(mode, mask, c.anchor(chosen), "profile_polish", v, True))
            chosen = take(rows)
        raw.extend(rows)
        points[v] = chosen
        return chosen

    anchor = float(best["x_buffer_beta"])
    grid = {0.,.08,.15,.23,.33,.46,.60,.72,.78538713,.88,.96}
    for sign in (-1,1):
        warm = best
        for v in sorted((v for v in grid if sign*(v-anchor) > 1e-7), reverse=sign<0):
            warm = evaluate(v,warm,True)
    candidate = take([best,*points.values()])
    if float(candidate["loss"]) < float(best["loss"]) - .005:
        new = solution(mode,mask,c.anchor(candidate),"profile_restart",polish=True)
        if float(new["loss"]) < float(best["loss"]):
            best = new
    anchor = float(best["x_buffer_beta"])
    winner = float(best["loss"])
    outside_type = []
    ends = []
    for sign in (-1,1):
        boundary = 0. if sign < 0 else 1.
        if abs(anchor-boundary) < 1e-8:
            ends.append(boundary)
            outside_type.append("box")
            continue
        for iteration in range(10):
            eligible = {v:r for v,r in points.items() if sign*(v-anchor)>1e-8}
            inside = [v for v,r in eligible.items() if float(r["all_constraint_margin"])>=0 and
                      float(r["loss"])-winner<=THRESHOLD]
            edge = (max if sign>0 else min)([anchor,*inside])
            if edge == boundary:
                ends.append(edge)
                outside_type.append("box")
                break
            outer = [v for v in eligible if sign*(v-edge)>0.]
            if not outer:
                ends.append(np.nan)
                outside_type.append("unresolved")
                break
            nearest = min(outer,key=lambda v: abs(v-edge))
            inner_row = best if edge == anchor else points[edge]
            outer_row = points[nearest]
            li, lo = (float(r["loss"])-winner for r in (inner_row,outer_row))
            if abs(edge-nearest)<.003 and stable(outer_row) and stable(inner_row):
                ends.append(edge+(nearest-edge)*(THRESHOLD-li)/(lo-li))
                outside_type.append("crossing")
                break
            if iteration==9:
                ends.append(edge+(nearest-edge)*(THRESHOLD-li)/(lo-li) if lo>li else np.nan)
                outside_type.append("crossing_recheck" if not (stable(outer_row) and stable(inner_row)) else "crossing")
                break
            if abs(edge-nearest)<.003:
                for cur in (edge,nearest):
                    if cur!=anchor and not stable(points[cur]):
                        new = solution(mode, mask, c.anchor(points[cur]),"edge_check",cur,True)
                        raw.append(new)
                        points[cur] = take([points[cur],new])
            evaluate((edge+nearest)/2.,inner_row,False)
    pvalue = lambda u: c.unmap(np.where(np.arange(len(c.names))==j,u,c.anchor(best)))[0][c.active[j]]
    at10_lr = float(at10["loss"])-winner
    for r in raw:
        r["lr"]=float(r["loss"])-winner
    write(OUT/f"profile_points_{prefix}.tsv",raw)
    write(OUT/f"profile_ref_{prefix}.tsv",[best,at10])
    result=dict(mode=mode,mask=mask,branches=",".join(n for i,n in enumerate(NAMES) if mask & 1<<i),
                n_kary_cells=int(k.kary[2:].sum()), beta=float(best["buffer_beta"]),
                lower=pvalue(ends[0]),upper=pvalue(ends[1]),
                lower_type=outside_type[0],upper_type=outside_type[1],
                x_lower=ends[0],x_upper=ends[1],deviance=winner,
                beta10_deviance=float(at10["loss"]),beta10_lr=at10_lr,beta10_supported=at10_lr<=THRESHOLD,
                free_stationarity=float(best["stationarity"]),beta10_stationarity=float(at10["stationarity"]))
    write(OUT/f"interval_{prefix}.tsv",[result])
    return result


def candidates(stage, mode, mask):
    if stage == "fit":
        return fit_mask(mode,mask)
    return profile(mode,mask)


def audit():
    out=[]
    for mode in MODES:
        ref, comparison, old = archived_rows(mode)
        km,cm,ym,_ = masked(mode,0)
        k,c,y = case(mode)
        assert km.kary.tolist()==k.kary[:2].tolist()
        assert np.array_equal(ym[:k.n],y[:k.n]) and km.n==221
        for mask in (0,15,1,2,4,8,14,13,11,7):
            sub,sc,sy,_ = masked(mode,mask)
            p,_ = sc.unmap(sc.anchor(comparison))
            ar, f, _, _, b = sub.parts(p,sy)
            allp,_=c.unmap(c.anchor(comparison))
            ra,ff, _,_,bf=k.parts(allp,y)
            assert abs(ar[0]-ra[0])<1e-10 and abs(b-bf)<1e-10
            pro, fr = deviance_rows(k,allp,y)
            indexed=[j for j in range(4) if mask & 1<<j]
            assert abs(f-b-sum(fr[j+2] for j in indexed))<1e-8
            if mask==15:
                assert abs(f-float(comparison["loss"]))<1e-8
            out.append(dict(mode=mode,mask=mask,rows=int(sub.kary.shape[0]),
                            cells=int(sub.kary[2:].sum()),loss=f,base=b))
        p,_=cm.unmap(cm.anchor(ref))
        assert abs(km.parts(p,ym)[1]-float(ref["loss"]))<1e-8
    write(OUT/"mask_replay.tsv",out)
    print("Reproduced common and full loss and all retained per-branch terms.",flush=True)


def deviance_rows(k,p,y):
    _,states = k.readouts(p)
    pr=(1.-EPS)*states/states.sum(axis=1,keepdims=True)+EPS/133.
    arr=y[k.n:].reshape(k.kary.shape)
    empirical=arr/arr.sum(axis=1,keepdims=True)
    return pr,2.*(xlogy(arr,empirical)-xlogy(arr,pr)).sum(axis=1)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("stage",choices=("audit","fit","profile"))
    p.add_argument("--mode",choices=MODES)
    p.add_argument("--masks",nargs="+",type=int,default=[1,2,4,8,14,13,11,7])
    p.add_argument("--workers",type=int,default=8)
    a=p.parse_args()
    if a.stage=="audit":
        audit()
    else:
        jobs=[(a.stage,m,s) for m in ([a.mode] if a.mode else MODES) for s in a.masks]
        if len(jobs)==1:
            print(candidates(*jobs[0]),flush=True)
        else:
            with pool(a.workers) as ex:
                for f in as_completed(ex.submit(candidates,*u) for u in jobs):
                    print(f.result(),flush=True)
