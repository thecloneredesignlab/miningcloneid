"""Choose registered endpoint probes and describe the sampled component around the fit."""

from math import isclose
import re
from scan import read
from endpoint_zero_review import SOURCE
from second_endpoint_prepare import FIRST

RESOLUTION = .02
PRIORITY = ("p_misseg","buffer_beta","lam_max","buffer_smax","p_mis_base","k_o_mis","n_O",
            "p_wgd","alpha_o2","buffer_n_exp","gamma_growth","mu_hp","gamma_mu","O2_crit")
SECOND = FIRST.with_name("20261124_joint_support_second_release_01")


def registered_value(conf,base,name,coordinate):
    temp = base.copy()
    temp[conf.idx[name]] = coordinate
    return conf.natural(temp,name)


def source_points(conf,previous):
    points = []
    for r in read(SECOND/"tested_values.tsv"):
        name = r["parameter"]
        v = float(r["fixed_value"])
        points.append(dict(parameter=name,value=v,coordinate=conf.target_coordinate(name,v),
                           decision=r["decision"],source=r["evidence"],packet=SECOND.name,
                           reason="Sealed pointwise calibration"))
    points.extend(previous)
    return points


def split_candidate(a,b,occupied):
    # Probe the widest untested gap. Retain both sides of an undecided probe.
    low,high = sorted((a,b))
    marks = [low,*sorted({x for x in occupied if low+1e-9 < x < high-1e-9}),high]
    left,right = max(zip(marks,marks[1:]),key=lambda r:r[1]-r[0])
    if right-left <= RESOLUTION/2 + 1e-12:
        return None
    return (left+right)/2


def evaluate_side(conf,base,reg,side,points,grid):
    name = reg["parameter"]
    center = conf.natural(base,name)
    cx = conf.target_coordinate(name,center)
    lower = side=="Lower"
    bound = float(reg["domain_lower" if lower else "domain_upper"])
    bx = conf.target_coordinate(name,bound)
    sign = -1 if lower else 1
    ordered = sorted((p for p in points if p["parameter"]==name and sign*(float(p["coordinate"])-cx)>1e-9),
                     key=lambda p:abs(float(p["coordinate"])-cx))
    retained = center
    excluded = None
    reentry = []
    unresolved = []
    for p in ordered:
        v = float(p["value"])
        if p["decision"]=="Retained":
            if excluded is None:
                retained = v
            else:
                reentry.append(v)
        elif p["decision"]=="Excluded" and excluded is None:
            excluded = v
        elif p["decision"] not in ("Excluded","Retained"):
            unresolved.append(v)
    rx = conf.target_coordinate(name,retained)
    ex = conf.target_coordinate(name,excluded) if excluded is not None else None
    width = abs(rx-ex) if ex is not None else None
    candidate = None
    occupied = [float(p["coordinate"]) for p in ordered]
    def unused(coord):
        return min([abs(coord-cx)]+[abs(v-coord) for v in occupied])>1e-8
    def in_search(coord):
        return sign*(coord-rx)>1e-9 and (ex is None or sign*(ex-coord)>1e-9)
    if abs(bx-cx)<1e-9:
        status = "Fitted domain boundary"
    elif width is not None and width<=RESOLUTION+1e-12:
        status = "Statistical endpoint bracket"
    elif abs(rx-bx)<1e-9:
        status = "Supported to domain boundary"
    else:
        status = "Unresolved"
        if ex is not None:
            candidate = split_candidate(rx,ex,occupied)
        else:
            # First use the actual archived values behind the rounded search register.
            search = reg["lower_search" if lower else "upper_search"]
            usable = [p for p in grid if p["parameter"]==name and p["status"]=="verified"]
            for text in re.findall(r"(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?",search):
                approx = float(text)
                options = [p for p in usable if isclose(float(p["value"]),approx,rel_tol=.005,abs_tol=2e-7)]
                if options:
                    actual = min(options,key=lambda p:abs(float(p["value"])-approx))
                    loc = float(actual["normalized"])
                    if in_search(loc) and unused(loc):
                        candidate = loc
                        break
            if candidate is None and unused(bx):
                # Prefer a feasible outer point when the domain face already failed the observed grid.
                face = [p for p in grid if p["parameter"]==name and abs(float(p["normalized"])-bx)<1e-8]
                if face and face[0]["status"]!="verified":
                    outer = [float(p["normalized"]) for p in usable if in_search(float(p["normalized"]))
                             and unused(float(p["normalized"]))]
                    if outer:
                        candidate = max(outer,key=lambda v:abs(v-cx))
                if candidate is None:
                    candidate = bx
            if candidate is None and abs(bx-rx)>RESOLUTION:
                candidate = split_candidate(rx,bx,occupied)
        if candidate is not None:
            matching = [p for p in grid if p["parameter"]==name and abs(float(p["normalized"])-candidate)<1e-10]
            value = float(matching[0]["value"]) if matching else registered_value(conf,base,name,candidate)
            candidate = value if unused(conf.target_coordinate(name,value)) else None
    stop = ""
    if status=="Unresolved":
        outer = bx if ex is None else ex
        low,high = sorted((rx,outer))
        undecided = [p for p in ordered if low-1e-9<=float(p["coordinate"])<=high+1e-9
                     and p["decision"] not in ("Retained","Excluded")]
        for p in undecided:
            x = float(p["coordinate"])
            if any(float(t["coordinate"])<x-1e-9 for t in undecided) and any(float(t["coordinate"])>x+1e-9 for t in undecided):
                candidate = None
                stop = "Undecided probe and both flanking probes remain open; precision review required"
                break
        if candidate is None and not stop:
            stop = "No further independent registered probe at the target resolution"
    return dict(parameter=name,side=side,estimate=center,domain_bound=bound,supported=retained,
                excluded="" if excluded is None else excluded,width="" if width is None else width,status=status,
                unresolved_points=" ".join(format(v,".16g") for v in unresolved),
                possible_reentry=" ".join(format(v,".16g") for v in reentry),
                candidate="" if candidate is None else candidate,stop_reason=stop,
                sealed_new_anchors=sum(p.get("generation")=="V6" for p in ordered))


def inventory(run,points):
    base = run.baseline_x()
    reg = read(FIRST/"endpoint_register.tsv")
    grid = read(SOURCE/"profile_grid.tsv")
    assert len(reg)==14 and {r["parameter"] for r in reg}==set(PRIORITY)
    return [evaluate_side(run.conf,base,r,side,points,grid) for r in reg for side in ("Lower","Upper")]


def queue(sides):
    # Round-robin coverage across all open directions; prioritize the registered focal
    # mechanisms within each pass. Keep different sides of the same parameter independent.
    options = sorted((r for r in sides if r["candidate"]!=""),
                     key=lambda r:(r["sealed_new_anchors"],PRIORITY.index(r["parameter"]),r["side"]))
    chosen = []
    for r in options:
        if not any(p["parameter"]==r["parameter"] for p in chosen):
            chosen.append(r)
        if len(chosen)==4:
            break
    return chosen
