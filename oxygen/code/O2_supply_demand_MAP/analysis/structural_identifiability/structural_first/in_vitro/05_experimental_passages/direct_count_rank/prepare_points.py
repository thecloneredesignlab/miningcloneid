"""Select exact parameter points from frozen fits and recorded profile solutions."""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
INHERITANCE = HERE.with_name("20261118_culture_experimental_inheritance_01")
PRODUCTION = HERE.with_name("20261117_culture_actual_schedule_step02")
COUNTS = RESULTS / "synthetic_identifiability/20261118_experimental_counts_practical_01"
sys.path.insert(0, str(COUNTS))
from count_model import MAPPING, OLS, SNAPSHOT, model
from practical import Configuration, MECHANISMS, read, row_x, write

PASSED = ("lam_max", "buffer_smax", "buffer_beta", "buffer_n_exp", "n_O")
IMPROVED = ("p_misseg", "k_o_mis", "gamma_growth", "mu_hp", "O2_crit")
TARGETS = PASSED + IMPROVED
VISIBILITY = ("eta_H_2N", "eta_M_2N")
PROFILE_FILES = tuple(COUNTS / f"profiles_{suffix}.tsv"
                      for suffix in ("direct", "supplement", "edges", "cross", "polish"))
PROFILE_THRESHOLD = 3.841458820694124


def registered():
    source=read(PRODUCTION/"interior_384/exact_parameters.tsv")
    bounds={row["parameter"]:(Fraction(row["lower"]),Fraction(row["upper"]))
            for row in source if row["parameter"] in MECHANISMS}
    bounds.update({name:(Fraction(0),Fraction(1)) for name in VISIBILITY})
    return bounds


def exact_from_saved(row, limits):
    all_values={}
    for key, (lo,hi) in limits.items():
        value=Fraction(row[key])
        # Clamp only a decimal serialization that rounds to a registered bound.
        if value < lo and abs(float(value-lo))<1e-12:
            value=lo
        if value > hi and abs(float(value-hi))<1e-12:
            value=hi
        assert lo<=value<=hi, (key,value,lo,hi)
        all_values[key]=value
    return all_values


def natural_array(c, exact_values):
    result=c.root.copy()
    result[c.active]=[float(exact_values[n]) for n in c.names]
    result[16:18]=result[14:16]
    return result


def from_coordinate(c,x,limits,face=None):
    p=c.unmap(x)[0]
    values=exact_from_saved(dict(zip(c.names,(repr(float(p[i])) for i in c.active))),limits)
    if face:
        values[face]=limits[face][0]
    return values


def search_near_fit(k,c,real_x,limits,face=None):
    mu=c.names.index(face) if face else None
    for step in (1e-8,3e-8,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4,3e-4,1e-3):
        baseline=np.clip(real_x,step,1-step)
        if mu is not None:
            baseline[mu]=0.
        options=[baseline]
        for j in range(len(baseline)):
            if j==mu:
                continue
            for direction in (-1,1):
                b=baseline.copy()
                b[j]=np.clip(b[j]+direction*4*step,step,1-step)
                options.append(b)
        feasible=[]
        for candidate in options:
            values=from_coordinate(c,candidate,limits,face)
            values_p=natural_array(c,values)
            loss,_,margin=k.constrained(values_p,k.y,c,False)
            if margin>5e-9:
                deviation=float(max(abs(c.normal(values_p[c.active])-real_x)))
                feasible.append((deviation,loss,values))
        if feasible:
            return min(feasible,key=lambda row:row[:2])[2]
    raise ValueError(f"No strict supply candidate close to fit for physical face {face}")


def profile_points(k,c,real,limits):
    base=float(real["loss"])
    collected=defaultdict(list)
    for path in PROFILE_FILES:
        for row in read(path):
            if row["feasible"]!="True" or row["parameter"] not in TARGETS:
                continue
            if float(row["min_aggregate_donor_margin_fraction"])<=5e-10:
                continue
            if 221*math.log(float(row["loss"])/base)>PROFILE_THRESHOLD:
                continue
            collected[row["parameter"]].append((path,row))
    for name in TARGETS:
        ref=float(real[f"x_{name}"])
        eligible=collected[name]
        sides=("low","high") if name in IMPROVED else ("farthest",)
        for side in sides:
            pool=eligible
            if side=="low":
                pool=[z for z in eligible if float(z[1]["coordinate"])<ref-1e-4]
                pool=sorted(pool,key=lambda z:(float(z[1]["coordinate"]),float(z[1]["loss"])))
            elif side=="high":
                pool=[z for z in eligible if float(z[1]["coordinate"])>ref+1e-4]
                pool=sorted(pool,key=lambda z:(-float(z[1]["coordinate"]),float(z[1]["loss"])))
            else:
                pool=[z for z in eligible if abs(float(z[1]["coordinate"])-ref)>1e-4]
                pool=sorted(pool,key=lambda z:(-abs(float(z[1]["coordinate"])-ref),float(z[1]["loss"])))
            for path,row in pool:
                values=exact_from_saved(row,limits)
                p=natural_array(c,values)
                loss,_,margin=k.constrained(p,k.y,c,False)
                if margin>5e-10 and 221*math.log(loss/base)<=PROFILE_THRESHOLD:
                    yield f"profile_{name}_{side}", "nominal_profile_support", name, path, values, row["coordinate"]
                    break


def record_provenance():
    registered_stage=HERE.with_name("20261118_culture_experimental_OLS_rank_01")
    calendar=HERE.with_name("20261118_culture_experimental_OLS_step02")
    archive=HERE.with_name("20261109_seed144_wp0")/"source/seed144"
    finite=HERE.parents[2]/"code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/culture_stageB_finite"
    historical=INHERITANCE/"hist/mechanism_witness.tsv"
    paths=(COUNTS/"real_best.tsv",*PROFILE_FILES,COUNTS/"count_inputs.tsv",COUNTS/"parameter_comparison.tsv",
           COUNTS/"count_model.py",COUNTS/"count_kernel.cpp",OLS/"practical.py",OLS/"ols_kernel.cpp",
           MAPPING,SNAPSHOT,calendar/"ols_calendar_map.tsv",calendar/"experimental_ols_reproduction.tsv",
           registered_stage/"certify_experimental.py",archive/"invitro_observed_kary.tsv",
           INHERITANCE/"inherit.py",INHERITANCE/"rank_cases.py",historical,INHERITANCE/"rank_summary.tsv",
           INHERITANCE/"hist/2N_counts_derivatives.tsv",INHERITANCE/"hist/4N_counts_derivatives.tsv",
           INHERITANCE/"hist/2N_counts_replay.tsv",INHERITANCE/"hist/4N_counts_replay.tsv",
           INHERITANCE/"hist/2N_ols_derivatives.tsv",INHERITANCE/"hist/4N_ols_derivatives.tsv",
           PRODUCTION/"interior_384/exact_parameters.tsv",PRODUCTION/"certify.py",
           finite/"model.py",finite/"algebra.py")
    write(HERE/"source_provenance.tsv",[dict(path=str(p),sha256=sha256(p.read_bytes()).hexdigest()) for p in paths])


def run():
    k=model()
    c=Configuration(k)
    limits=registered()
    source=COUNTS/"real_best.tsv"
    real=read(source)[0]
    real_x=row_x(c,real)
    data=[]
    historical=INHERITANCE/"hist/mechanism_witness.tsv"
    assert next(r["exact_value"] for r in read(historical) if r["parameter"]=="k_clear")=="1e-04"
    values={r["parameter"]:Fraction(r["exact_value"]) for r in read(historical) if r["parameter"] in MECHANISMS}
    values.update(eta_H_2N=Fraction(1,3),eta_M_2N=Fraction(2,3))
    data.append(("original_interior","original_ols_certificate","",historical,values,""))
    data.append(("archived_real_fit","saved_numeric_best_fit","",source,exact_from_saved(real,limits),""))
    data.append(("near_fit_interior","strict_surplus_interior_near_fit","",source,search_near_fit(k,c,real_x,limits),""))
    face=search_near_fit(k,c,real_x,limits,"mu_hp")
    data.append(("mu_hp_lower_face","exact_mu_hp_lower_bound_near_fit","mu_hp",source,face,""))
    for step in (Fraction(1,10**8),Fraction(1,10**9),Fraction(1,10**10),Fraction(1,10**11)):
        probe=face | {"mu_hp":face["mu_hp"]+step}
        if k.constrained(natural_array(c,probe),k.y,c,False)[2]>5e-9:
            data.append(("mu_hp_inward","one_sided_mu_hp_direction","mu_hp",source,probe,""))
            break
    else:
        raise ValueError("No finite one-sided mu_hp check passed initial screening")
    data.extend(profile_points(k,c,real,limits))
    informative={a[2] for a in data if a[1]=="nominal_profile_support"}
    assert informative==set(TARGETS),set(TARGETS)-informative
    parameters=[]
    details=[]
    for case,kind,target,path,exact_values,coord in data:
        p=natural_array(c,exact_values)
        loss,_,supply=k.constrained(p,k.y,c,False)
        faces=[f"{n}:{'lower' if exact_values[n]==lo else 'upper'}" for n,(lo,hi) in limits.items()
               if exact_values[n] in (lo,hi)]
        if kind not in ("saved_numeric_best_fit",):
            assert supply>0,(case,supply)
        details.append(dict(case=case,scope=kind,profile_parameter=target,source=str(path),
                            source_pinned_coordinate=coord,loss=loss,
                            count_profile_lr=221*math.log(loss/float(real["loss"])),
                            max_normalized_distance_from_fit=max(abs(c.normal(p[c.active])-real_x)),
                            native_min_supply_fraction=supply,declared_box_faces=",".join(faces)))
        parameters.extend(dict(case=case,parameter=key,exact_value=str(exact_values[key]),
                               lower=str(limits[key][0]),upper=str(limits[key][1])) for key in c.names)
    write(HERE/"point_selection.tsv",details)
    write(HERE/"exact_points.tsv",parameters)
    record_provenance()
    for row in details:
        print(row["case"],row["max_normalized_distance_from_fit"],row["native_min_supply_fraction"],
              row["declared_box_faces"],flush=True)


if __name__=="__main__":
    run()
