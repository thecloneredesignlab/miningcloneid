"""Review only the twelve sealed J bootstrap datasets and their saved fit arms."""

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from math import ceil, isclose
from pathlib import Path
from time import monotonic
from hashlib import sha256

from scan import (ARCHIVES, MATERIAL_SEPARATION, SUSPENDED, Run, apply_fixed, check_sources,
                  cpu, digest, exact_margins, fingerprint, independent_score, normalized_constraint,
                  read, stationarity, write)
from bootstrap_pilot import NEAR_EQUAL, bound_observations, verify_likelihood_scale
from review_support import BudgetStop, refine_observation_block
from scipy.optimize import minimize
import numpy as np

SOURCE = ARCHIVES / "20261123_experimental_joint_all_mechanism_profiles_01"
TARGETS = (
    (202611230203, "restricted"),
    (202611230302, "restricted"),
    (202611230102, "unrestricted"),
    (202611230201, "unrestricted"),
    (202611230303, "restricted"),
    (202611230303, "unrestricted"),
)
CONTROLS = (
    (202611230002, "restricted"),
    (202611230003, "unrestricted"),
    (202611230103, "unrestricted"),
)
ARMS = ("restricted", "unrestricted")
AGREE = .02
FIT_CAP = 31
CPU_CAP = 7200.
REVIEW_NAME = "Frozen twelve-pair refit review"


class BudgetEnd(RuntimeError):
    """The approved local CPU or optimization budget has been used."""


@dataclass(frozen=True)
class Case:
    prior: dict
    observations: tuple
    directory: Path

    @property
    def seed(self):
        return int(self.prior["seed"])

    @property
    def target(self):
        return self.prior["parameter"]

    @property
    def value(self):
        return float(self.prior["value"])

    def fit_target(self, arm):
        return (self.target,self.value) if arm=="restricted" else (None,None)


def load_case(model, prior):
    d = SOURCE/"bootstrap_pilot"/prior["parameter"]/("replicate_"+prior["replicate"])
    generator, = read(d.parent/"generator.tsv")
    assert prior["parameter"]==generator["parameter"] and float(prior["value"])==float(generator["value"])
    arrays=[]
    inputs=(("log_count",model.indices,False),("dead_objects",model.classifier_slots,False),
            ("karyotype_cells",model.kary_rows,True))
    for name,ids,kary in inputs:
        raw=read(d/(name+".tsv"))
        bins=133 if kary else 1
        assert len(raw)==len(ids)*bins
        expected=[(model.slot_ids[slot],str(i),str(22+j) if kary else "") for i,slot in enumerate(ids) for j in range(bins)]
        assert [(r["event"],r["row"],r["autosome_bin"]) for r in raw]==expected, (d,name)
        arr=np.asarray([float(row["value"]) for row in raw],dtype=np.float64)
        assert np.all(np.isfinite(arr))
        arrays.append(arr.reshape(-1,bins) if kary else arr)
    y,dead,kary=arrays
    assert y.shape==(225,) and dead.shape==(195,) and kary.shape==(6,133)
    assert np.all(dead>=0) and np.all(dead<=model.eligible) and np.all(dead==np.floor(dead))
    assert np.all(kary>=0) and np.all(kary==np.floor(kary))
    assert kary.sum()==104 and np.array_equal(kary.sum(axis=1),model.kary.sum(axis=1))
    assert np.array_equal(model.indices[model.death_log_index],model.classifier_slots)
    assert sha256(b"".join(a.tobytes() for a in arrays)).hexdigest()==prior["observations_sha256"], d
    return Case(prior,tuple(arrays),d)


def verify_packet():
    for line in (SOURCE/"SHA256SUMS").read_text().splitlines():
        h,f=line.split(maxsplit=1)
        assert digest(SOURCE/f)==h, f


def active_distance(conf, a, b, fixed):
    active = [i for i in range(len(a)) if i not in fixed]
    for mechanism, carrier in SUSPENDED.items():
        i,j=conf.idx[mechanism],conf.idx[carrier]
        if a[i]==b[i]==0 and j in active:
            active.remove(j)
    return float(np.linalg.norm((a-b)[active]))


def feasible_rows(rows):
    return sorted([r for r in rows if r["loss"]!="" and float(r["minimum_margin"])>=0.],
                  key=lambda r:float(r["loss"]))


def review_arm(run, rows, fixed):
    fits=feasible_rows(rows)
    if not fits:
        return None,[]
    floor=fits[0]
    separate=[]
    for row in fits:
        if row["outcome"]!="qualified" or float(row["loss"])>float(floor["loss"])+AGREE:
            continue
        x=run.start(row)
        if all(active_distance(run.conf,x,run.start(q),fixed)>MATERIAL_SEPARATION for q in separate):
            separate.append(row)
    return floor,separate


class Review:
    def __init__(self,out):
        self.out=out
        self.run=Run(out)
        self.model=self.run.model
        self.conf=self.run.conf
        self.begin=cpu()
        self.deadline=self.begin+CPU_CAP-20.
        self.wall=monotonic()
        self.fits=0
        self.new=[]
        self.cases={}
        self.rows=defaultdict(list)
        self.prior_replay=[]
        self.steps=defaultdict(int)
        self.stop_reason="Scheduled review finished"

    def guard(self):
        if cpu()>=self.deadline:
            raise BudgetEnd("CPU allocation exhausted")

    def preflight(self):
        verify_packet()
        check_sources(self.out)
        prior = read(SOURCE/"bootstrap_pilot_pairs.tsv")
        assert len(prior)==12 and len({r["seed"] for r in prior})==12
        assert Counter(r["parameter"] for r in prior)==Counter(p_misseg=3,buffer_beta=3,k_o_mis=3,n_O=3)
        cases=[load_case(self.model,r) for r in prior]
        self.cases={c.seed:c for c in cases}
        assert sum(c.prior["pair_status"]=="verified" for c in cases)==7
        records=[]
        for c in cases:
            records.append(dict(seed=c.seed,target=c.target,fixed_value=c.value,
               observations_sha256=c.prior["observations_sha256"], counts=225,death=195,
               specimens=6,autosome_bins=133,karyotyped_cells=104,binding="Verified frozen TSV"))
            with bound_observations(self.model,c.observations):
                for arm in ARMS:
                    for row in read(c.directory/arm/"attempts.tsv"):
                        if not row["loss"]:
                            continue
                        name,value=c.fit_target(arm)
                        fixed=self.conf.constraints(name,value)
                        x=self.conf.vector(row)
                        p=self.conf.unmap(x)[0]
                        assert np.array_equal(apply_fixed(x,fixed),x), row["attempt"]
                        baseline=independent_score(self.model,p,row)
                        margin=float(min(exact_margins(self.model,p)))
                        assert isclose(margin,float(row["minimum_margin"]),abs_tol=1e-11,rel_tol=0),row["attempt"]
                        rechecked = dict(row)
                        check = "Independent score and exact feasibility"
                        if row["outcome"]=="qualified":
                            metric=dict(cpu_seconds=0.,wall_seconds=0.,calls=0,iterations=0,
                                solver_rounds=0,solver_message="Replay of sealed endpoint",observation_polish_steps=0)
                            rechecked,_=self.run.accept(x,self.run.start(row),fixed,name,value,
                                row["attempt"],row["source"],"sealed replay",metric,saved=row)
                            assert rechecked["outcome"]==row["outcome"],row["attempt"]
                            check = "Independent score, exact feasibility and local audit"
                        self.rows[c.seed,arm].append(rechecked)
                        self.prior_replay.append(dict(seed=c.seed,arm=arm,attempt=row["attempt"],
                            loss=baseline["loss"],minimum_supply_margin=margin,result=check))
        write(self.out/"frozen_observations.tsv",records)
        write(self.out/"prior_replay.tsv",self.prior_replay)
        controls=[]
        for c in cases:
            if c.prior["pair_status"]!="verified":
                continue
            r,u=[next(r for r in self.rows[c.seed,arm] if r["attempt"]==c.prior[arm+"_attempt"]) for arm in ARMS]
            assert r["outcome"]==u["outcome"]=="qualified"
            with bound_observations(self.model,c.observations):
                delta,score_delta=verify_likelihood_scale(self.model,self.conf.vector(u),self.conf.vector(r),self.conf)
                previous=float(c.prior["statistic"])
                assert abs(max(0.,score_delta)-previous)<=1e-6 and score_delta>=-NEAR_EQUAL
                controls.append(dict(seed=c.seed,prior_statistic=previous,replayed_signed_rise=score_delta,
                                     independently_replayed_rise=delta,**self.check_nesting(c,r)))
        write(self.out/"prior_verified_pair_replay.tsv",controls)
        assert len(controls)==7
        print("Preflight passed for all twelve frozen datasets.",flush=True)

    def status(self,c,arm):
        return review_arm(self.run,self.rows[c.seed,arm],self.conf.constraints(*c.fit_target(arm)))

    def ready(self,c,arm):
        return len(self.status(c,arm)[1])>=2

    def candidates(self,c,arm,kind,step):
        conf,model=self.conf,self.model
        fixed=conf.constraints(*c.fit_target(arm))
        rows=feasible_rows(self.rows[c.seed,arm])
        other_arm = "unrestricted" if arm=="restricted" else "restricted"
        other=feasible_rows(self.rows[c.seed,other_arm])
        af,_=self.status(c,arm)
        best=conf.vector(af)
        starts=[self.run.start(r) for r in self.rows[c.seed,arm]]
        ordered=[]
        if kind=="main" and step==0:
            # A verified interior or boundary counterpart can be close to the newly fixed optimum.
            if arm=="restricted" and c.seed in (202611230203,202611230302):
                ordered.extend(("Projected unrestricted low endpoint",conf.vector(r)) for r in other[:2])
            ordered.extend(("Saved feasible endpoint",conf.vector(r)) for r in rows[:3])
        else:
            ordered.extend(("Lowest qualified endpoint on the other arm",conf.vector(r))
                           for r in other if r["outcome"]=="qualified")
            ordered.extend(("Other-arm feasible endpoint",conf.vector(r)) for r in other[:2])
        if kind=="return":
            ordered=[("Embedded restricted solution",conf.vector(r)) for r in other[:2]]+ordered
        ordered.extend(("Within-arm feasible endpoint",conf.vector(r)) for r in rows[:3])
        # Observation coordinates do not affect living supply. Their perturbations
        # release a distinct initial point even when the two arms share a mode.
        expanded=[]
        for label,x in ordered:
            trial=apply_fixed(x,fixed)
            expanded.append((label,trial))
        for name in ("rho","class_log_odds","shared_count_loading","misseg_dead_share"):
            j=conf.idx[name]
            for size in (.002,.012,.03):
                x=best.copy()
                x[j]=np.clip(x[j]+(size if x[j]<.85 else -size),0.,1.)
                expanded.append((f"Feasible disturbance in {name}",apply_fixed(x,fixed)))
        for label,x in expanded:
            self.guard()
            if min(exact_margins(model,conf.unmap(x)[0]))<0:
                continue
            if all(active_distance(conf,x,v,fixed)>MATERIAL_SEPARATION for v in starts):
                return label,x
        return None

    def scale_at(self, start, active):
        """Estimate diagonal curvature on unchanged native transition support."""
        model,conf=self.model,self.conf
        p,j=conf.unmap(start)
        base=fingerprint(model,p)
        deriv=conf.gradient(model.score(p,True)["gradient"],j)
        scales=np.ones(len(active))
        for k,i in enumerate(active):
            self.guard()
            q=start.copy()
            h=2e-5 if q[i]<.99 else -2e-5
            q[i]+=h
            pp,jj=conf.unmap(q)
            if fingerprint(model,pp)==base:
                g=conf.gradient(model.score(pp,True)["gradient"],jj)
                curvature=abs((g[i]-deriv[i])/h)
                scales[k]=np.clip(np.sqrt(80./max(20.,curvature)),.035,1.2)
        return scales

    def optimize(self, c, arm, kind, step):
        self.guard()
        if self.fits>=FIT_CAP:
            raise BudgetEnd("Optimization allocation exhausted")
        begun,wall=cpu(),monotonic()
        chosen=self.candidates(c,arm,kind,step)
        if chosen is None:
            write(self.out/"start_failures.tsv",[dict(seed=c.seed,arm=arm,phase=kind,
                     reason="Distinct exactly feasible start search exhausted",cpu_seconds=cpu()-begun)],append=True)
            return False
        source,start=chosen
        model,conf=self.model,self.conf
        target,value=c.fit_target(arm)
        fixed=conf.constraints(target,value)
        p,_=conf.unmap(start)
        m=exact_margins(model,p)
        assert len(m)==110 and min(m)>=0
        active=np.array([i for i in range(len(start)) if i not in fixed])
        scales=self.scale_at(start,active)
        # The reserve works with a gentler scaling and longer iteration allowance.
        if kind=="reserve":
            scales=np.sqrt(scales)
        offset=model.score(p)["loss"]
        best=[offset,start.copy()]
        evaluations=0
        total_constraints=set(np.argsort(m)[:4])
        total_constraints.update(np.flatnonzero(m<.025))
        for saved in feasible_rows(self.rows[c.seed,arm])[:3]:
            total_constraints.update(int(label.split("_",1)[1]) for label in saved["active_constraints"].split()
                                     if label.startswith("donor_"))
        ids=np.array(sorted(total_constraints))
        bounds=[(-start[i]/s,(1.-start[i])/s) for i,s in zip(active,scales)]
        divisor=12. if kind=="main" else 6.
        constraint_scale=4.
        def from_scaled(t):
            x=start.copy()
            x[active]=np.clip(start[active]+scales*t,0.,1.)
            return x

        def objective(t):
            nonlocal evaluations
            self.guard()
            x=from_scaled(t)
            a,b=conf.unmap(x)
            score=model.score(a,True)
            evaluations+=1
            if score["loss"]<best[0] and min(exact_margins(model,a))>=0.:
                best[:]=[float(score["loss"]),x.copy()]
            return (score["loss"]-offset)/divisor, conf.gradient(score["gradient"],b)[active]*scales/divisor

        def values(t):
            self.guard()
            a=conf.unmap(from_scaled(t))[0]
            return constraint_scale*(exact_margins(model,a)[ids]-1e-9)

        def jacobian(t):
            self.guard()
            a,b=conf.unmap(from_scaled(t))
            return normalized_constraint(model,conf,a,b,active,ids,[])[1]*scales*constraint_scale

        self.fits+=1
        self.steps[c.seed,arm,kind]+=1
        label=f"{c.seed}_{arm[0]}_{kind}_{step+1}"
        interrupted=False
        try:
            result=minimize(objective,np.zeros(len(active)),jac=True,method="SLSQP",bounds=bounds,
                constraints=[dict(type="ineq",fun=values,jac=jacobian)],
                options=dict(maxiter=650 if kind=="reserve" else 420,ftol=3e-11))
            candidate=from_scaled(result.x)
            val=model.score(conf.unmap(candidate)[0])
            if min(exact_margins(model,conf.unmap(candidate)[0]))>=0 and val["loss"]<=best[0]+1e-8:
                best[:]=[val["loss"],candidate]
            x,steps=refine_observation_block(model,conf,best[1],fixed,self.deadline)
            msg=str(result.message)
            iterations=result.nit
        except (BudgetEnd, BudgetStop):
            interrupted=True
            x,steps=best[1],0
            msg="CPU allocation exhausted during optimization"
            iterations=""
        except (ArithmeticError, np.linalg.LinAlgError) as error:
            # Retain the best exactly feasible iterate for the current frozen sample.
            x,steps=best[1],0
            msg=f"Numerical solver error: {error}"
            iterations=""
        row,audit=self.run.accept(x,start,fixed,target,value,label,source,"refit",
            dict(cpu_seconds=0,wall_seconds=0,calls=evaluations,iterations=iterations,solver_rounds=1,
                 solver_message=msg,observation_polish_steps=steps))
        row.update(cpu_seconds=cpu()-begun,wall_seconds=monotonic()-wall)
        self.rows[c.seed,arm].append(row)
        write(self.out/"refit_endpoints.tsv",[dict(seed=c.seed,arm=arm,phase=kind,**row)],append=True)
        write(self.out/"direction_checks.tsv",[dict(seed=c.seed,arm=arm,**a) for a in audit],append=True)
        self.new.append(dict(seed=c.seed,arm=arm,phase=kind,attempt=label,optimizer_invocations=1,
                  cpu_seconds=row["cpu_seconds"],evaluations=evaluations,numerical_outcome=row["outcome"]))
        write(self.out/"refit_budget_ledger.tsv",self.new)
        print(f"{label}: {row['outcome']}, J={row['loss']:.9f}, KKT={row['stationarity']:.4g}, "
              f"CPU={row['cpu_seconds']:.2f}s, {msg}",flush=True)
        if interrupted:
            raise BudgetEnd(msg)
        return True

    def recover(self):
        for seed,arm in TARGETS+CONTROLS:
            c=self.cases[seed]
            with bound_observations(self.model,c.observations):
                for trial in range(2):
                    if self.ready(c,arm):
                        break
                    self.optimize(c,arm,"main",trial)
        for seed,arm in TARGETS+CONTROLS:
            c=self.cases[seed]
            with bound_observations(self.model,c.observations):
                if not self.ready(c,arm):
                    self.optimize(c,arm,"reserve",0)
        for seed in sorted({seed for seed,_ in TARGETS}):
            c=self.cases[seed]
            with bound_observations(self.model,c.observations):
                r,rs=self.status(c,"restricted")
                u,us=self.status(c,"unrestricted")
                if r and u and (float(u["loss"])>float(r["loss"])+NEAR_EQUAL or
                                (rs and us and float(us[0]["loss"])>float(rs[0]["loss"])+NEAR_EQUAL)):
                    if sum(n for (_,_,kind),n in self.steps.items() if kind=="return")<4:
                        self.optimize(c,"unrestricted","return",0)

    def check_nesting(self, c, restricted):
        x=self.conf.vector(restricted)
        assert np.array_equal(x,apply_fixed(x,self.conf.constraints(c.target,c.value)))
        p_r=self.conf.unmap(x)[0]
        p_u=self.conf.unmap(apply_fixed(x,self.conf.constraints(None,None)))[0]
        r,u=(independent_score(self.model,p) for p in (p_r,p_u))
        mr,mu=(exact_margins(self.model,p) for p in (p_r,p_u))
        d=max(float(np.max(abs(r[k]-u[k]))) for k in ("predicted_logs","channels","fraction","probability"))
        delta=abs(r["loss"]-u["loss"])
        dm=float(np.max(abs(mr-mu)))
        assert d<=1e-10 and dm<=1e-12 and delta<=1e-7 and min(mu)>=0.
        return dict(embedded_restricted_loss=r["loss"],same_vector_unrestricted_loss=u["loss"],
                    absolute_score_error=delta,maximum_prediction_error=d,maximum_margin_error=dm)

    def finish(self):
        pairs=[]
        nesting=[]
        for c in self.cases.values():
            with bound_observations(self.model,c.observations):
                r,rs=self.status(c,"restricted")
                u,us=self.status(c,"unrestricted")
                assert r and u
                picked_r,picked_u=(rs[0] if rs else r),(us[0] if us else u)
                chosen_delta=float(picked_r["loss"])-float(picked_u["loss"])
                floor_delta=float(r["loss"])-float(u["loss"])
                ell,score=verify_likelihood_scale(self.model,self.conf.vector(picked_u),self.conf.vector(picked_r),self.conf)
                error=abs(ell-score)
                nested=chosen_delta>=-NEAR_EQUAL and floor_delta>=-NEAR_EQUAL
                accept=len(rs)>=2 and len(us)>=2 and nested
                record=dict(seed=c.seed,target=c.target,fixed_value=c.value,prior_status=c.prior["pair_status"],
                    restricted_attempt=picked_r["attempt"],unrestricted_attempt=picked_u["attempt"],
                    restricted_lowest_feasible_score=r["loss"],unrestricted_lowest_feasible_score=u["loss"],
                    restricted_accepted_score=picked_r["loss"],unrestricted_accepted_score=picked_u["loss"],
                    restricted_KKT=picked_r["stationarity"],unrestricted_KKT=picked_u["stationarity"],
                    restricted_endpoint_qualified=bool(rs),unrestricted_endpoint_qualified=bool(us),
                    restricted_qualified_starts=len(rs),unrestricted_qualified_starts=len(us),
                    signed_score_rise=chosen_delta,signed_floor_rise=floor_delta,
                    calibrated_statistic_input=max(0.,chosen_delta) if accept else "",
                    independent_likelihood_error=error,nesting_pass=nested,
                    final_status="verified" if accept else "unresolved")
                pairs.append(record)
                for role,pair in (("Lowest feasible",r),("Selected qualified" if rs else "Feasible candidate",picked_r)):
                    nesting.append(dict(seed=c.seed,endpoint=role,attempt=pair["attempt"],
                                        **self.check_nesting(c,pair)))
        write(self.out/"pair_review.tsv",pairs)
        write(self.out/"nesting_receipts.tsv",nesting)
        return pairs


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("out",type=Path)
    args=cli.parse_args()
    out=args.out.resolve()
    assert out.is_relative_to(ARCHIVES.resolve()) and not out.is_relative_to(SOURCE.resolve())
    out.mkdir(exist_ok=False)
    review=Review(out)
    review.preflight()
    try:
        review.recover()
    except BudgetEnd as error:
        review.stop_reason=str(error)
    result=review.finish()
    incremental=cpu()-review.begin
    found=sum(r["final_status"]=="verified" for r in result)
    prior_verified_stable=[r for r in result if r["prior_status"]=="verified"]
    gate = found==12 and len(prior_verified_stable)==7 and all(r["final_status"]=="verified" for r in prior_verified_stable)
    assert review.fits<=FIT_CAP and incremental<=CPU_CAP
    assert all(v<=1 for (_,_,kind),v in review.steps.items() if kind=="reserve")
    assert sum(n for (_,_,kind),n in review.steps.items() if kind=="return")<=4
    times=defaultdict(float)
    for r in review.new:
        times[r["seed"]]+=float(r["cpu_seconds"])
    old={c.seed:float(c.prior["total_cpu_seconds"]) for c in review.cases.values()}
    pooled=np.array([t+times[k] for k,t in old.items()])
    rate=1.25*np.quantile(pooled,.75)
    estimate=ceil(rate*199*4/3600.)
    write(out/"batch_gate.tsv",[dict(verified_pairs=found,frozen_datasets=12,numeric_gate_pass=gate,
        separate_batch_approval_required=True,observed_data_calibration_started=False,calibration_batch_started=False,
        median_recovery_inclusive_cpu_seconds_per_pair=np.median(pooled),
        planning_cpu_seconds_per_pair=rate,projected_core_hours_4_by_199=estimate)])
    write(out/"resources.tsv",[dict(numerical_threads=1,optimizer_invocations=review.fits,
        reserve_invocations=sum(n for (_,_,kind),n in review.steps.items() if kind=="reserve"),
        cross_arm_return_invocations=sum(n for (_,_,kind),n in review.steps.items() if kind=="return"),
        cpu_seconds=incremental,core_hours=incremental/3600,elapsed_seconds=monotonic()-review.wall,
        allowed_cpu_seconds=CPU_CAP,allowed_invocations=FIT_CAP,stop_reason=review.stop_reason)])
    print(f"Verified {found}/12; {review.fits} new optimizations; {incremental:.2f} CPU s; "
          f"provisional batch expense {estimate} core-hours.",flush=True)


if __name__=="__main__":
    main()
