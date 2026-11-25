"""Small fixed-limit test of J parametric bootstrap generation and paired refits."""

import argparse
from contextlib import contextmanager
import hashlib
import math
from pathlib import Path
from time import monotonic

from scan import PILOT, Run, apply_fixed, cpu, independent_score, known, read, write
from runtime import DEATH, RHO
from profile_model import EPS, independently
from scipy.linalg import solve_triangular
from scipy.special import expit, gammaln, xlogy
from scipy.stats import betabinom
import numpy as np

REPLICATES = 3
NEAR_EQUAL = 1e-6


def conditional_death(model, channels, p, y):
    live, hypoxia, misseg = channels[:, model.classifier_slots]
    a, loading, attribution, _ = p[DEATH:DEATH+4]
    total_log = np.log(np.maximum(channels[:,model.indices].sum(axis=0)*1e12,1.))
    residual = (y-total_log)[model.death_log_index] * model.death_shared
    log_odds = a + loading*residual + np.log(hypoxia+attribution*misseg+EPS)
    log_odds -= np.log(live+(1.-attribution)*misseg+EPS)
    return EPS +(1.-2.*EPS)*expit(log_odds)


def full_loglikelihood(model,p):
    """Independent actual log probability, including Gaussian and PMF constants."""
    mean = independently(model,p)
    chol = np.linalg.cholesky(model.noise_covariance(p[RHO]))
    z = solve_triangular(chol,model.y-mean["predicted_logs"],lower=True)
    variance = float(z @ z)/model.n
    count = -.5*model.n*(1.+np.log(2.*np.pi*variance))-np.log(np.diag(chol)).sum()
    pi = conditional_death(model, mean["channels"],p,model.y)
    assert np.max(abs(pi-mean["fraction"]))<=1e-12
    kappa = 10.**p[DEATH+3]
    death = betabinom.logpmf(model.dead,model.eligible,kappa*pi,kappa*(1.-pi)).sum()
    kary = np.sum(xlogy(model.kary,mean["probability"])) + np.sum(gammaln(model.kary.sum(axis=1)+1.))
    kary -= np.sum(gammaln(model.kary+1.))
    return float(count+death+kary)


def verify_likelihood_scale(model,x0,x1,conf):
    p0,p1 = (conf.unmap(x)[0] for x in (x0,x1))
    j0,j1 = (independent_score(model,p)["loss"] for p in (p0,p1))
    delta = -2.*(full_loglikelihood(model,p1)-full_loglikelihood(model,p0))
    assert abs(delta-(j1-j0))<2e-7, (delta,j1-j0)
    return delta, j1-j0


def simulate(model,p,seed):
    rng = np.random.default_rng(seed)
    frozen = independently(model,p)
    chol = np.linalg.cholesky(model.noise_covariance(p[RHO]))
    innovation = rng.standard_normal(model.n)
    noise = frozen["sigma"]*(chol @ innovation)
    y = np.ascontiguousarray(frozen["predicted_logs"]+noise,dtype=np.float64)
    assert np.allclose(y-frozen["predicted_logs"],noise,atol=4e-15,rtol=0.)
    pi = conditional_death(model,frozen["channels"],p,y)
    kappa = 10.**p[DEATH+3]
    n = model.eligible.astype(np.int64)
    assert np.array_equal(n,model.eligible) and np.all(pi>0) and np.all(pi<1)
    propensity = rng.beta(kappa*pi,kappa*(1.-pi))
    dead = np.ascontiguousarray(rng.binomial(n,propensity),dtype=np.float64)
    cells = model.kary.sum(axis=1).astype(np.int64)
    kary = np.ascontiguousarray([rng.multinomial(m,prob)
                                 for m,prob in zip(cells,frozen["probability"])],dtype=np.float64)
    assert np.all((dead>=0.)&(dead<=n))
    assert np.array_equal(kary.sum(axis=1),cells)
    assert len(y)==225 and len(dead)==195 and kary.shape==(6,133) and kary.sum()==104
    return (y,dead,kary),pi


@contextmanager
def bound_observations(model, observations):
    names = ("y","dead","kary")
    old = [getattr(model,n) for n in names]
    def assign(arrays):
        for name,value in zip(names,arrays):
            setattr(model,name,value)
        for cache in ("kary_const","death_const"):
            model.__dict__.pop(cache,None)
    assign(observations)
    try:
        yield
    finally:
        assign(old)


def choose_anchor(run,name):
    profile = run.out/"profiles"/name
    verified = [r for r in read(profile/"grid.tsv") if r["status"]=="verified"]
    floor = run.conf.natural(run.baseline_x(),name)
    if name in ("p_misseg","buffer_beta"):
        selected = [r for r in verified if float(r["value"])==0.]
    elif name=="k_o_mis":
        selected = sorted(verified,key=lambda r:float(r["normalized"]))[:1]
    else:
        separated = [r for r in verified if 1e-5<float(r["normalized"])<1-1e-5 and
                     abs(run.conf.target_coordinate(name,floor)-float(r["normalized"]))>.01]
        selected = sorted(separated,key=lambda r:abs(float(r["rise"])))[:1]
    if not selected:
        raise RuntimeError(f"No qualified pilot generating anchor is available for {name}.")
    row,=selected
    fit,=[r for r in read(profile/"attempts.tsv") if r["attempt"]==row["attempt"]]
    return float(row["value"]),run.conf.vector(fit),row


def fit_arm(run,folder,tag,name,value,bank,reserve):
    used=[]
    for i in range(4 if reserve else 3):
        preferred=bank
        if i==3:
            _,floor,agree=run.verdict(read(folder/"attempts.tsv"))
            candidates=[r for r in [floor,*agree] if r is not None]
            preferred=[]
            for candidate in candidates:
                x=run.conf.vector(candidate)
                j=run.conf.idx["rho"]
                x[j]=np.clip(x[j]+(.001 if x[j]<.9 else -.001),0.,1.)
                preferred.append(("Reserve near the simulated feasible endpoint",x))
            preferred.extend(bank)
        r=run.execute(folder,f"{tag}_{i+1}",name,value,preferred,used,i)
        if r["start_"+run.conf.names[0]]!="":
            used.append(run.start(r))
        verdict,floor,agree=run.verdict(read(folder/"attempts.tsv"))
        if verdict=="verified":
            return verdict,floor,agree[0]
    return verdict,floor,agree[0] if agree else None


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("out",type=Path)
    parser.add_argument("parameter",nargs="?",choices=PILOT)
    parser.add_argument("--reserve",action="store_true")
    args=parser.parse_args()
    run=Run(args.out)
    base=run.baseline_x()
    for target in (args.parameter,) if args.parameter else PILOT:
        value,null,receipt=choose_anchor(run,target)
        model,conf=run.model,run.conf
        anchor=conf.unmap(null)[0]
        j_true=model.score(anchor)["loss"]
        direct,discrep=verify_likelihood_scale(model,base,null,conf)
        bank=known(conf)
        top=run.out/"bootstrap_pilot"/target
        write(top/"generator.tsv",[dict(parameter=target,value=value,attempt=receipt["attempt"],
            objective=j_true,full_loglikelihood_difference=direct,score_difference=discrep,
            log_count_scale=independently(model,anchor)["sigma"],replicates=REPLICATES,
            seed_policy="Prefix 20261123 and parameter index and replicate")])
        paired=[]
        for j in range(REPLICATES):
            begin,wall=cpu(),monotonic()
            seed=202611230000+100*PILOT.index(target)+j+1
            obs,pi=simulate(model,anchor,seed)
            repeat,_=simulate(model,anchor,seed)
            assert all(np.array_equal(a,b) for a,b in zip(obs,repeat))
            observed_digest=hashlib.sha256(b"".join(x.tobytes() for x in obs)).hexdigest()
            directory=top/f"replicate_{j+1}"
            sources=(("log_count", model.indices),("dead_objects",model.classifier_slots),
                     ("karyotype_cells",model.kary_rows))
            for (name,indices),array in zip(sources,obs):
                records=np.asarray(array).reshape(len(indices), -1)
                write(directory/f"{name}.tsv",[dict(event=model.slot_ids[slot],row=i,
                      autosome_bin=k+22 if name=="karyotype_cells" else "",value=v)
                      for i,slot in enumerate(indices) for k,v in enumerate(records[i])])
            with bound_observations(model,obs):
                check=independent_score(model,anchor)
                assert np.max(abs(check["fraction"]-pi))<=1e-12
                score_delta, independent_delta=verify_likelihood_scale(model,base,null,conf)
                ls,lf,lc=fit_arm(run,directory/"restricted",f"{target}_boot_{j+1}_r",target,value,
                                  [("generating constrained endpoint",null),("verified measured-data baseline",base)]+bank,args.reserve)
                need = [("simulated restricted minimum",conf.vector(lf))] if lf else []
                us,uf,uc=fit_arm(run,directory/"unrestricted",f"{target}_boot_{j+1}_u",None,None,
                                  need+[("verified measured-data baseline",base),("generating endpoint",null)]+bank,args.reserve)
                lr=float(lc["loss"])-float(uc["loss"]) if lc and uc else math.nan
                nested=(lf is not None and uf is not None and float(uf["loss"])<=float(lf["loss"])+NEAR_EQUAL)
                status="verified" if ls==us=="verified" and nested and lr>=-NEAR_EQUAL else "unresolved"
                rr=read(directory/"restricted"/"attempts.tsv")
                ur=read(directory/"unrestricted"/"attempts.tsv")
                pair_seconds=sum(float(r["cpu_seconds"]) for r in rr+ur)
                actual=cpu()-begin
                paired.append(dict(parameter=target,value=value,replicate=j+1,seed=seed,
                                  observations_sha256=observed_digest,conditional_death_replay=True,
                                  likelihood_scale_agreement=abs(score_delta-independent_delta),
                                  restricted_status=ls,unrestricted_status=us,
                                  restricted_attempt=lc["attempt"] if lc else "",
                                  unrestricted_attempt=uc["attempt"] if uc else "",
                                  restricted_floor=lf["loss"] if lf else "",
                                  unrestricted_floor=uf["loss"] if uf else "",
                                  signed_difference=lr if math.isfinite(lr) else "",
                                  statistic=max(0.,lr) if status=="verified" else "",
                                  nested_check=nested,pair_status=status,
                                  optimization_attempts=sum(r["origin"]=="fit" for r in rr+ur),
                                  evaluations=sum(int(r["calls"]) for r in rr+ur if r["calls"]),
                                  fit_cpu_seconds=pair_seconds,total_cpu_seconds=max(pair_seconds,actual),
                                  wall_seconds=monotonic()-wall))
            assert abs(model.score(anchor)["loss"]-j_true)<1e-7
            write(top/"pairs.tsv",paired)
            print(f"{target} replicate {j+1}: {status}, {pair_seconds:.2f} CPU s.",flush=True)


if __name__=="__main__":
    main()
