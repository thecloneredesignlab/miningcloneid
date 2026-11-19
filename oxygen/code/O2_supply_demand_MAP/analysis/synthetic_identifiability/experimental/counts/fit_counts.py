"""Noiseless recovery and constrained multistart fitting of observed counts."""
import argparse
import concurrent.futures as cf

import numpy as np
from scipy.optimize import least_squares

from count_model import HERE, OLS, model, optimize
from inference import likelihood_ratio, pool, SLACK
from practical import Configuration, read, result_row, row_x, write


def fit_noiseless(job):
    label, truth, start = job
    k=model()
    c=Configuration(k)
    logs,diag = k.log_predict(c.unmap(truth)[0])
    assert diag[0] == 0
    target=logs[k.capture]

    def residual(x):
        logs,diag = k.log_predict(c.unmap(x)[0])
        return logs[k.capture] - target + (15*(1+diag[0]) if diag[0]>0 else 0.)
    opt = least_squares(residual,start,bounds=(0,1),max_nfev=75,x_scale="jac",
                        gtol=1e-12,xtol=1e-12,ftol=1e-12,diff_step=1e-5)
    f,_,margin = k.constrained(c.unmap(opt.x)[0],target,c,False)
    fit = dict(x=opt.x,loss=f,supply=margin,gradient=opt.optimality,calls=opt.nfev,
               iteration=opt.nfev,scipy_success=opt.success,termination=str(opt.message),seconds="")
    return result_row(c,fit,origin=label,max_log_count_error=max(abs(residual(opt.x))),
                      max_mechanistic_error=max(abs(truth[:14]-opt.x[:14])))


def solve_real(job):
    label, x, limit, *scaling = job
    k = model()
    c = Configuration(k)
    fit=optimize(k,c,k.y,x,limit,divisor=scaling[0] if scaling else 10.)
    return result_row(c,fit,origin=label)


def sorted_feasible(rows):
    return sorted((r for r in rows if float(r["min_aggregate_donor_margin_fraction"]) >= -SLACK),
                  key=lambda r:float(r["loss"]))


def record_real(k,c,rows):
    feasible=sorted_feasible(rows)
    assert feasible
    best=feasible[0]
    close=[r for r in feasible if float(r["max_projected_gradient"]) <= .025
           and likelihood_ratio(float(r["loss"]),float(best["loss"]),k.n) < .001]
    if close:
        best=close[0]
    write(HERE/"real_best.tsv",[best])
    logs,_=k.log_predict(c.unmap(row_x(c,best))[0])
    write(HERE/"real_residuals.tsv", [dict(r,predicted_count=float(np.exp(p)),
                                          log_residual=float(np.log(r["corrected_count"])-p))
                                      for r,p in zip(k.observations,logs) if r["fitted_count"]])
    print("Real reference",best["origin"],best["loss"],best["max_projected_gradient"],
          "sigma",np.sqrt(2*float(best["loss"])/k.n),flush=True)
    return best


def real_starts(c):
    rng=np.random.default_rng(153312)
    real = read(OLS/"main_constrained_best.tsv")
    base = row_x(c,real[0])
    previous=read(OLS/"main_constrained_refinement.tsv")+read(OLS/"profile_reference_refits.tsv")
    points=[base,c.x_anchor]
    for r in sorted_feasible(previous):
        a=row_x(c,r)
        if min(np.linalg.norm(a-x) for x in points)>.22:
            points.append(a)
        if len(points)==6:
            break
    return points + [np.clip(base+rng.normal(0,s,len(base)),0,1) for s in (.008,.02,.04,.075)]


def run(mode, workers):
    k=model()
    c=Configuration(k)
    rng=np.random.default_rng(210188)
    if mode=="noiseless":
        truth=np.clip(c.x_anchor,.015,.985)
        starts=[np.clip(c.x_anchor+rng.normal(0,d,len(truth)),.001,.999)
                for d in (.02,.04,.07,.10,.12,.20)]
        worker=fit_noiseless
        jobs=[(f"control_{i}",truth,x) for i,x in enumerate(starts,1)]
    elif mode=="real":
        worker=solve_real
        jobs=[(f"real_{i}",x,115) for i,x in enumerate(real_starts(c),1)]
    else:
        incumbent=sorted_feasible(read(HERE/"real_fits.tsv"))
        best=row_x(c,incumbent[0])
        worker=solve_real
        starts=[(best,divisor) for divisor in (1.,8.,30.,85.)]
        starts.extend((np.clip(best+rng.normal(0,s,len(best)),0,1),d)
                      for s,d in ((.0002,10.),(.001,30.),(.005,5.),(.015,10.)))
        starts.extend((row_x(c,r),12.) for r in incumbent[1:3])
        jobs=[(f"refine_{i}",x,115,divisor) for i,(x,divisor) in enumerate(starts,1)]
    path=HERE/("noiseless_fits.tsv" if mode=="noiseless" else "real_fits.tsv")
    records=read(path) if mode=="refine" else []
    with pool(workers) as work:
        for f in cf.as_completed([work.submit(worker,j) for j in jobs]):
            r=f.result()
            records.append(r)
            print(mode,r["origin"],r["loss"],r["max_projected_gradient"],
                  r.get("max_mechanistic_error",""),flush=True)
            write(path,sorted(records,key=lambda r:float(r["loss"])))
    if mode=="noiseless":
        best=sorted_feasible(records)[0]
        write(HERE/"noiseless_best.tsv",[best])
        errors=abs(row_x(c,best)-truth)
        write(HERE/"noiseless_recovery.tsv",[dict(parameter=p,normalized_error=e)
                                            for p,e in zip(c.names,errors)])
        assert max(errors[:14])<=.10, "Noiseless control did not recover all mechanisms"
    elif mode=="real":
        best=record_real(k,c,records)
        # Two local restarts check the lowest observed valley before profiling.
        x=row_x(c,best)
        jobs=[(f"verify_{i}",np.clip(x+rng.normal(0,s,len(x)),0,1),125)
              for i,s in enumerate((.00015,.0015,.01),1)]
        with pool(min(3,workers)) as work:
            for f in cf.as_completed([work.submit(worker,j) for j in jobs]):
                r=f.result()
                records.append(r)
                print(mode,r["origin"],r["loss"],r["max_projected_gradient"],flush=True)
                write(path,sorted(records,key=lambda r:float(r["loss"])))
        record_real(k,c,records)
    else:
        record_real(k,c,records)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("noiseless","real","refine"))
    parser.add_argument("--workers",type=int,default=8)
    args=parser.parse_args()
    run(args.mode,args.workers)
