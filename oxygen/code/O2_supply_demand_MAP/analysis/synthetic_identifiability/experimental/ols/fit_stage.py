"""First-stage recovery controls and observed-data fit multistart."""
import argparse
import concurrent.futures
import multiprocessing
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from practical import (HERE, Configuration, Kernel, optimize, result_row,
                       row_x, starting_points, write)

KERNEL = None


def model():
    global KERNEL
    if KERNEL is None:
        KERNEL = Kernel()
    return KERNEL


def fit_worker(task):
    label, mode, target, start, limit = task
    kernel = model()
    config = Configuration(kernel)
    if mode == "noiseless":
        # Trust-region sum of squares with scale-sensitive finite differences
        # provides an optimizer independent of the main adjoint-based fits.
        def residual(x):
            p,_ = config.unmap(x)
            pred,_,diag,_,_ = kernel.predict(p)
            if diag[0] > 0:
                return (pred - target)/kernel.scale + 15*(1+diag[0])
            return (pred - target)/kernel.scale
        res = least_squares(residual, start, bounds=(0,1), max_nfev=limit, x_scale="jac",
                            xtol=1e-12, ftol=1e-12, gtol=1e-13, diff_step=1e-5)
        p,_ = config.unmap(res.x)
        y,_,diag,_,_ = kernel.predict(p)
        f = dict(x=res.x, loss=0.5*float(np.dot(res.fun,res.fun)), supply=diag[1], gradient=res.optimality,
                 calls=res.nfev, iteration=res.nit if hasattr(res,"nit") else res.nfev,
                 scipy_success=bool(res.success), termination=str(res.message), seconds="")
    else:
        f = optimize(kernel, config, target, start, limit)
        y,raw,diag,_,_ = kernel.predict(config.unmap(f["x"])[0])
    return result_row(config,f,origin=label,mode=mode,
                      maximum_absolute_slope_residual=float(max(abs(y-target))))


def run(mode, workers):
    kernel = model()
    config = Configuration(kernel)
    job = []
    if mode == "noiseless":
        x = config.x_anchor.copy()
        x = np.clip(x,.015,.985)
        p,_ = config.unmap(x)
        y,_,diag,_,_ = kernel.predict(p)
        assert diag[0] == 0
        check,f,s=kernel.evaluate(p,y,config)
        assert check<1e-20
        write(HERE/"noiseless_truth.tsv",[dict(parameter=q,unit_coordinate=float(v),natural_value=p[j])
                                         for j,q,v in zip(config.active,config.names,x)])
        write(HERE/"noiseless_truth_predictions.tsv",[dict(passage_id=k,slope=z) for k,z in zip(kernel.ids,y)])
        # Distinct starting points, none set to the noiseless truth.
        rng=np.random.default_rng(210188)
        for n,d in enumerate((.02,.04,.07,.10,.12,.20)):
            start=np.clip(config.x_anchor+rng.normal(0,d,16),.001,.999)
            job.append((f"control_{n+1}",mode,y,start,60))
    elif mode == "real":
        for n,start in enumerate(starting_points(kernel,config,401,n=10)):
            job.append((f"real_{n+1}",mode,kernel.y,start,320))
    receipts=[]
    path=HERE/f"{mode}_fits.tsv"
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context("spawn")) as pool:
        results={pool.submit(fit_worker, item):item[0] for item in job}
        for future in concurrent.futures.as_completed(results):
            row=future.result()
            print(mode,row["origin"],row["loss"],row["maximum_absolute_slope_residual"],
                  row["max_projected_gradient"],row["termination"],flush=True)
            receipts.append(row)
            write(path,sorted(receipts,key=lambda q:q["loss"]))
    fits=[r for r in receipts if r["min_aggregate_donor_margin_fraction"]>=0 and r["loss"]<100]
    assert fits
    best=min(fits,key=lambda r:r["loss"])
    write(HERE/f"{mode}_best.tsv",[best])
    if mode=="noiseless":
        err=np.abs(row_x(config,best)-x)
        write(HERE/"noiseless_summary.tsv",[dict(parameter=name,normalized_absolute_error=error)
                                            for name,error in zip(config.names,err)])
    else:
        pred,_,diag,_,_=kernel.predict(config.unmap(row_x(config,best))[0])
        write(HERE/"real_residuals.tsv",[dict(passage_id=key,y=y,prediction=z,error=y-z,variance_multiplier=v)
                                        for key,y,z,v in zip(kernel.ids,kernel.y,pred,kernel.var)])
        write(HERE/"model_parameters.tsv",kernel.rows)
    return best


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=["real","noiseless"])
    parser.add_argument("--workers",type=int,default=5)
    a=parser.parse_args()
    run(a.mode,a.workers)
