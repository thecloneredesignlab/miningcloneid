"""Replay every simulated estimate and audit pairing and scoring."""
from hashlib import sha256

import numpy as np
from scipy.stats import binomtest

from recovery import MODE, generate, target_names
from joint_model import BASE, HERE, ALIGNMENT, read, row_x, write
from refit import case


def validate():
    k,c,_ = case(MODE)
    data = read(HERE/"recovery_results.tsv")
    previous = {(r["scenario"],int(r["dataset"])):float(r["innovation_rms"])
                for r in read(BASE/"draw_manifest.tsv")}
    matches = {(int(r["case_id"]),int(r["dataset"])) for r in data}
    assert len(data)==len(matches)==200
    assert matches=={(j,i) for j in (1,2) for i in range(1,101)}
    checks=[]
    counts={}
    names=target_names()
    for r in data:
        j,i=int(r["case_id"]),int(r["dataset"])
        t, truth, y=generate(j,i)
        assert t["scenario"]==r["scenario"]
        innovations=np.random.default_rng(np.random.SeedSequence([20261118,221,j,i])).normal(size=k.n)
        assert abs(previous[t["scenario"],i]-float(np.std(innovations))) < 1e-14
        x=row_x(c,r)
        p=c.unmap(x)[0]
        out,_,f,_=k.evaluate_joint(p,y)
        assert abs(f-float(r["loss"]))<1e-9
        assert abs(out[1]-float(r["donor_margin"]))<1e-10
        ok=r["qualified"]=="True"
        if ok:
            assert out[1]>=0 and float(r["stationarity"])<=.025
            assert float(r["independent_lr_gap"])<=.10
        checks.append(dict(scenario=t["scenario"],case_id=j,dataset=i,qualified=ok,
                           final_loss=f,donor_margin=out[1],noise_seed_pair_verified=True))
        for name in names:
            ix=c.names.index(name)
            err=abs(x[ix]-truth[ix])
            result=counts.setdefault((name,j),[0,0])
            result[0]+=ok
            result[1]+=ok and err<=.10
    write(HERE/"recovery_replay.tsv",checks)
    metrics=[]
    for (name,j),(verified,hits) in sorted(counts.items()):
        ci=binomtest(hits,100).proportion_ci()
        metrics.append(dict(parameter=name,case_id=j,assigned=100,numerically_verified=verified,
                            recovered=hits,rate=hits/100, binomial_95_lower=ci.low,
                            binomial_95_upper=ci.high))
    write(HERE/"recovery_metrics.tsv",metrics)
    inputs=[BASE/name for name in ("ar1_model.py","ar1_kernel.cpp","truths.tsv",
                                    "observed_best.tsv","practical_verdicts.tsv",
                                    "draw_manifest.tsv")]
    inputs.extend(ALIGNMENT/name for name in ("cell_alignment.tsv","sample_alignment.tsv"))
    inputs.extend(HERE/name for name in ("joint_kernel.cpp","joint_model.py","refit.py",
                                          "profiles.py","recovery.py"))
    write(HERE/"input_hashes.tsv",[dict(path=str(p),sha256=sha256(p.read_bytes()).hexdigest())
                                  for p in inputs])
    print(metrics,flush=True)


if __name__=="__main__":
    validate()
