"""Audit the eligible-object denominator, forward map and objective gradients."""
from pathlib import Path
import numpy as np
from scipy.stats import betabinom
from modal import HERE,PLAN,STAGE3,case,read,write


def run():
    ref=read(STAGE3/"best_autosomes_joint.tsv")[0]
    pilot=Path.home()/("Documents/CloneID/Deployments/sum159-red-pilot-diagnostic-001/"
                       "operator_project/calibration")
    denominator={r["source_id"]:r for r in read(pilot/"plan/rds_morphology_source_audit.tsv")}
    count={r["source_record_id"]:r for r in read(
        pilot/"prepared_inputs/count_growth_observations.tsv")}
    checks=[]
    for r in read(PLAN/"paired_death_readout_audit.tsv"):
        i=r["record_id"]
        m=denominator.get(i)
        if m is not None:
            assert sum(int(m[j]) for j in ("live_cell","dead_cell","multinucleated_cell")) == int(m["n_eligible_cells"])
            assert int(m["dead_cell"])==int(r["dead_count"])
            assert int(m["n_eligible_cells"])==int(r["eligible_objects"])
        p=count.get(i,{})
        checks.append(dict(r,object_partition_replayed=m is not None,
                           available_classifier_fields=m["image_bundle_count"] if m else "",
                           count_semantics=p.get("count_semantics",""),
                           exact_field_overlap="unverified",
                           primary=r["classification_status"]=="classified_full_promoted"))
    assert len(checks)==195 and sum(r["object_partition_replayed"] for r in checks)==46
    write(HERE/"image_denominator_audit.tsv",checks)
    rows=[]
    for mode in ("baseline","a7k","death"):
        k,c,y=case(mode)
        x=c.anchor(ref)
        p,j=c.unmap(x)
        ar,f,g,ky,de=k.parts(p,y,True)
        if mode=="baseline":
            assert np.isclose(f,k.n*np.log(float(ref["loss"])),atol=1e-9,rtol=0.)
        if mode=="death":
            channels,_=k.readouts(p)
            prob,kappa,*_=k.death_prediction(p,channels,y)
            n=k.eligible
            d=k.dead
            from scipy.special import gammaln
            logchoose=gammaln(n+1)-gammaln(d+1)-gammaln(n-d+1)
            assert abs(-de+logchoose.sum()-np.sum(betabinom.logpmf(d,n,kappa*prob,kappa*(1-prob))))<1e-8
        x=np.clip(x,.007,.993)
        h=3e-6
        z,j=c.unmap(x)
        analytic=c.gradient(k.parts(z,y,True)[2],j)
        numerical=np.array([(k.parts(c.unmap(x+u*h)[0],y)[1]
                            -k.parts(c.unmap(x-u*h)[0],y)[1])/(2*h) for u in np.eye(len(x))])
        err=max(abs(analytic-numerical)/np.maximum(1.,abs(numerical)))
        assert err<5e-4,(err,analytic,numerical)
        rows.append(dict(mode=mode,stage3_anchor_loss=f,stage3_anchor_margin=ar[1],
                         counts=k.n,kary_cells=k.kary.sum(),death_harvests=len(k.dead),
                         max_scaled_grad_error=err))
        print(rows[-1],flush=True)
    write(HERE/"likelihood_audit.tsv",rows)


if __name__=="__main__":
    run()
