"""Replay the old likelihood, both strict profiles and 200 recovery decisions."""
import numpy as np

from tilt_model import HERE,PREV,TARGETS,case,read,receipt,write
from draws import draw


def run():
    k,c,y=case("strict")
    old=read(PREV/"best_stock_a7k.tsv")[0]
    x=c.anchor(old)
    p,_=c.unmap(x)
    ar,f,g,_,_=k.parts(p,y,True)
    assert abs(f-float(old["loss"]))<1e-9
    assert abs(ar[1]-float(old["donor_margin"]))<1e-9
    calc=k.composition(p,True)
    assert max(calc[0])<1e-15
    expected=[float(a["fixed_recorded_seed"]) for a in k.a7]
    assert max(abs(calc[3]-expected))<1e-7
    result=receipt(k,c,y,dict(x=x,gradient=float(old["stationarity"]),success=True,
                   seconds=0),mode="strict",origin="exact_strict_replay")
    write(HERE/"best_strict.tsv",[result])
    verified=[]
    for r in read(PREV/"profile_points_stock_a7k.tsv"):
        if r["parameter"] not in TARGETS:
            continue
        xx=c.anchor(r)
        pp,_=c.unmap(xx)
        ar,loss,*_=k.parts(pp,y)
        assert abs(loss-float(r["loss"]))<1e-8
        assert abs(ar[1]-float(r["donor_margin"]))<1e-9
        verified.append(dict(parameter=r["parameter"],coordinate=r["coordinate"],
                              loss=loss,lr=loss-f,donor_margin=ar[1]))
    write(HERE/"strict_profile_replay.tsv",verified)
    write(HERE/"intervals_strict.tsv",[dict(r,mode="strict") for r in
          read(PREV/"intervals_stock_a7k.tsv") if r["parameter"] in TARGETS])
    source=PREV/"arrays_stock_a7k_expanded/frozen.npz"
    with np.load(source) as hist_draw:
        exact={tuple(key):v for key,v in zip(hist_draw["indices"],hist_draw["targets"])}
    check=[]
    for r in read(PREV/"recovery_stock_a7k.tsv"):
        j=int(r["case_id"]); i=int(r["replicate"])
        draw_y,truth,_=draw("strict",j,i)
        if (j,i) in exact:
            assert np.max(np.abs(draw_y-exact[j,i]))<1e-11
        pp,_=c.unmap(c.anchor(r))
        ar,loss,*_=k.parts(pp,draw_y)
        delta=abs(loss-float(r["loss"]))
        assert delta<1e-7 and ar[1]>=-3e-12 and r["validated"]=="True"
        z=c.names.index("p_misseg")
        check.append(dict(case_id=j,replicate=i,loss_error=delta,
                          p_recovered=abs(c.anchor(r)[z]-truth[z])<=.10))
    assert len(check)==200
    assert [sum(r["case_id"]==j and r["p_recovered"] for r in check) for j in (1,2)]==[94,88]
    write(HERE/"strict_recovery_replay.tsv",check)

    rows=[]
    for mode in ("group_010","group_025","group_050","branch_050"):
        model,config,goal=case(mode)
        z=np.clip(config.anchor(old),.003,.997)
        z[config.dim:]=.5+np.linspace(-.006,.009,len(z)-config.dim)
        pp,jac=config.unmap(z)
        calc=model.composition(pp,True)
        no_tilt=pp.copy()
        no_tilt[46:54]=0.
        assert np.array_equal(model.predict_logs(pp),model.predict_logs(no_tilt))
        assert np.array_equal(model.readouts(pp)[1][:2],model.readouts(no_tilt)[1][:2])
        assert np.max(abs(calc[3]-expected))<1e-7
        assert np.max(abs(calc[0]-.5*np.abs(calc[2]-calc[1]).sum(axis=1)))<1e-13
        h=1e-6
        g=config.gradient(model.parts(pp,goal,True)[2],jac)
        selected=np.arange(110,114,dtype=np.int32)
        vals,d=model.supply(pp,pick=selected)
        dev=np.stack([config.gradient(a,jac) for a in d])
        ng=[]; nd=[]
        for v in np.eye(len(z)):
            up=config.unmap(z+h*v)[0]
            down=config.unmap(z-h*v)[0]
            ng.append((model.parts(up,goal)[1]-model.parts(down,goal)[1])/(2*h))
            nd.append((model.supply(up)[0][selected]-model.supply(down)[0][selected])/(2*h))
        ng=np.array(ng); nd=np.array(nd).T
        g_error=np.max(abs(g-ng)/np.maximum(1.,abs(ng)))
        tv_error=np.max(abs(dev-nd)/np.maximum(1.,abs(nd)))
        rows.append(dict(mode=mode,scaled_likelihood_gradient_error=g_error,
                         scaled_tv_gradient_error=tv_error,max_seed_error=max(abs(calc[3]-expected))))
        assert g_error<8e-5 and tv_error<8e-5,(mode,g_error,tv_error)
    write(HERE/"numerical_audit.tsv",rows)
    print("Strict replay:",f,len(verified),"profile points, 200 simulations, rates 94/88",flush=True)
    print("Audited TV gradients:",rows,flush=True)


if __name__=="__main__":
    run()
