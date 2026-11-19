"""Assess the shared count loading and summarize the extra supply requirement."""
from concurrent.futures import as_completed
import numpy as np

from modal import HERE,PLAN,STAGE3,case,read,pool,write
from fit import run_one,take

mode="stock_death"


def audit_supply():
    k,c,y=case("stock_baseline")
    ref=read(STAGE3/"best_autosomes_joint.tsv")[0]
    z=c.anchor(ref)
    ph,_=c.unmap(z)
    channels,_=k.readouts(ph)
    ids={r["record_id"]:i for i,r in enumerate(k.slots)}
    out=[]
    for a in read(PLAN/"A7K_matched_readouts.tsv"):
        donor=k.slots[ids[a["donor_record"]]]
        main=sum(float(r["fixed_live_cells"]) for r in k.graph
                 if r["parent_id"]==donor["passage_id"] and int(r["donor_step"])==100)
        supplement=float(a["fixed_recorded_seed"])
        predicted=float(channels[0,ids[a["donor_record"]]]*1e12)
        out.append(dict(branch=a["physical_passage"],donor=a["donor_record"],
                        fixed_branch_seed=supplement,registered_main_demand=main,
                        stage3_predicted_donor_viable=predicted,
                        combined_margin=(predicted-main-supplement)/(main+supplement)))
    write(HERE/"A7K_supply_sensitivity.tsv",out)


def audit_shared():
    k,c,y=case(mode)
    base=read(HERE/f"best_{mode}.tsv")[0]
    anchor=c.anchor(base)
    options=[(mode,"loading_fixed_zero",anchor,{c.names.index("shared_count_loading"):.5},100,50.),
             (mode,"loading_fixed_zero_alt",anchor,{c.names.index("shared_count_loading"):.5},108,20.)]
    with pool(2) as ex:
        rows=[f.result() for f in as_completed(ex.submit(run_one,j) for j in options)]
    r=take(rows)
    r["lr_compared_free_loading"]=float(r["loss"])-float(base["loss"])
    write(HERE/"shared_loading_sensitivity.tsv",[r])


if __name__=="__main__":
    audit_supply()
    audit_shared()
