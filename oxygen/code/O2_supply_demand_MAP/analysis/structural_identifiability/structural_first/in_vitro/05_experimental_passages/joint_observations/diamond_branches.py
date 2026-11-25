"""Cross-check the fitted point's native and interval pruning branches."""
import argparse
import os
from pathlib import Path
from time import process_time
import sys

from paths import FROZEN_DESIGN, PROFILES, output_directory


def native(out):
    os.environ["JOINT_KERNEL"]="1"
    sys.path.insert(0,str(PROFILES))
    from profile_model import ProfileModel
    from tv_review import fingerprint
    from runtime import read, write
    import numpy as np
    model=ProfileModel("J")
    p=np.array([float(r["value"]) for r in read(out/"native_vector.tsv")])
    flags=fingerprint(model,p)
    write(out/"native_masks.tsv",[dict(oxygen=str(o),signature=f"{m:x}") for o,m in zip(model.oxygens,flags)])


def audit(out):
    import joint_rank as rank
    from certify_diamond import parameters
    from flint import ctx
    from algebra import exact
    ctx.prec=512
    p={n:exact(v) for n,v in parameters(out).items()}
    native_masks={float(r["oxygen"]):r["signature"] for r in rank.read(out/"native_masks.tsv")}
    graph=rank.read(FROZEN_DESIGN/"passage_graph.tsv")
    oxygen=sorted({r["oxygen_pct"] for r in graph},key=float)
    assert {float(o) for o in oxygen}==native_masks.keys()
    math=rank.inherit.cert.model
    math.rates=rank.inherit.cert.fixed_rates
    audited=[]
    for o in oxygen:
        op=math.operator(p,exact(o))
        mask=1469598103934665603
        for j,r in enumerate(op["masks"]):
            for m,retained in enumerate(r["mask"]):
                if retained=="1":
                    mask=((mask^(1+m+j*154))*1099511628211) % 2**64
        gap=min(r["minimum_margin"] for r in op["masks"])
        lower=min((1-exact("0.05")*(v["lam"]+v["mu"])).lower() for v in op["values"])
        upper=max(v["survival"].upper() for v in op["values"])
        sha=op["mask_sha256"]
        assert gap>0 and lower>0 and upper<1 and f"{mask:x}"==native_masks[float(o)]
        for bits in (384,512):
            passage_masks={r["passage"]:r["operator_mask"] for r in rank.read(out/f"arb_{bits}"/"branch_checks.tsv")}
            assert all(passage_masks[r["passage_id"]]==sha for r in graph if r["oxygen_pct"]==o)
        audited.append(dict(oxygen=o,state_count=len(op["values"]),minimum_pruning_distance=str(gap),
              positive_diagonal_lower=str(lower),maximum_survival_upper=str(upper),structural_mask=sha,
              native_mask=native_masks[float(o)],branch_match=True))
    rank.write(out/"operator_branches.tsv",audited)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("stage",choices=("native","interval"))
    cli.add_argument("out",type=output_directory)
    a=cli.parse_args()
    start=process_time()
    (native if a.stage=="native" else audit)(a.out)
    (a.out/f"branch_{a.stage}_resources.tsv").write_text(f"cpu_seconds\n{process_time()-start}\n")
