import math
from time import monotonic

from flint import arb, ctx
import numpy as np

from practical import DESIGN, HERE, Configuration, read, write


def value(r):
    return math.ldexp(int(r["mid_man"]), int(r["mid_exp"]))


def check(k):
    ctx.prec = 192
    rng = np.random.default_rng(2026111801)
    comp = rng.normal(size=len(k.ids))
    positions = {p: i for i,p in enumerate(k.ids)}
    receipts = []
    for root in ("hist", "gaussian"):
        c = Configuration(k, gaussian=root=="gaussian", split=True, tilt=True)
        f, actual, diag, mask, g = k.predict(k.anchor, c.gaussian, comp)
        assert diag[0] == 0 and diag[1] > 0
        ref_counts = {}
        expected_gradient = np.zeros(46)
        column_l1 = np.zeros(46)
        for cohort in ("2N","4N"):
            for r in read(DESIGN/root/f"{cohort}_counts_replay.tsv"):
                ref_counts[r["record_id"]] = float(arb(r["expected_visible_witness_cells"]).mid())
            for r in read(DESIGN/root/f"{cohort}_ols_derivatives.tsv"):
                i = [x["parameter"] for x in k.rows].index(r["column"])
                v = comp[positions[r["output_id"]]]*value(r)
                expected_gradient[i] += v
                column_l1[i] += abs(v)
        relative = float(max(abs(actual[i]/ref_counts[s["record_id"]]-1) for i,s in enumerate(k.slots)))
        err = float(np.max(np.abs(g-expected_gradient)/np.maximum(column_l1,1e-20)))
        assert relative < 2e-13 and err < 1e-6, (root, relative, err)
        p,_ = c.unmap(c.x_anchor)
        loss,d,supply=k.evaluate(p,k.y,c)
        check_loss = 0.5*np.sum((f-k.y)**2/k.var)
        assert abs(loss-check_loss)<1e-10
        receipts.append(dict(root=root, max_count_relative_error=relative,
                             max_gradient_discrepancy_over_L1=err, cost_abs_error=abs(loss-check_loss),
                             donor_surplus_fraction=supply))
        print(receipts[-1],flush=True)
    write(HERE/"kernel_validation.tsv",receipts)
    t=monotonic()
    for _ in range(100):k.evaluate(p,k.y,c)
    print("gradient seconds", (monotonic()-t)/100,flush=True)
    return receipts
