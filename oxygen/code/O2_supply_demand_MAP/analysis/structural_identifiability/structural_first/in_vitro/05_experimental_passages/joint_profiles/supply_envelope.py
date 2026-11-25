"""Check a necessary donor-supply inequality for the lam_max scan."""
from collections import defaultdict
from fractions import Fraction
from math import expm1, log
from pathlib import Path
import sys

from scan import DESIGN, read, write

DT = Fraction(1, 20)


def demand_groups():
    graph=read(DESIGN/"passage_graph.tsv")
    seed={r["passage_id"]:Fraction(r["fixed_seed_live_objects"]) for r in graph}
    demand=defaultdict(Fraction)
    for row in graph:
        parent=row["donor_passage"]
        if not parent or "_A7K" in row["passage_id"]:
            continue
        steps=Fraction(row["donor_model_day"])/DT
        assert steps.denominator==1 and steps>0
        demand[parent,int(steps)]+=seed[row["passage_id"]]
    assert len(demand)==110
    return seed,demand


def collect(out):
    seed,demand=demand_groups()
    grid=read(out/"profiles"/"lam_max"/"grid.tsv")
    detailed=[]
    for row in grid:
        rate=Fraction(row["value"])
        assert 0<=rate*DT<=1
        # Dropping death and truncation gives at most 1 + dt*lambda_max
        # per Euler update. The column bound also holds under clipping.
        candidates=[(seed[parent]*(1+DT*rate)**steps/required,parent,steps,required)
                    for (parent,steps),required in demand.items()]
        bound,parent,steps,required=min(candidates)
        gap=bound-1
        fails=gap < -Fraction(1, 10**6)
        if fails:
            assert not row["loss"], row
        lower=expm1(log(float(required/seed[parent]))/steps)/float(DT)
        detailed.append(dict(parameter="lam_max",value=row["value"],donor_passage=parent,euler_steps=steps,
            starting_live_count=float(seed[parent]),required_live_count=float(required),
            upper_live_count=float(bound*required),best_case_relative_margin=float(gap),
            donor_necessary_minimum_rate=lower,violates_live_supply_envelope=fails,
            numerical_grid_status=row["status"]))
    write(out/"lam_max_supply_envelope.tsv",detailed)


if __name__=="__main__":
    collect(Path(sys.argv[1]))
