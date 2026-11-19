"""Refine two open branch interfaces with exact rational witnesses and outward guards."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction

from verify_intervals import (arb, ctx, exact, read, table, count_or_none, OUT,
    SCHEDULE, parameter_points, forward, hill)
from verify_switch_boxes import branch_growth, branch_live


def do_case(args):
    row, points, changes = args
    group=row['category']
    coords_L={k:Fraction(v) for k,v in points[row['left']].items()}
    coords_R={k:Fraction(v) for k,v in points[row['right']].items()}
    ctx.prec=384
    a={k:exact(v) for k,v in coords_L.items()}
    b={k:exact(v) for k,v in coords_R.items()}
    ay,A,am,_,_=forward(a)
    by,B,bm,_,_=forward(b)
    if group=='mask':
        assert len(am^bm)==1
        o,n,m=next(iter(am^bm))
        from math import comb
        def guard(p):
            H=hill(p,exact(o));u=arb(n)/44
            mu=p['mu_hp']*H*u**p['gamma_mu']
            mis=p['p_mis_base']+p['p_misseg']*mu/(p['k_o_mis']+mu)
            return arb(comb(n,m))*mis**m*(1-mis)**(n-m)-exact('1e-8')
        guard_id=f'O={o};N={n};m={m}'
        rounds=24
    else:
        changed=[k for k in A if A[k][0]!=B[k][0]]
        assert len(changed)==1 and am==bm
        cohort,segment=changed[0]
        sr=next(r for r in SCHEDULE if (r['cohort'],r['segment_id'])==(cohort,segment))
        pid=sr['data_ids'].split(';')[0]
        days=[A[changed[0]][0],B[changed[0]][0]]
        target=count_or_none(sr['target_live_cells'])
        assert target is not None
        def guard(p):
            totals,_=branch_live(p,A,am,pid,days)
            return abs(totals[days[0]]-target)-abs(totals[days[1]]-target)
        guard_id=f'{pid}; days={days}'
        rounds=15

    def side(x):
        assert not x.contains(0), f'Indeterminate exact-guard midpoint {x}'
        return x>0
    left_sign=side(guard(a))
    assert side(guard(b))!=left_sign
    for _ in range(rounds):
        mid={k:(coords_L[k]+coords_R[k])/2 for k in coords_L}
        u={k:exact(v) for k,v in mid.items()}
        if side(guard(u))==left_sign:
            coords_L=mid
        else:
            coords_R=mid
    aa={k:exact(v) for k,v in coords_L.items()}
    bb={k:exact(v) for k,v in coords_R.items()}
    ra,A2,am2,_,_=forward(aa)
    rb,B2,bm2,_,_=forward(bb)
    assert (A,am)==(A2,am2) and (B,bm)==(B2,bm2)
    box={k:aa[k].union(bb[k]) for k in aa}
    keys=sorted((k for k in ra if k[0]=='growth'),
                key=lambda k:abs(rb[k]-ra[k]).abs_lower(),reverse=True)[:3]
    for key in keys:
        l=branch_growth(box,A,am,key[1])
        r=branch_growth(box,B,bm,key[1])
        gap=r-l
        if not gap.contains(0):
            break
    parameters=[dict(category=group,side=side_name,coordinate=k,exact_natural_value=str(v))
                for side_name,vec in [('left',coords_L),('right',coords_R)]
                for k,v in vec.items()]
    result=dict(category=group,guard=guard_id,rational_bisections=rounds,
          result='branches_separated_in_common_box' if not gap.contains(0) else 'box_overlap',
          passage=key[1],left_growth_interval=str(l),right_growth_interval=str(r),
          right_minus_left_interval=str(gap),pruning_flips=len(am^bm),
          endpoint_flips=sum(A[k][0]!=B[k][0] for k in A),
          fallback_flips=sum(A[k][2]!=B[k][2] for k in A))
    return result,parameters


def main():
    rows=[r for r in read(OUT/'narrow_pairs.tsv') if r['category'] in ('endpoint','mask')]
    points=parameter_points('narrow_pair_parameters.tsv')
    changes=read(OUT/'paired_branch_changes.tsv')
    with ProcessPoolExecutor(max_workers=2) as pool:
        records=list(pool.map(do_case,((r,points,changes) for r in rows)))
    table(OUT/'refined_branch_box_certificates.tsv',[t[0] for t in records])
    table(OUT/'refined_exact_box_endpoints.tsv',[v for t in records for v in t[1]])
    for r,_ in records:
        print(r['category'],r['result'],r['right_minus_left_interval'],flush=True)


if __name__=='__main__':
    main()
