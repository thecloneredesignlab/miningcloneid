"""Disjoint output enclosures for two fixed branch extensions in a common tight box."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction

from verify_intervals import (arb, arb_mat, ctx, exact, read, table, count_or_none,
    OUT, N, GRID, DT, CAPACITY, FALLBACK, SCHEDULE, operator, parameter_points, forward)


def branch_live(p, branches, mask, passage, terminal_days=None):
    dt=exact(DT)
    by_key={(r['cohort'],r['segment_id']):r for r in SCHEDULE}
    r=next(r for r in SCHEDULE if passage in r['data_ids'].split(';'))
    if terminal_days is None:
        terminal_days = [branches[(r['cohort'],r['segment_id'])][0]]
    final_steps = {int(Fraction(day)/Fraction(DT)):day for day in terminal_days}
    final_totals = {}
    path=[(r,max(terminal_days,key=Fraction))]
    while r['parent_segment_id']!='NA':
        same_key=(r['cohort'],r['segment_id'])
        day=branches[same_key][1]
        r=by_key[(r['cohort'],r['parent_segment_id'])]
        path.append((r,day))
    opt={}
    x=None
    size=None
    for path_index,(r,day) in enumerate(reversed(path)):
        size=count_or_none(r['initial_cells'])
        if size is None:
            size=arb(FALLBACK)/CAPACITY
        if x is None:
            cohort=r['cohort']
            mean,sd=p['init_mean_'+cohort],p['init_sd_'+cohort]
            gaussian=[(-((n-mean)/sd)**2/2).exp() for n in GRID]
            x=arb_mat([[size*v/sum(gaussian)] for v in gaussian])
        else:
            x=x*(size/sum(x.entries()))
        o=r['oxygen_pct']
        if o not in opt:
            G,death,_,valid=operator(p,o,fixed_mask=mask)
            if not valid:
                raise ValueError('Positive live update not separated in box')
            diag=[G[i,i] for i in range(N)]
            for i in range(N):
                G[i,i]=arb(0)
            opt[o]=(G,diag,death)
        G,diag,death=opt[o]
        for step in range(1,1+int(Fraction(day)/Fraction(DT))):
            c=1-sum(x.entries())
            if c<0:
                c=arb(0)
            if not c>=0:
                raise ValueError('Logistic sign unresolved in common box')
            inflow=G*x
            # Combine the diagonal loss terms before interval multiplication.
            x=arb_mat([[(1+dt*(c*diag[i]-death[i]))*x[i,0]+dt*c*inflow[i,0]]
                        for i in range(N)])
            if path_index==len(path)-1 and step in final_steps:
                final_totals[final_steps[step]]=sum(x.entries())
    return final_totals,size


def branch_growth(p, branches, mask, passage):
    totals,size = branch_live(p,branches,mask,passage)
    day,total = next(iter(totals.items()))
    return (total/size).log()/exact(day)


def analyze(args):
    r, left, right, expected=args
    ctx.prec=384
    a={k:exact(v) for k,v in left.items()}
    b={k:exact(v) for k,v in right.items()}
    try:
        ay,A,am,_,_=forward(a)
        by,B,bm,_,_=forward(b)
        e=sum(A[k][0]!=B[k][0] for k in A)
        z=sum(A[k][2]!=B[k][2] for k in A)
        q=len(am^bm)
        if e!=int(r['endpoint']) or z!=int(r['reset_mode']) or q!=int(r['mask']):
            raise ValueError('Tight real-arithmetic branch difference does not match native claim')
        if A!=expected[r['left']] or B!=expected[r['right']]:
            raise ValueError('Tight real-arithmetic endpoint or boundary branch does not match native')
        box={k:a[k].union(b[k]) for k in a}
        candidates=sorted((k for k in ay if k[0]=='growth'),
                          key=lambda k:abs(by[k]-ay[k]).abs_lower(),reverse=True)[:3]
        best=None
        for key in candidates:
            left_band=branch_growth(box,A,am,key[1])
            right_band=branch_growth(box,B,bm,key[1])
            diff=right_band-left_band
            best=dict(category=r['category'],left=r['left'],right=r['right'],
                result='branches_separated_in_common_box' if not diff.contains(0) else 'box_overlap',
                passage=key[1],left_growth_interval=str(left_band),right_growth_interval=str(right_band),
                right_minus_left_interval=str(diff),normalized_unit_distance=r['unit_distance'],
                endpoint_changes=e,mask_changes=q,fallback_changes=z,reason='')
            if not diff.contains(0):
                break
        return best
    except (ValueError,AssertionError,ArithmeticError) as error:
        return dict(category=r['category'],left=r['left'],right=r['right'],result='unresolved',
            passage='',left_growth_interval='',right_growth_interval='',right_minus_left_interval='',
            normalized_unit_distance=r['unit_distance'],endpoint_changes='',mask_changes='',fallback_changes='',
            reason=str(error))


def main():
    rows=read(OUT/'narrow_pairs.tsv')
    points=parameter_points('narrow_pair_parameters.tsv')
    expected={}
    native=read(OUT/'narrow_native_branches.tsv')
    if (OUT/'targeted_pairs.tsv').is_file():
        rows += [r for r in read(OUT/'targeted_pairs.tsv') if r['kind']=='targeted_fallback_straddle']
        points.update(parameter_points('targeted_pair_parameters.tsv'))
        native += read(OUT/'targeted_native_branches.tsv')
    for row in native:
        expected.setdefault(row['id'],{})[(row['cohort'],row['segment'])]=(
            row['endpoint_day'],row['parent_day'],row['inoculum_fallback']=='TRUE')
    with ProcessPoolExecutor(max_workers=3) as pool:
        evidence=list(pool.map(analyze,((r,points[r['left']],points[r['right']],expected) for r in rows)))
    table(OUT/'branch_box_certificates.tsv',evidence)
    for r in evidence:
        print(r['category'],r['result'],r['right_minus_left_interval'] or r['reason'],flush=True)


if __name__=='__main__':
    main()
