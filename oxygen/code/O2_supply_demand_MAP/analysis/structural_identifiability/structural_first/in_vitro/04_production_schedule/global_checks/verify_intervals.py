"""Rigorous forward replay and pair discrimination for the Step 03 endpoint study."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from pathlib import Path
import sys

OUT = Path(__file__).resolve().parent
STEP2 = OUT.with_name('20261117_culture_actual_schedule_step02')
sys.path.insert(0, str(STEP2))
from certify import (arb, arb_mat, ctx, exact, nonzero, encode, read, table, count_or_none,
                     GRID, N, DT, CAPACITY, FALLBACK, CONTRACT, STEP1)

COLS = tuple(r['parameter'] for r in read(STEP2 / 'anchor.tsv') if r['parameter'] != 'sigma_growth')
BOUNDS = {r['parameter']:(Fraction(r['lower']), Fraction(r['upper']))
          for r in read(STEP2 / 'anchor.tsv')}
SCHEDULE = read(CONTRACT / 'experimental_schedule.tsv')
PASSAGES = {r['passage_id']: r for r in read(STEP1 / 'passage_crosswalk.tsv')}


def parameter_points(filename):
    points = {}
    for r in read(OUT / filename):
        value = Fraction(float.fromhex(r['binary64']))
        for bound in BOUNDS[r['parameter']]:
            if float(bound) == float(value):
                value = bound
                break
        assert BOUNDS[r['parameter']][0] <= value <= BOUNDS[r['parameter']][1]
        points.setdefault(r['id'], {})[r['parameter']] = str(value)
    assert all(set(p) == set(COLS) for p in points.values())
    return points


def hill(p, O):
    if O.is_zero():
        assert p['n_O'] > 0 and p['O2_crit'] > 0
        return arb(1)
    C, n = p['O2_crit'], p['n_O']
    return C**n/(C**n+O**n)


def operator(p, oxygen, fixed_mask=None):
    h = hill(p, exact(oxygen))
    G = arb_mat(N, N)
    death = []
    mask = []
    eps = exact('1e-8')
    for j, n in enumerate(GRID):
        ratio = arb(n)/44
        lam = p['lam_max'] / (1+p['alpha_o2']*h*ratio**p['gamma_growth'])
        mu = p['mu_hp']*h*ratio**p['gamma_mu']
        mis = p['p_mis_base']+p['p_misseg']*mu/(p['k_o_mis']+mu)
        assert lam > 0 and mu > 0 and mis > 0 and mis < 1
        S = p['buffer_smax'] * (-p['buffer_beta'] * (1/ratio)**p['buffer_n_exp']).exp()
        assert S.is_zero() or 0 < S < 1
        death.append(mu)
        G[j,j] -= lam
        if 2*n <= GRID[-1]:
            G[2*n-GRID[0],j] += lam*p['p_wgd']
        pr = (1-mis)**n
        for m in range(n+1):
            if m:
                pr = pr * (n-m+1) * mis / (m*(1-mis))
            keep = pr > eps if fixed_mask is None else (oxygen,n,m) in fixed_mask
            if keep:
                mask.append((oxygen,n,m))
                survival = S**m if m else arb(1)
                for dest in (n-m,n+m):
                    if dest in GRID:
                        G[dest-GRID[0],j] += lam*(1-p['p_wgd'])*pr*survival
            elif fixed_mask is None and not pr < eps:
                raise ValueError(f'Cutoff ambiguity at O={oxygen} N={n} m={m}')
    scale = exact(DT)
    sufficient_live_positivity = all(1-scale*(p['lam_max']+mu) > 0 for mu in death)
    return G, death, set(mask), sufficient_live_positivity


def select(states, target, required=None):
    days = [d for d in states if Fraction(d)>0] or list(states)
    if target is None:
        selected = days[-1]
    else:
        errors = {day:abs(states[day][1]-target) for day in days}
        selected = min(days,key=lambda day:errors[day].mid())
        for day in days:
            if day != selected and not errors[day] > errors[selected]:
                raise ValueError(f'Nearest-day branch unresolved, {day} versus {selected}')
    if required is None:
        return selected, False
    excess = states[selected][1]-required
    if excess > 0 or excess.is_zero():
        return selected, False
    if not excess < 0:
        raise ValueError('Target-day supply ambiguity')
    eligible = []
    for day in days:
        excess = states[day][1]-required
        if excess > 0 or excess.is_zero():
            eligible.append(day)
        elif not excess < 0:
            raise ValueError(f'Eligibility ambiguity at {day}')
    if not eligible:
        raise ValueError('Protocol infeasible')
    fallback = min(eligible,key=lambda day:states[day][1].mid())
    for day in eligible:
        if day != fallback and not states[day][1] > states[fallback][1]:
            raise ValueError(f'Eligible-state branch unresolved, {day} versus {fallback}')
    return fallback, True


def kernel(width):
    # Stable representation of K-I. The tiny diagonal correction is accumulated
    # through off-diagonal mass, avoiding subtraction of two O(1) numbers.
    result = arb_mat(N,N)
    kernel_dist = [(-arb(j*j)/(2*width**2)).exp() for j in range(1,N)]
    for j in range(N):
        off = {i:kernel_dist[abs(i-j)-1] for i in range(N) if i!=j}
        mass = sum(off.values())
        for i,v in off.items():
            result[i,j]=v/(1+mass)
        result[j,j]=-mass/(1+mass)
    return result


def forward(p, oxygen_override=None):
    dt = exact(DT)
    operators = {o:operator(p,o) for o in sorted(
        {r['oxygen_pct'] if oxygen_override is None else oxygen_override for r in SCHEDULE},key=Fraction)}
    correction = kernel(p['sigma_kary'])
    results, branches, values, latent, corrections = {}, {}, {}, {}, {}
    for j,r in enumerate(SCHEDULE):
        cohort,seg = r['cohort'],r['segment_id']
        oxygen = r['oxygen_pct'] if oxygen_override is None else oxygen_override
        size = count_or_none(r['initial_cells'])
        if size is None:
            size = arb(FALLBACK)/CAPACITY
        parent = cohort,r['parent_segment_id']
        parent_day, fallback = 'NA', False
        if parent in results:
            states,target = results[parent]
            parent_day,fallback = select(states,target,size)
            v,L = states[parent_day]
            x=v*(size/L)
        else:
            assert r['parent_segment_id']=='NA'
            m,s = p['init_mean_'+cohort],p['init_sd_'+cohort]
            values0 = [(-((arb(n)-m)/s)**2/2).exp() for n in GRID]
            x=arb_mat([[size*a/sum(values0)] for a in values0])
        G,death,_,positive_live = operators[oxygen]
        wanted = {int(Fraction(day)/Fraction(DT)):day for day in r['candidate_days'].split(';')}
        stop = int(Fraction(r['nominal_duration_days'])/Fraction(DT))
        states = {}
        for step in range(stop+1):
            L = sum(x.entries())
            if step in wanted:
                states[wanted[step]]=x,L
            if step==stop:
                break
            raw_crowd=1-L
            if raw_crowd>0:
                crowd=raw_crowd
            elif raw_crowd<0 or raw_crowd.is_zero():
                crowd=arb(0)
            else:
                raise ValueError('Logistic branch ambiguous')
            rates=G*x*crowd
            new=[x[i,0]+dt*(rates[i,0]-death[i]*x[i,0]) for i in range(N)]
            if not positive_live:
                for i,v in enumerate(new):
                    if v < 0:
                        new[i]=arb(0)
                    elif not v>0 and not v.is_zero():
                        raise ValueError('Live clamp ambiguous')
            x=arb_mat([[v] for v in new])

        target = count_or_none(r['target_live_cells'])
        selected,_ = select(states,target)
        results[(cohort,seg)]=(states,target)
        branches[(cohort,seg)]=(selected,parent_day,fallback)
        x,L = states[selected]
        if L <= 0:
            raise ValueError('Nonpositive endpoint live total')
        for pid in r['data_ids'].split(';'):
            if PASSAGES[pid]['recorded_ols_growth_per_day']!='NA':
                minimum = exact('1e-12')/CAPACITY
                assert L>minimum and size>minimum
                values[('growth',pid,'')]=(L/size).log()/exact(selected)
            if int(PASSAGES[pid]['karyotype_cell_count']):
                frac=x/L
                k_delta=correction*frac
                latent[pid]=frac
                corrections[pid]=k_delta
                for k,n in enumerate(GRID):
                    values[('kary',pid,str(n))]=frac[k,0]+k_delta[k,0]
    assert len(values)==114+12*133
    mask = set().union(*(x[2] for x in operators.values()))
    return values, branches, mask, latent, corrections


def bounded(p):
    return all(
        (exact(lo) < p[k] or (p[k]-exact(lo)).contains(0))
        and (p[k] < exact(hi) or (p[k]-exact(hi)).contains(0))
        for k,(lo,hi) in BOUNDS.items() if k in p
    )


def controls():
    ctx.prec=384
    points=parameter_points('control_parameters.tsv')
    p={k:exact(v) for k,v in points['anchor'].items()}
    O=exact('20.5')
    q=p.copy()
    q['O2_crit'] *= exact('0.99')
    q['n_O'] += exact('0.01')
    scale=hill(p,O)/hill(q,O)
    q['alpha_o2'] *= scale
    q['mu_hp'] *= scale
    assert bounded(q) and scale>1
    test=forward(p,'20.5')
    assert len(test[1])==122
    source=[dict(control='single_O_20.5',status='exact_symbolic',oxygen='20.5',
         feasible_branches=len(test[1]),ratio_H_over_Hprime=str(scale),
         explanation='Both alpha*h and mu*h are preserved identically; same states and branch selectors by induction')]

    face_a={k:exact(v) for k,v in points['survival_face_a'].items()}
    face_b={k:exact(v) for k,v in points['survival_face_b'].items()}
    assert face_a['buffer_smax'].is_zero() and face_b['buffer_smax'].is_zero()
    assert {k for k in COLS if points['survival_face_a'][k]!=points['survival_face_b'][k]}=={
         'buffer_beta','buffer_n_exp'}
    f=forward(face_a)
    source.append(dict(control='actual_multi_O_zero_smax_face',status='exact_symbolic',
         oxygen='actual',feasible_branches=len(f[1]),ratio_H_over_Hprime='NA',
         explanation='m>0 survival identically zero and m=0 survival identically one'))

    # Exact real-number distinction of the two widths despite identical binary64 output.
    w={k:exact(v) for k,v in points['width_plateau'].items()}
    a=forward(p)
    pid='SUM-159_NLS_2N_A5_seed'
    bin45=45-GRID[0]
    width_delta=(kernel(w['sigma_kary'])*a[3][pid])[bin45,0]-a[4][pid][bin45,0]
    assert nonzero(width_delta)
    source.append(dict(control='actual_multi_O_width_plateau',status='exact_binary64_only',
         oxygen='actual',feasible_branches=len(a[1]),ratio_H_over_Hprime='NA',
         explanation='Real difference at 2N A5 chromosome 45: '+str(width_delta)))
    delta_actual=forward(q)[0]
    gmax=max((abs(delta_actual[key]-a[0][key]).abs_lower() for key in a[0]
              if key[0]=='growth'))
    assert gmax>0
    source.append(dict(control='single_O_transform_applied_to_actual_multi_O',status='strictly_distinct',
         oxygen='actual',feasible_branches=122,ratio_H_over_Hprime=str(scale),
         explanation='Certified lower bound for largest growth gap: '+str(gmax)))
    table(OUT/'rigorous_controls.tsv',source)
    for row in source:
        print(row['control'],row['status'],flush=True)


def one_point(inp):
    ident,params=inp
    ctx.prec=384
    try:
        value=forward({k:exact(v) for k,v in params.items()})
        # Serialise outward enclosures using exact dyadic representations.
        return ident,({k:encode(v) for k,v in value[0].items()},value[1],value[2]),None
    except (ValueError,AssertionError) as exc:
        return ident,None,str(exc)


def covers(a,b):
    return a.overlaps(b)


def pair_review(workers):
    ctx.prec=384
    pairs=read(OUT/'selected_pairs.tsv')
    from algebra import decode
    points=parameter_points('selected_pair_parameters.tsv')
    expected=read(OUT/'paired_native_branches.tsv')
    for prefix in ('additional','targeted'):
        if (OUT/f'{prefix}_pairs.tsv').exists():
            pairs += read(OUT/f'{prefix}_pairs.tsv')
            points.update(parameter_points(f'{prefix}_pair_parameters.tsv'))
            expected += read(OUT/f'{prefix}_native_branches.tsv')
    native={}
    for row in expected:
        native.setdefault(row['id'],{})[(row['cohort'],row['segment'])]=(
            row['endpoint_day'],row['parent_day'],row['inoculum_fallback']=='TRUE')
    native_value_rows=read(OUT/'paired_native_outputs.tsv')
    for prefix in ('additional','targeted'):
        if (OUT/f'{prefix}_native_outputs.tsv').exists():
            native_value_rows += read(OUT/f'{prefix}_native_outputs.tsv')
    native_values={}
    for row in native_value_rows:
        native_values.setdefault(row['id'],{})[(row['stream'],row['passage'],row['chr'])]=float(row['value'])
    selected_ids={row[k] for row in pairs for k in ('left','right')}
    stats={}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for name,payload,error in pool.map(one_point,((k,points[k]) for k in sorted(selected_ids))):
            if payload is None:
                stats[name]=(None,error)
                print(name,'UNRESOLVED',error,flush=True)
                continue
            numbers,branches,mask=payload
            vals={k:decode(v) for k,v in numbers.items()}
            state_ok=branches==native[name]
            gap=max(abs(float(vals[k].mid())-native_values[name][k]) for k in vals)
            stats[name]=((vals,branches,mask,gap),None)
            print(name,'replayed',state_ok,'native output max gap',gap,flush=True)
            assert state_ok and gap<1e-10
    results=[]
    for r in pairs:
        left,err_a=stats[r['left']]
        right,err_b=stats[r['right']]
        if left is None or right is None:
            results.append(dict(category=r['category'],kind=r['kind'],left=r['left'],right=r['right'],
                result='unresolved',distinguishing_stream='',passage='',chromosome='',
                signed_right_minus_left='',growth_gap_lower='',growth_gap_upper='',kary_tv_lower='',
                kary_tv_upper='',endpoint_switch='',mask_switch='',fallback_switch='',reason=err_a or err_b))
            continue
        a,A,am,_=left; b,B,bm,_=right
        delta={k:b[k]-a[k] for k in a}
        strongest=max((k for k in delta if k[0]=='growth'),key=lambda k:abs(delta[k]).abs_lower())
        if not nonzero(delta[strongest]):
            strongest=max(delta,key=lambda k:abs(delta[k]).abs_lower())
        growth=[abs(v) for k,v in delta.items() if k[0]=='growth']
        all_tvs=[]
        for pid in {k[1] for k in delta if k[0]=='kary'}:
            all_tvs.append(sum(abs(delta['kary',pid,str(n)]) for n in GRID)/2)
        e=sum(A[k][0]!=B[k][0] for k in A)
        z=sum(A[k][2]!=B[k][2] for k in A)
        q=len(am^bm)
        assert e==int(r['endpoint']) and z==int(r['reset_mode']) and q==int(r['mask'])
        results.append(dict(category=r['category'],kind=r['kind'],left=r['left'],right=r['right'],
          result='certified_distinct' if nonzero(delta[strongest]) else 'unresolved',
          distinguishing_stream=strongest[0],passage=strongest[1],chromosome=strongest[2],
          signed_right_minus_left=str(delta[strongest]),
          growth_gap_lower=str(max(v.abs_lower() for v in growth)),
          growth_gap_upper=str(max(v.abs_upper() for v in growth)),
          kary_tv_lower=str(max(v.abs_lower() for v in all_tvs)),
          kary_tv_upper=str(max(v.abs_upper() for v in all_tvs)),
          endpoint_switch=e,mask_switch=q,fallback_switch=z,reason=''))
    table(OUT/'rigorous_pairs.tsv',results)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('phase',choices=('controls','pairs'))
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    if args.phase=='controls':
        controls()
    else:
        pair_review(args.workers)
