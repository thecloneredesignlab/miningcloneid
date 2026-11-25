"""Certified derivatives on the experimental culture schedule."""

import argparse
from collections import Counter, defaultdict
import csv
from fractions import Fraction
from math import ceil
from pathlib import Path
import resource
import sys
from time import process_time, perf_counter

from flint import arb, arb_mat, ctx
import numpy as np
from scipy.linalg import qr

from paths import (
    ALIGNMENT as ALIGN, FINITE_STATE, FIT, FROZEN_DESIGN, IN_VITRO,
    STRUCTURAL as STRUCT, output_directory,
)

for directory in (
    FINITE_STATE,
    IN_VITRO / "04_production_schedule/local_rank",
    IN_VITRO / "05_experimental_passages/unknown_seeding",
    IN_VITRO / "05_experimental_passages/donor_inheritance",
):
    sys.path.insert(0, str(directory))
import inherit
from algebra import decode, encode, exact, initial

N = 133
M = 14
CHROMOSOMES = tuple(range(22, 155))
CAPACITY = 10**12
MECHANISMS = tuple(inherit.MECHANISMS)
OBSERVATION_COORDINATES = ("eta_H", "eta_M", "class_log_odds", "misseg_dead_share")
INITIALS = tuple(inherit.INITIALS)
DT = Fraction(1, 20)
COMMON_COLS = MECHANISMS + OBSERVATION_COORDINATES
FULL_CASES = ("C", "CD", "CK", "DK", "CDK")
MODALITY = {"corrected_count": "C", "death_fraction": "D", "karyotype": "K"}
ROOT_FILE = STRUCT / "20261118_culture_experimental_fixed_inoculum_01/counts_shared/mechanism_witness.tsv"


def read(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, records):
    records = list(records)
    if not records:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def ticks(value):
    return int(round(float(value) * 20))


def design(case, design_dir=FROZEN_DESIGN):
    graph = read(design_dir / "passage_graph.tsv")
    observations = defaultdict(list)
    for row in read(design_dir / "event_mask.tsv"):
        observations[row["passage"]].append(row)
    ages = {r["seed"]: Fraction(r["elapsed_age_days"]) for r in read(FIT / "a7k_event_contract.tsv")}
    roots, demand, wanted = {}, defaultdict(Counter), defaultdict(set)
    for row in graph:
        node, parent = row["passage_id"], row["donor_passage"]
        row["parent"] = parent
        if case in ("free", "elapsed_free", "free_seeds") and row["mixture_rule_free_mixtures"] == "free_simplex":
            row["parent"] = ""
        roots[node] = roots[row["parent"]] if row["parent"] else node
        row["root"] = roots[node]
        row["seed"] = int(float(row["fixed_seed_live_objects"]))
        row["dt_last"] = DT
        row["final"] = None
        if parent:
            row["donor_tick"] = ticks(row["donor_model_day"])
            wanted[parent].add(row["donor_tick"])
            if node not in ages:
                demand[parent][row["donor_tick"]] += row["seed"]
        if case.startswith("elapsed") and node in ages:
            age = ages[node]
            row["final"] = ceil(age / DT)
            row["dt_last"] = age - (row["final"] - 1) * DT
        for obs in observations[node]:
            obs["tick"] = ticks(obs["model_day"]) if row["final"] is None else row["final"]
            wanted[node].add(obs["tick"])
    assert sum(map(len, demand.values())) == 110
    return graph, observations, demand, wanted


def witness(case, name):
    raw = {r["parameter"]: r["exact_value"] for r in read(ROOT_FILE)}
    raw.update(eta_H="1/3", eta_M="2/3", class_log_odds="0", misseg_dead_share="1/4",
               init_mean_2N="47", init_sd_2N="2", init_mean_4N="88", init_sd_4N="4")
    for cohort in ("2N", "4N"):
        raw[f"tilt_mean_{cohort}"] = "0"
        raw[f"tilt_shape_{cohort}"] = "0"
    if name == "second":
        raw.update(p_misseg="11/1000", buffer_beta="3/2", class_log_odds="1/10",
                   eta_H="2/5", eta_M="3/5", misseg_dead_share="3/10")
    if name == "archived_fit":
        fit = read(FIT / "selected_primary.tsv")[0]
        for name in MECHANISMS + ("class_log_odds", "misseg_dead_share"):
            raw[name] = fit[name]
        raw["eta_H"], raw["eta_M"] = fit["eta_H_2N"], fit["eta_M_2N"]
    return raw


def proxy_states():
    counts = defaultdict(Counter)
    for r in read(ALIGN / "cell_alignment.tsv"):
        for cohort in ("2N", "4N"):
            if r["assay_origin_id"] == f"SUM-159_NLS_{cohort}_MRCA_harvest":
                counts[cohort][int(r["observed_N_22_autosomes"])] += 1
    assert all(sum(counts[c].values()) == 20 for c in ("2N", "4N"))
    return dict((c, [Fraction(counts[c][n], 20) for n in CHROMOSOMES]) for c in counts)


def free_profile(cohort):
    center = 46 if cohort == "2N" else 88
    weights = [Fraction(9**4, 9**4 + (n - center)**4) for n in CHROMOSOMES]
    return [w / sum(weights) for w in weights]


def rank_context():
    return arb_mat([[1]*N]), arb_mat([[int(i == j) for j in range(N)] for i in range(N)])


def normalize(x, tangent, ones):
    total = (ones*x)[0,0]
    assert total > 0
    pi = x/total
    return pi, (tangent - pi*(ones*tangent))/total


def init_state(row, cols, p, samples, proxies, ones, case):
    node, cohort, parent = row["passage_id"], row["cohort"], row["parent"]
    seed = exact(Fraction(row["seed"], CAPACITY))
    if parent:
        state, tangent = samples[parent][row["donor_tick"]]
        pi, jac = normalize(state, tangent, ones)
    elif case == "gaussian":
        pi, jac = initial(p, cohort, cols, grid=CHROMOSOMES, N0=1)
    else:
        values = free_profile(cohort) if case in ("free", "elapsed_free", "free_seeds") else proxies[cohort]
        pi = arb_mat([[exact(q)] for q in values])
        jac = arb_mat(N, len(cols))
        if case in ("free", "elapsed_free", "free_seeds"):
            for j in range(N-1):
                at = len(COMMON_COLS)+j
                jac[j,at], jac[N-1,at] = arb(1), arb(-1)
    if case == "tv" and "_A7K_" in node:
        offset = 44 if cohort == "2N" else 88
        z = [exact(Fraction(n-offset,44)) for n in CHROMOSOMES]
        # The two cohort bases match the archived wrapper, up to a constant that normalizes away.
        for j, basis in zip((f"tilt_mean_{cohort}", f"tilt_shape_{cohort}"), (z, [q*q for q in z])):
            dim = cols.index(j)
            expectation = sum(pi[i,0]*q for i,q in enumerate(basis))
            for i,q in enumerate(basis):
                jac[i,dim] += pi[i,0]*(q-expectation)
    h = m = arb(0)
    dh = arb_mat(1,len(cols))
    dm = arb_mat(1,len(cols))
    if case == "free_seeds" and row["donor_passage"]:
        # The coordinate is the viable fraction at seeding. This is an invertible scaling of L0.
        h = m = seed / 200
        seed *= exact("0.99")
        slot = cols.index("seed:" + node)
        scale = exact(Fraction(row["seed"], CAPACITY))
        for i in range(N):
            jac[i,slot] += pi[i,0]*scale/seed
        dh[0,slot] = dm[0,slot] = -scale/2
    return seed*pi, seed*jac, h, m, dh, dm


def with_params(x, m):
    """Place the mechanistic derivative block at the left of a wider row."""
    return arb_mat([[x[i,j] if j < x.ncols() else 0 for j in range(m)] for i in range(x.nrows())])


def observe(row, x, J, H, V, dH, dV, total, L, p, values, jac):
    label = row["event"]
    kind = MODALITY[row["modality"]]
    if kind == "C":
        c = L + H + V
        assert c*CAPACITY > 1
        values[kind,label] = (c*CAPACITY).mid()
        jac[kind,label] = ((total+dH+dV)/c).entries()
    elif kind == "D":
        eps = exact("1e-12")
        q = p["misseg_dead_share"]
        dead = H + q*V + eps
        live = L + (1-q)*V + eps
        assert dead > 0 and live > 0, (label, dead, live)
        g = p["class_log_odds"] + (dead/live).log()
        odds = g.exp()
        values[kind,label] = (eps+(1-2*eps)*odds/(1+odds)).mid()
        dg = (dH+q*dV)/dead - (total+(1-q)*dV)/live
        dg[0,M+2] += 1
        dg[0,M+3] += V/dead + V/live
        jac[kind,label] = dg.entries()
    else:
        pi = x/L
        jp = (J - pi*total)/L
        for i in range(N):
            values[kind,label,i+22] = ((1-exact("1e-12"))*pi[i,0]+exact("1e-12")/N).mid()
            if i < N-1:
                jac[kind,label,i+22] = [jp[i,j] for j in range(J.ncols())]


def propagate(case, name, wanted_roots=None, design_dir=FROZEN_DESIGN, *,
              parameter_values=None, supply_receipts=None):
    graph, obs, demand, wanted = design(case, design_dir)
    raw = witness(case, name) if parameter_values is None else parameter_values
    p = {key: exact(value) for key,value in raw.items()}
    proxies = proxy_states()
    ones, eye = rank_context()
    rates, samples, blocks, checks, predictions = {}, {}, {}, [], {}
    uses = Counter((r["parent"],r["donor_tick"]) for r in graph if r["parent"]
                   and (not wanted_roots or r["root"] in wanted_roots))
    for row in graph:
        node, root, o = row["passage_id"], row["root"], row["oxygen_pct"]
        if wanted_roots and root not in wanted_roots:
            continue
        if root not in blocks:
            local = COMMON_COLS
            if case in ("free","elapsed_free","free_seeds"):
                local += tuple(f"pi_{n}" for n in CHROMOSOMES[:-1])
            if case == "free_seeds":
                local += tuple("seed:"+r["passage_id"] for r in graph if r["root"]==root and r["donor_passage"])
            if case == "gaussian":
                local += INITIALS
            if case == "tv":
                local += tuple(f"tilt_{kind}_{c}" for c in ("2N","4N") for kind in ("mean","shape"))
            blocks[root] = (local, {})
        cols, jac = blocks[root]
        m = len(cols)
        if o not in rates:
            G, dG, diagonal, Hrate, dHrate, Vrate, dVrate, mask = inherit.operator(p,o)
            transport = arb_mat([[diagonal[i] if i==j else 0 for j in range(N)] for i in range(N)])
            rates[o] = G,dG,transport,Hrate,dHrate,Vrate,dVrate,mask
        G,dG,A,BH,dBH,BV,dBV,mask = rates[o]
        x,J,H,V,dH,dV = init_state(row,cols,p,samples,proxies,ones,case)
        if row["parent"]:
            reference = (row["parent"],row["donor_tick"])
            uses[reference] -= 1
            if not uses[reference]:
                del samples[reference[0]][reference[1]]
        assert wanted[node], node
        by_tick = defaultdict(list)
        for r in obs[node]:
            by_tick[r["tick"]].append(r)
        samples[node] = {}
        last = max(wanted[node])
        min_margin = None
        min_density_gap = None
        for t in range(last+1):
            L = (ones*x)[0,0]
            assert L > 0 and L < 1, (node,t,L)
            crowd = 1-L
            if min_density_gap is None or crowd.lower() < min_density_gap:
                min_density_gap = crowd.lower()
            if uses[node,t]:
                samples[node][t] = x,J
            if t in demand[node]:
                margin = L*CAPACITY-demand[node][t]
                if supply_receipts is not None:
                    supply_receipts[node, t] = margin
                assert margin > 0, (case,name,node,t,margin)
                if min_margin is None or margin.lower() < min_margin:
                    min_margin = margin.lower()
            total = ones*J
            for observation in by_tick[t]:
                observe(observation,x,J,H,V,dH,dV,total,L,p,predictions,jac)
            if t == last:
                break
            ratio = exact(row["dt_last"]/DT) if row["final"] is not None and t+1 == last else arb(1)
            decay = 1 - exact(DT)*ratio*p["k_clear"]
            hyp = (BH*x)[0,0]
            mis = (BV*x)[0,0]
            dHyp = (BH*J) + with_params((dBH*x).transpose(),m)
            dMis = (BV*J) + with_params((dBV*x).transpose(),m)
            eh, em = p["eta_H"], p["eta_M"]
            H = decay*H + ratio*eh*hyp
            V = decay*V + ratio*em*crowd*mis
            dH = decay*dH + ratio*eh*dHyp
            dH[0,M] += ratio*hyp
            dV = decay*dV + ratio*em*(crowd*dMis-mis*total)
            dV[0,M+1] += ratio*crowd*mis
            gx = G*x
            forcing = dG*x
            F = arb_mat([[crowd*forcing[i*M+j,0] - x[i,0]*dBH[j,i] if j<M else 0
                          for j in range(m)] for i in range(N)])
            transport = (1-ratio)*eye + ratio*(A+crowd*G)
            J = transport*J + ratio*(F - gx*total)
            x = transport*x
        checks.append(dict(passage=node, root=root, minimum_supply_margin_cells=str(min_margin) if min_margin is not None else "",
                           minimum_logistic_gap=str(min_density_gap), operator_mask=mask, steps=last))
        if len(checks)%15 == 0:
            print(f"{case} {name} {ctx.prec} bits: {len(checks)} passages",flush=True)
    return blocks, checks, predictions, raw


def entries_at(rows, ri, ci):
    return arb_mat([[rows[i][j] for j in ci] for i in ri])


def precondition(Mx):
    """Check an approximate dyadic inverse using outward arithmetic."""
    size = Mx.nrows()
    I = arb_mat([[int(i==j) for j in range(size)] for i in range(size)])
    try:
        R = Mx.mid().inv().mid()
        E = I - R*Mx
        norm = max(sum(abs(E[i,j]) for j in range(size)).upper() for i in range(size))
    except (ValueError, ZeroDivisionError):
        return None, arb("inf")
    return R, norm


def qr_take(rows, n, columns=False):
    values = np.array([[float(v.mid()) for v in row] for row in rows])
    if columns:
        values = values.T
    scale = np.max(np.abs(values),axis=0)
    scale[scale==0] = 1
    values /= scale
    scale = np.max(np.abs(values),axis=1)
    scale[scale==0] = 1
    values /= scale[:,None]
    return list(map(int, qr(values.T, mode="economic",pivoting=True,check_finite=True)[2][:n]))


def choose_and_check(rows, ri, ci):
    mat = entries_at(rows,ri,ci)
    pre, bound = precondition(mat)
    return mat, pre, bound


def matrix_rows(label, matrix, row_ids, columns):
    for i,event in enumerate(row_ids):
        for j,name in enumerate(columns):
            yield dict(certificate=label, row=i, column=j, output="|".join(map(str,event)), coordinate=name,
                       **encode(matrix[i,j]))


def rank_case(blocks, mask, case):
    global_count = len(COMMON_COLS) + (4 if case in ("gaussian","tv") else 0)
    global_cols = list(range(global_count))
    if "D" not in mask:
        global_cols.remove(M+2)
        global_cols.remove(M+3)
    effective, tags, certificates, slopes = [], [], [], []
    root_rank, unresolved = 0, []
    for root,(cols,J) in blocks.items():
        keys = [k for k in J if k[0] in mask]
        if not keys:
            continue
        rows = [J[k] for k in keys]
        local_cols = list(range(global_count,len(cols)))
        if not local_cols:
            effective += [[r[j] for j in global_cols] for r in rows]
            tags += keys
            continue
        nuisance = [[r[j] for j in local_cols] for r in rows]
        n = min(len(keys),len(local_cols))
        ri = list(range(len(keys))) if len(keys)==n else qr_take(nuisance,n)
        ci = list(range(len(local_cols))) if len(local_cols)==n else qr_take(nuisance,n,columns=True)
        mat, _, bound = choose_and_check(nuisance,ri,ci)
        if not bound < 1 and len(keys)>n:
            one_sample = [i for i,k in enumerate(keys) if len(k)==3 and k[1]==next(v[1] for v in keys if len(v)==3)]
            if len(one_sample)==n:
                ri = one_sample
                mat, _, bound = choose_and_check(nuisance,ri,ci)
        if not bound < 1:
            unresolved.append(root)
            continue
        certificates.append((root,mat,[keys[i] for i in ri],[cols[local_cols[j]] for j in ci],bound))
        root_rank += n
        if len(keys)==n:
            continue
        rem = [i for i in range(len(keys)) if i not in set(ri)]
        B = entries_at(rows,ri,global_cols)
        K = entries_at(nuisance,rem,ci)
        A = entries_at(rows,rem,global_cols)
        P = mat.solve(B)
        effective += (A - K*P).tolist()
        tags += [keys[i] for i in rem]
        slopes.append((root, [keys[i] for i in rem]))
    return global_cols,effective,tags,certificates,root_rank,unresolved,slopes


def verify_minor_rows(rows):
    matrix = arb_mat(max(int(r["row"]) for r in rows)+1, max(int(r["column"]) for r in rows)+1)
    for r in rows:
        matrix[int(r["row"]),int(r["column"])] = decode(r)
    _, bound = precondition(matrix)
    return matrix, bound


def count_compensation(blocks):
    """Follow one verified count-preserving direction in the longest 2N block."""
    root = "SUM-159_NLS_2N_dp_seed"
    cols, jac = blocks[root]
    count_keys = [key for key in jac if key[0]=="C"]
    local = list(range(len(COMMON_COLS),len(cols)))
    count = [[jac[key][j] for j in local] for key in count_keys]
    chosen = qr_take(count,len(count_keys),columns=True)
    pivot, _, bound = choose_and_check(count,range(len(count)),chosen)
    if not bound < 1:
        return []
    chosen = [local[j] for j in chosen]
    result = []
    for target in ("p_misseg","buffer_beta"):
        t = cols.index(target)
        dpi = pivot.solve(arb_mat([[-jac[key][t]] for key in count_keys]))
        scale = max(1.,max(float(abs(v).upper()) for v in dpi.entries()))
        for kind in ("C","D","K"):
            keys = [key for key in jac if key[0]==kind]
            restored = [(key,(jac[key][t]+sum(jac[key][j]*dpi[i,0] for i,j in enumerate(chosen)))/exact(scale))
                        for key in keys]
            maximum = max(restored,key=lambda kv:float(abs(kv[1]).lower()))
            result.append(dict(root=root,target=target,output_kind=kind,rows=len(keys),direction_normalization=scale,
                               largest_certified_effect="|".join(map(str,maximum[0])),
                               derivative=str(maximum[1]),all_enclosures_contain_zero=all(v.contains(0) for _,v in restored)))
    return result


def run(args):
    ctx.prec = args.bits
    start, wall = process_time(), perf_counter()
    blocks, check, prediction, raw = propagate(args.case,args.witness,args.root,args.design)
    out = args.out
    out.mkdir(parents=True,exist_ok=True)
    write(out/"witness.tsv",(dict(parameter=k,value=v) for k,v in raw.items()))
    write(out/"branch_checks.tsv",check)
    write(out/"predictions.tsv",(dict(output="|".join(map(str,key)),value=float(value)) for key,value in prediction.items()))
    propagation_cpu = process_time()-start
    if args.validate:
        from validate_joint import compare
        # Avoid importing this file twice when it is executing as the entry point.
        validation = compare(args.case,raw,args.root,blocks,prediction,args.design,all_local=bool(args.root))
        write(out/"derivative_validation.tsv",validation)
        failed = [row for row in validation if not row["passed"]]
        if failed:
            raise ArithmeticError(f"Independent derivative gates failed: {failed[:6]}")
    summary = []
    for mask in args.masks:
        global_cols, eff, keys, certs, root_rank, unresolved, _ = rank_case(blocks,mask,args.case)
        names = list(next(iter(blocks.values()))[0])
        col_ids = [names[j] for j in global_cols]
        n_global = len(global_cols)
        eff_rank = 0
        if eff:
            order = qr_take(eff,min(n_global,len(eff)))
            for n in range(min(n_global,len(eff)),0,-1):
                ci = list(range(n_global)) if n==n_global else qr_take([eff[i] for i in order[:n]],n,columns=True)
                mat, _, bound = choose_and_check(eff,order[:n],ci)
                if bound < 1:
                    eff_rank = n
                    certs.append(("global_after_root_elimination",mat,[keys[i] for i in order[:n]],
                                  [col_ids[i] for i in ci],bound))
                    break
        for label,mat,keylist,colnames,bound in certs:
            print(f"{args.case} {mask} {label} {mat.nrows()}: bound {float(bound):.3g}",flush=True)
        write(out/f"{mask}_minors.tsv",(r for label,mat,keylist,colnames,_ in certs for r in matrix_rows(label,mat,keylist,colnames)))
        write(out/f"{mask}_certificate_bounds.tsv",(dict(certificate=label,order=mat.nrows(),bound=str(bound))
                                                   for label,mat,_,_,bound in certs))
        write(out/f"{mask}_effective.tsv",(dict(row="|".join(map(str,key)),coordinate=col,**encode(v))
                                            for key,row in zip(keys,eff) for col,v in zip(col_ids,row)))
        n_rows = sum(k[0] in mask for _,D in blocks.values() for k in D)
        summary.append(dict(case=args.case,witness=args.witness,mask=mask,precision=args.bits,rows=n_rows,
                            global_columns=n_global,root_column_rank_lower=root_rank,
                            effective_global_rank_lower=eff_rank,global_rank_target=n_global,
                            root_blocks_without_full_pivot="|".join(unresolved),effective_row_count=len(eff),
                            propagation_cpu_seconds=propagation_cpu,total_cpu_seconds=process_time()-start,
                            wall_seconds=perf_counter()-wall,
                            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        print(summary[-1],flush=True)
    if args.case=="free" and "C" in args.masks and "SUM-159_NLS_2N_dp_seed" in blocks:
        write(out/"count_compensation.tsv",count_compensation(blocks))
    write(out/"summary.tsv",summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",choices=("free","inherit","gaussian","tv","elapsed_free","elapsed_inherit","free_seeds"),required=True)
    parser.add_argument("--witness",choices=("original","second","archived_fit"),default="original")
    parser.add_argument("--bits",type=int,default=384)
    parser.add_argument("--root",action="append")
    parser.add_argument("--validate",action="store_true")
    parser.add_argument("--masks",nargs="+",choices=FULL_CASES,default=list(FULL_CASES))
    parser.add_argument("--design",type=Path,default=FROZEN_DESIGN)
    parser.add_argument("--out",type=output_directory,required=True)
    run(parser.parse_args())
