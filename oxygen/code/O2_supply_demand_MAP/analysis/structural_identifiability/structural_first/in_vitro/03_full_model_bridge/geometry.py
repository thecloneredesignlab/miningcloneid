"""Finite full-D secants; no optimizer or structural-rank claim."""
from pathlib import Path
import sys
from common import read,table
from runtime import boot

def main():
    study=Path(sys.argv[1]);task=sys.argv[2]
    sp,flint,np,out=boot(study,task)
    assert task=='G01'
    def vector(path):
        rows=read(path);return [r['output'] for r in rows],np.array([float(r['value']) for r in rows])
    ids,base=vector(study/'verification/R01/mean_vector.tsv')
    noise=float(read(study/'verification/R02/anchor_replay_checks.tsv')[0]['mean_vector_discrepancy'])
    points=read(study/'verification/M01/bridge_points.tsv');branches={r['eval_id']:r for r in read(study/'verification/R03/selection_branch_audit.tsv')}
    b=read(study/'verification/R01/branch.tsv')[0]
    names=[];matrix=[];diagnostics=[]
    for name in [r['parameter'] for r in read(study/'inputs/parameter_roles.tsv') if r['parameter']!='sigma_growth']:
        samples={};tags=[]
        if name=='sigma_kary':
            for i,r in enumerate(read(study/'verification/M01/observer_points.tsv'),1):
                path=study/'verification/O02'/('sigma_'+str(i)+'_vector.tsv')
                if path.exists():
                    ii,v=vector(path);assert ids==ii;samples[float(r['offset'])]=v
            tags=['observer_only; propagation unchanged']
        else:
            for r in points:
                if r['parameter']!=name or r['eval_id'] not in branches:continue
                ii,v=vector(study/'verification/R03'/(r['eval_id']+'_vector.tsv'));assert ii==ids
                samples[float(r['offset'])]=v;branch=branches[r['eval_id']]
                for key in ('selection_sha','pruning_sha','initial_zero_sha','score_floor_sha'):
                    if branch[key]!=b[key]:tags.append(key+'_switch')
                if branch['clamp_branch_certified']!='TRUE':tags.append('clamp_branch_not_certified')
        derivatives=[];signals=[]
        for h in (0.01,0.005):
            if h in samples and -h in samples:d=(samples[h]-samples[-h])/(2*h)
            elif h in samples:d=(samples[h]-base)/h
            elif -h in samples:d=(samples[-h]-base)/(-h)
            else:continue
            derivatives.append(d)
        for value in samples.values():signals.append(float(np.linalg.norm(value-base)))
        signal=max(signals,default=0.)
        agreement=float(np.linalg.norm(derivatives[-1]-derivatives[0])/max(np.linalg.norm(derivatives[-1]),np.linalg.norm(derivatives[0]),1e-300)) if len(derivatives)==2 else float('inf')
        d=derivatives[-1] if derivatives else np.zeros(len(base))
        branch_switch=any('switch' in t or 'not_certified' in t for t in tags)
        status='NUMERICALLY_FLAT_AT_REGISTERED_PROBES' if signal<=10*noise else 'STABLE_LOCAL_DIAGNOSTIC' if len(derivatives)==2 and agreement<=.1 and not branch_switch else 'PIECEWISE_OR_SCALE_UNRESOLVED_SECANT'
        diagnostics.append(dict(parameter=name,samples=len(samples),signal_l2=signal,repeat_noise_max=noise,two_scale_relative_gap=agreement,observed_branch_flags=';'.join(sorted(set(tags))),status=status,structural_classification='not_inferred_from_numeric_rank'))
        names.append(name);matrix.append(d)
    M=np.column_stack(matrix);assert M.shape==(len(ids),19) and np.isfinite(M).all()
    norms=np.linalg.norm(M,axis=0);scales=np.array([1/v if v>0 else 1 for v in norms]);E=M*scales
    singular=[];directions=[]
    for kind,A in (('registered_block_units',M),('declared_column_equilibration',E)):
        U,S,V=np.linalg.svd(A,full_matrices=False)
        for i,value in enumerate(S):singular.append(dict(scaling=kind,index=i,singular_value=float(value),relative_value=float(value/S[0]) if S[0]>0 else 0,nominal_rank_1e6=int(np.sum(S>S[0]*1e-6)),nominal_rank_1e8=int(np.sum(S>S[0]*1e-8)),nominal_rank_1e10=int(np.sum(S>S[0]*1e-10)),interpretation='finite secant geometry only; branch/scale diagnostics take precedence'))
        for direction in range(max(0,len(S)-3),len(S)):
            for j,name in enumerate(names):directions.append(dict(scaling=kind,direction=direction,parameter=name,coefficient=float(V[direction,j]),column_multiplier=float(scales[j]) if kind=='declared_column_equilibration' else 1,scope='nominal weak secant direction; not a global identifiable combination'))
    table(out/'local_observation_geometry.tsv',diagnostics);table(out/'singular_values.tsv',singular);table(out/'weak_secant_directions.tsv',directions)
    table(out/'column_scaling.tsv',[dict(parameter=n,unscaled_column_norm=float(norms[i]),equilibration_multiplier=float(scales[i]),parameter_unit='normalized declared transformed coordinate') for i,n in enumerate(names)])
    np.save(out/'secant_matrix.npy',M,allow_pickle=False)
    (out/'clamp_bound_argument.md').write_text('''# Sufficient clamp-branch bounds

For valid probabilities, G has nonnegative off-diagonals and G_NN>=-lambda_N.
With c in[0,1], live Euler steps are nonnegative if
1-DT*(lam_max+mu_hp*(N_MAX/44)^gamma_mu)>0; retained steps are nonnegative if
1-DT*k_clear>0. Both margins are recorded at every evaluated vector.

Between observations, total live growth is bounded above by lam_max*L. Thus
Euler growth over a gap is bounded by exp(lam_max*gap). If the maximum observed
L times that factor is below K, logistic crowding stays below capacity within
the gaps. Every segment includes its endpoints. This is a sufficient bound,
not a new simulation. If a bound fails, clamp branches remain uncertified;
unchanged daily zero masks alone would not prove absence of intervening clamps.
Floating underflow and initial-zero/floor masks remain separate flags.
''')
    table(out/'gate.tsv',[dict(status='PASS_ENGINEERING_DIAGNOSTICS',columns=19,mean_outputs=len(ids),structural_rank_claim=False,optimizer_calls=0)])
    table(out/'application_result.tsv',[dict(status='PASS',task=task,scope='registered secants; no fit')])
    print('FULL_D_FINITE_GEOMETRY_RETURNED',len(diagnostics),flush=True)

if __name__=='__main__':main()
