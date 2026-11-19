"""Finite candidate nomination and independently reconstructed minor evidence."""
from flint import arb,arb_mat,fmpq_mat
from algebra import encode,decode,dim,independent_minor_check,nonzero
from protocol import table

def exponent(v):
    if v.mid().is_zero():return None
    m,e=v.mid().man_exp();return abs(int(m)).bit_length()+int(e)

def scaled_rows(rows):
    n=len(rows[0]);cs=[]
    for j in range(n):cs.append(-max((exponent(r[j]) for r in rows if exponent(r[j]) is not None),default=0))
    first=[[v*arb((1,cs[j])) for j,v in enumerate(r)] for r in rows]
    rs=[-max((exponent(v) for v in r if exponent(v) is not None),default=0) for r in first]
    return [[v*arb((1,rs[i])) for v in r] for i,r in enumerate(first)],rs,cs

def pivot_columns(q):
    rref,rank=q.rref()
    return [next(j for j in range(rref.ncols()) if rref[i,j]!=0) for i in range(rank)]

def nominate(rows,labels,columns,np):
    scaled,rs,cs=scaled_rows(rows);n=len(columns)
    targets=[j for j,name in enumerate(columns) if name not in ('O2_crit','n_O')] if labels[0]['case']=='B01' else list(range(n))
    def choose(pool,cols):
        rational=fmpq_mat([[scaled[i][j].mid().fmpq() for j in cols] for i in pool])
        ids=pivot_columns(rational.transpose())
        if len(ids)<len(cols):return None
        return [pool[i] for i in ids[:len(cols)]]
    allrows=list(range(len(rows)));r1=choose(allrows,targets)
    result=[]
    if r1:result.append(dict(candidate='primary',rows=r1,columns=targets))
    nominal=np.array([[float(v.mid()) for v in r] for r in scaled],dtype=float)
    assert np.isfinite(nominal).all()
    sv=np.linalg.svd(nominal,compute_uv=False)
    rn=int(sum(v>1e-10*sv[0] for v in sv)) if len(sv) and sv[0]>0 else 0
    rn=min(rn,len(targets))
    pool=[]
    for i,l in enumerate(labels):
        if l['N']=='':pool.append(i)
        elif l['output']=='flow_density' and 1.5<=float(l['N'])<=6.5:pool.append(i)
        elif l['output']!='flow_density' and 33<=int(l['N'])<=120:pool.append(i)
    if rn>0:
        q=fmpq_mat([[scaled[i][j].mid().fmpq() for j in targets] for i in pool])
        cids=pivot_columns(q)[:rn];c2=[targets[j] for j in cids]
        r2=choose(pool,c2) if c2 else None
        if r2 and not any(r2==r['rows'] and c2==r['columns'] for r in result):result.append(dict(candidate='stratified',rows=r2,columns=c2))
    return result[:2],rs,cs,[float(v) for v in sv],rn

def submatrix(rows,selection,rs,cs):
    return arb_mat([[rows[i][j]*arb((1,rs[i]+cs[j])) for j in selection['columns']] for i in selection['rows']])

def save_matrix(path,M):
    n,m=dim(M)
    table(path,[dict(row=i,column=j,**encode(M[i,j])) for i in range(n) for j in range(m)])

def validate_and_save(out,selection,M,reference,precision):
    assert dim(M)==dim(reference)
    for i in range(M.nrows()):
        for j in range(M.ncols()):assert M[i,j].overlaps(reference[i,j]),'Independent selected derivative enclosures disagree'
    determinant=M.det()
    positive=nonzero(determinant)
    verified,bound=independent_minor_check(reference) if positive else (False,arb('inf'))
    record=dict(candidate=selection['candidate'],precision=precision,size=M.nrows(),nonzero_determinant=positive,
                independent_neumann_verified=verified,neumann_bound=str(bound),determinant=str(determinant),
                status='CERTIFIED' if positive and verified else 'UNRESOLVED')
    if positive and verified:
        prefix=selection['candidate']+'_'+str(precision)
        save_matrix(out/(prefix+'_matrix.tsv'),M);save_matrix(out/(prefix+'_reference.tsv'),reference)
        table(out/(prefix+'_determinant.tsv'),[dict(precision=precision,**encode(determinant))])
    return record
