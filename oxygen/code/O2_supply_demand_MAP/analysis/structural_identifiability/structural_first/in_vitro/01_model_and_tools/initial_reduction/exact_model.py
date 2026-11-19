"""Factorized rational five-class recurrent model over exact prime fields.

No finite differences or floating point enter a certificate. A nonzero minor
modulo a prime proves the corresponding rational-function minor is not zero.
The mathematical upper bound and output-field argument are separate evidence.
"""
from fractions import Fraction
from math import comb,factorial
NODES=(22,44,66,88,110)
PARAMETERS=('bmax','p0','a','k','w','alpha','dH','K','s2','s4')
POINT=dict(zip(PARAMETERS,('3/4','1/10000','1/100','1/50','1/1000','2/5','1/50','4/5','3/10','4/5')))
PRIMES=(2147483647,2147483629)
class Field:
 def __init__(self,p,n):self.p=p;self.n=n
 def __call__(self,value=0,index=None):
  if isinstance(value,Dual):return value
  q=Fraction(value);v=q.numerator*pow(q.denominator,-1,self.p)%self.p
  return Dual(self,v,tuple(int(i==index) for i in range(self.n)))
class Dual:
 def __init__(self,field,value,gradient):self.f=field;self.v=value%field.p;self.d=gradient
 def __add__(self,other):
  o=self.f(other);p=self.f.p;return Dual(self.f,self.v+o.v,tuple((a+b)%p for a,b in zip(self.d,o.d)))
 __radd__=__add__
 def __neg__(self):return Dual(self.f,-self.v,tuple(-a%self.f.p for a in self.d))
 def __sub__(self,other):return self+-self.f(other)
 def __rsub__(self,other):return self.f(other)+-self
 def __mul__(self,other):
  o=self.f(other);p=self.f.p;return Dual(self.f,self.v*o.v,tuple((a*o.v+b*self.v)%p for a,b in zip(self.d,o.d)))
 __rmul__=__mul__
 def __truediv__(self,other):
  o=self.f(other);p=self.f.p;iv=pow(o.v,-1,p);v=self.v*iv%p
  return Dual(self.f,v,tuple((a-v*b)*iv%p for a,b in zip(self.d,o.d)))
 def __rtruediv__(self,other):return self.f(other)/self
 def __pow__(self,n):
  assert isinstance(n,int) and n>=0
  if n==0:return self.f(1)
  p=self.f.p;factor=n*pow(self.v,n-1,p)%p
  return Dual(self.f,pow(self.v,n,p),tuple(a*factor%p for a in self.d))
def dest(n):
 assert n>0
 return 0 if n<=32 else 1 if n<=54 else 2 if n<=76 else 3 if n<=98 else 4
def active(variant):return [n for n in PARAMETERS if not(variant=='a_zero' and n in ('a','k')) and not(variant=='flat' and n=='s4')]
def build(field,theta,oxygen):
 z=field(0);G=[[z for j in range(7)] for i in range(7)];M=[[z for j in range(7)] for i in range(7)]
 H=theta['K']/(theta['K']+field(oxygen));error=[];ordinary=[];retained=[];direct=[];ledger=[]
 for j,N in enumerate(NODES):
  ratio=field(Fraction(N,44));b=theta['bmax']/(1+theta['alpha']*H*ratio);d=theta['dH']*H*ratio
  p=theta['p0']+theta['a']*d/(theta['k']+d)
  weight=field(max(0,min(1,Fraction(N-44,44))));s=(1-weight)*theta['s2']+weight*theta['s4']
  B=[z for i in range(5)];viable=z;altered=z
  for m in range(N+1):
   prob=comb(N,m)*p**m*(1-p)**(N-m);mass=prob*s**m
   for daughter in (N-m,N+m):
    if daughter>0:
     B[dest(daughter)]=B[dest(daughter)]+mass;viable=viable+mass
     if m:altered=altered+mass
  dead=2-viable
  for i in range(5):G[i][j]=b*(1-theta['w'])*B[i]
  G[dest(2*N)][j]=G[dest(2*N)][j]+b*theta['w']
  G[j][j]=G[j][j]-b-d
  # Constant cell volume is the frozen production convention. Q is measured in
  # one fixed cell-volume equivalent; it is not chromosome-mass conservation.
  G[5][j]=d;G[6][j]=b*(1-theta['w'])*dead
  for i in range(5):M[i][j]=G[i][j]
  M[j][j]=M[j][j]+d;M[6][j]=G[6][j]
  error.append(b*(1-theta['w'])*(1-(1-p)**N));ordinary.append(b*(1-theta['w']));retained.append(b*(1-theta['w'])*altered);direct.append(d)
  ledger.append((N,b,d,p,s,viable,dead))
 G[5][5]=-field('1/10000');G[6][6]=-field('1/10000')
 return G,dict(error=error,ordinary=ordinary,retained=retained,direct=direct,division_matrix=M),ledger
def matvec(A,x):return [sum((a*b for a,b in zip(row,x)),x[0].f(0)) for row in A]
def series_div(a,b):
 result=[]
 for j in range(len(a)):result.append((a[j]-sum((b[k]*result[j-k] for k in range(1,j+1)),a[0].f(0)))/b[0])
 return result
def output_jets(G,flux,x0,kind,order=16,crowding=True):
 x=[x0];mx=[matvec(flux['division_matrix'],x0)];L=[sum(x0[:5],x0[0].f(0))];capacity=x0[0].f(10**12)
 for j in range(1,order+1):
  rhs=matvec(G,x[-1])
  if crowding:rhs=[v-sum((L[k]*mx[j-1-k][i] for k in range(j)),v.f(0))/capacity for i,v in enumerate(rhs)]
  x.append([v/j for v in rhs]);mx.append(matvec(flux['division_matrix'],x[-1]));L.append(sum(x[-1][:5],x0[0].f(0)))
 live=[sum(row[:5],row[0].f(0)) for row in x];outputs=[];labels=[]
 if 'L' in kind:outputs.append(live);labels.append('live')
 if 'C' in kind:
  for j in range(4):outputs.append(series_div([row[j] for row in x],live));labels.append('fraction_'+str(j+1))
 if 'Q' in kind:outputs.append([row[5]+row[6] for row in x]);labels.append('retained_dead')
 if 'E' in kind:
  for key in ('error','ordinary'):
   series=[sum((flux[key][j]*row[j] for j in range(5)),row[0].f(0)) for row in x]
   if crowding:series=[v-sum((L[k]*series[h-k] for k in range(h+1)),v.f(0))/capacity for h,v in enumerate(series)]
   outputs.append(series);labels.append(key+'_flux')
 return [(name+'[t^'+str(k)+']',value) for name,series in zip(labels,outputs) for k,value in enumerate(series)]
def rank_pivots(matrix,p):
 if not matrix:return 0,[],[],1
 a=[row[:] for row in matrix];m=len(a);n=len(a[0]);row_ids=list(range(m));pivot_rows=[];pivot_cols=[];r=0;det=1
 for c in range(n):
  found=next((i for i in range(r,m) if a[i][c]%p),None)
  if found is None:continue
  a[r],a[found]=a[found],a[r];row_ids[r],row_ids[found]=row_ids[found],row_ids[r]
  v=a[r][c]%p;det=det*v%p;inv=pow(v,-1,p)
  a[r]=[q*inv%p for q in a[r]]
  for i in range(r+1,m):
   factor=a[i][c]%p
   if factor:a[i]=[(u-factor*v)%p for u,v in zip(a[i],a[r])]
  pivot_rows.append(row_ids[r]);pivot_cols.append(c);r+=1
  if r==m:break
 return r,pivot_rows,pivot_cols,det
def cases():
 result=[];i=0
 for multi in (False,True):
  for unknown in (False,True):
   for kind in ('L','C','LC'):
    i+=1;result.append(dict(case='C%02d'%i,kind=kind,oxygens=('1/5','1','41/2') if multi else ('1/5',),unknown=unknown,cohorts=2,initial='generic_positive',variant='reference'))
 result.extend([
 dict(case='C13',kind='L',oxygens=('1/5',),unknown=False,cohorts=2,initial='pure',variant='reference'),
 dict(case='C14',kind='L',oxygens=('1/5',),unknown=True,cohorts=1,initial='generic_positive',variant='reference'),
 dict(case='C17',kind='LCQ',oxygens=('1/5',),unknown=True,cohorts=2,initial='generic_positive',variant='reference'),
 dict(case='C18',kind='LQ',oxygens=('1/5',),unknown=True,cohorts=2,initial='generic_positive',variant='reference'),
 dict(case='C19',kind='LCE',oxygens=('1/5',),unknown=True,cohorts=2,initial='generic_positive',variant='reference'),
 dict(case='C21',kind='LC',oxygens=('1/5',),unknown=False,cohorts=2,initial='generic_positive',variant='a_zero'),
 dict(case='C22',kind='LC',oxygens=('1/5','1','41/2'),unknown=True,cohorts=2,initial='generic_positive',variant='flat')])
 for case in result:case['crowding']=True;case['model']='R5D2-C'
 controls=[]
 for origin in ('C04','C05','C12'):
  old=next(c for c in result if c['case']==origin);controls.append(dict(old,case='L_'+origin,crowding=False,model='R5D2-L'))
 return result+controls
def evaluate_case(case,prime):
 names=active(case['variant']);experiments=len(case['oxygens'])*case['cohorts'];nparam=len(names);nvar=nparam+(5*experiments if case['unknown'] else 0)
 f=Field(prime,nvar);theta={name:f(POINT[name],names.index(name) if name in names else None) for name in PARAMETERS}
 if case['variant']=='a_zero':theta['a']=f(0);theta['k']=f(POINT['k'])
 if case['variant']=='flat':theta['s4']=theta['s2']
 rows=[];labels=[];columns=names[:];exp=0;hankels=[]
 for oxygen in case['oxygens']:
  G,flux,ledger=build(f,theta,oxygen)
  for cohort in range(case['cohorts']):
   initial=([0,1,0,0,0] if cohort==0 else [0,0,0,1,0]) if case['initial']=='pure' else ([1,8,2,1,1] if cohort==0 else [1,1,2,8,1])
   values=[]
   for j,value in enumerate(initial):
    index=nparam+5*exp+j if case['unknown'] else None;values.append(f(value,index))
    if case['unknown']:columns.append('initial_'+str(exp)+'_'+str(j))
   values.extend([f(0),f(0)])
   outputs=output_jets(G,flux,values,case['kind'],crowding=case['crowding'])
   for label,value in outputs:rows.append(list(value.d));labels.append('O='+oxygen+';cohort='+str(cohort)+';'+label)
   if case['kind']=='L' and not case['crowding']:
    sequence=[value.v for label,value in outputs if label.startswith('live[')];moments=[v*factorial(k)%prime for k,v in enumerate(sequence)]
    hankels.append([[moments[i+j] for j in range(5)] for i in range(5)])
   exp+=1
 rank,rr,cc,det=rank_pivots(rows,prime);nuisance=rank_pivots([r[nparam:] for r in rows],prime)[0] if case['unknown'] else 0
 symmetry=1 if len(case['oxygens'])==1 else 0
 scale=experiments if case['unknown'] and case['kind']=='C' and not case['crowding'] else 0
 upper=nvar-symmetry-scale
 if case['kind']=='L' and case['unknown'] and not case['crowding']:upper=min(upper,5*experiments+5*len(case['oxygens']))
 assert rank<=upper,'Mathematical upper bound inconsistent with exact witness'
 augmented={}
 for j,name in enumerate(names):
  unit=[int(k==j) for k in range(nvar)];augmented[name]=rank_pivots(rows+[unit],prime)[0]
 return dict(rank=rank,nuisance_rank=nuisance,parameter_rank=rank-nuisance,upper=upper,columns=columns,row_labels=labels,jacobian=rows,pivot_rows=rr,pivot_cols=cc,nonzero_pivot_product=det,augmented=augmented,hankels=hankels)
