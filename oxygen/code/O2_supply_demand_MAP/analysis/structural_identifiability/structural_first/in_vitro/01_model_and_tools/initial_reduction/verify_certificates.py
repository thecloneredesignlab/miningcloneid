"""Independent determinant calculation, without importing model or AD code."""
from pathlib import Path
import csv,sys
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def determinant(A,p):
 A=[r[:] for r in A];ans=1
 for k in range(len(A)):
  pivot=next((i for i in range(k,len(A)) if A[i][k]%p),None)
  if pivot is None:return 0
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];ans=-ans
  value=A[k][k]%p;ans=ans*value%p;inverse=pow(value,-1,p)
  for i in range(k+1,len(A)):
   scale=A[i][k]*inverse%p
   for j in range(k+1,len(A)):A[i][j]=(A[i][j]-scale*A[k][j])%p
   A[i][k]=0
 return ans%p
if __name__=='__main__':
 root=Path(sys.argv[1]);receipts=[]
 for p in sorted(root.glob('certificates/*/receipt.tsv')):
  r={x['key']:x['value'] for x in read(p)};prime=int(r['prime']);n=int(r['rank']);rows=read(p.parent/'minor.tsv');A=[[int(x['c'+str(j)]) for j in range(n)] for x in rows]
  d=determinant(A,prime);assert d!=0 and d==int(r['nonzero_pivot_product']),(p,d,r['nonzero_pivot_product']);receipts.append((str(p.parent.name),n,prime,d))
  for h in p.parent.glob('hankel_*.tsv'):
   H=[[int(x['c'+str(j)]) for j in range(5)] for x in read(h)];hd=determinant(H,prime);assert hd!=0;receipts.append((p.parent.name+'/'+h.name,5,prime,hd))
 assert receipts,'No project certificate'
 with (root/'independent_minor_verification.tsv').open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['certificate','rank','prime','determinant']);w.writerows(receipts)
 print('INDEPENDENT_NONZERO_MINORS',len(receipts))
