"""Small arithmetic and finite-allocation contract tests, no project analysis."""
import unittest
from fractions import Fraction
import exact_model as m
import workflow as w
from verify_certificates import determinant
class Contracts(unittest.TestCase):
 def test_dual_rational_derivative(self):
  f=m.Field(101,2);x=f('2/3',0);y=f('3/7',1);a=(x*x+y)/(1+x*y)
  X=Fraction(2,3);Y=Fraction(3,7)
  ds=[(2*X*(1+X*Y)-(X*X+Y)*Y)/(1+X*Y)**2,((1+X*Y)-(X*X+Y)*X)/(1+X*Y)**2]
  self.assertEqual(a.d,tuple(f(v).v for v in ds))
 def test_rank_minor_independent(self):
  A=[[0,1,2],[2,3,5],[7,11,13],[1,0,0]];rank,rr,cc,d=m.rank_pivots(A,101)
  self.assertEqual(rank,3);self.assertEqual(d,determinant([[A[i][j] for j in cc] for i in rr],101))
 def test_exact_binomial_and_wgd_accounting(self):
  f=m.Field(2147483647,0);theta={n:f(m.POINT[n]) for n in m.PARAMETERS};G,flux,ledger=m.build(f,theta,'1/5')
  for j,(N,b,d,p,s,live,dead) in enumerate(ledger):
   self.assertEqual((live+dead).v,f(2).v)
   self.assertEqual((sum(G[i][j] for i in range(7))+b*theta['w']-b).v,0)
 def test_endpoint_inactive_map(self):
  self.assertNotIn('k',m.active('a_zero'));self.assertNotIn('a',m.active('a_zero'));self.assertNotIn('s4',m.active('flat'))
 def test_budget_live_and_failed(self):
  row=dict(phase='tool',state='FAILED',actual='108',reservation='7200');live=dict(phase='tool',state='RUNNING',actual='',reservation='7200')
  p={r['phase']:r for r in w.position([row,live])};self.assertEqual(p['tool']['actual_core_seconds'],108);self.assertEqual(p['tool']['committed_core_seconds'],7200)
 def test_overcommit_refused(self):
  with self.assertRaises(ValueError):w.position([dict(phase='tool',state='RUNNING',actual='',reservation=str(9*3600))])
 def test_nonlinear_taylor_recursion(self):
  f=m.Field(2147483647,2);r=f('2/3',0);a=f(7,1);G=[[f(0) for j in range(7)] for i in range(7)];G[0][0]=r
  x=[a]+[f(0)]*6;out=m.output_jets(G,{'division_matrix':G},x,'L',order=2,crowding=True)
  expected=r*r*a*(1-a/f(10**12))*(1-2*a/f(10**12))/2
  self.assertEqual(out[2][1].v,expected.v);self.assertEqual(out[2][1].d,expected.d)
 def test_fixed_oxygen_symmetry(self):
  f=m.Field(2147483647,0);theta={n:f(m.POINT[n]) for n in m.PARAMETERS};h=theta['K']/(theta['K']+f('1/5'));changed=dict(theta);changed['K']=f('7/5');hp=changed['K']/(changed['K']+f('1/5'))
  for n in ('alpha','dH'):changed[n]=theta[n]*h/hp
  A,flux,_=m.build(f,theta,'1/5');B,other,_=m.build(f,changed,'1/5')
  self.assertEqual([[v.v for v in row] for row in A],[[v.v for v in row] for row in B]);self.assertEqual([[v.v for v in row] for row in flux['division_matrix']],[[v.v for v in row] for row in other['division_matrix']])
 def test_scale_symmetry_only_linear_control(self):
  f=m.Field(2147483647,0);theta={n:f(m.POINT[n]) for n in m.PARAMETERS};A,flux,_=m.build(f,theta,'1/5');x=[f(v) for v in (1,8,2,1,1,0,0)];xx=[2*v for v in x]
  lin=[v.v for _,v in m.output_jets(A,flux,x,'C',order=1,crowding=False)];lin2=[v.v for _,v in m.output_jets(A,flux,xx,'C',order=1,crowding=False)];self.assertEqual(lin,lin2)
  non=[v.v for _,v in m.output_jets(A,flux,x,'C',order=1,crowding=True)];non2=[v.v for _,v in m.output_jets(A,flux,xx,'C',order=1,crowding=True)];self.assertNotEqual(non,non2)
if __name__=='__main__':unittest.main()
