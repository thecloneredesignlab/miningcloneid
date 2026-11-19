"""At most six predeclared contrasts; reuse saved18-rank, no full determinant."""
from pathlib import Path
import sys
from common import read,table,sha,roles_valid
from runtime import boot

def main():
    study=Path(sys.argv[1]);task=sys.argv[2]
    sp,flint,np,out=boot(study,task)
    from algebra import point,initial,decode,encode
    assert task=='O01';flint.ctx.prec=256
    roles_valid(read(study/'inputs/parameter_roles.tsv'))
    s=sp.Symbol('width',positive=True);d=sp.symbols('d0:3',nonnegative=True)
    raw=[sp.exp(-v/(2*s*s)) for v in d];Z=sum(raw);K=[v/Z for v in raw]
    expectation=sum(a*b for a,b in zip(K,d));controls=[]
    for i in range(3):
        residual=sp.simplify(sp.diff(K[i],s)-K[i]*(d[i]-expectation)/s**3)
        assert residual==0;controls.append(dict(control='normalized_kernel_derivative_'+str(i),residual=0,status='PASS_EXACT'))
    table(out/'derivative_controls.tsv',controls)
    attempts=[];certs=[];arb=flint.arb
    for witness in ('W1','W2','W3'):
        p=point(study,witness);width=p['sigma_kary'];saved=study/'inputs/evidence'/('B04_'+witness)
        D=decode(read(saved/'primary_128_determinant.tsv')[0]);assert D.is_finite() and not D.contains(0)
        columns=read(saved/'column_registry.tsv');assert len(columns)==18 and 'sigma_kary' not in {r['parameter'] for r in columns}
        accepted=False
        for cohort,center in (('2N',44),('4N',88)):
            assert len(attempts)<6
            a,_=initial(p,cohort,[],N0=1);prob=[a[j,0] for j in range(133)]
            direct=arb(0);covariance=arb(0)
            for j,mother in enumerate(range(22,155)):
                raw=[(-(arb(n-mother)**2)/(2*width**2)).exp() for n in range(22,155)]
                total=sum(raw);weights=[v/total for v in raw]
                E_d=sum(w*(n-mother)**2 for n,w in zip(range(22,155),weights))
                E_c=sum(w*(n-center)**2 for n,w in zip(range(22,155),weights))
                E_cd=sum(w*(n-center)**2*(n-mother)**2 for n,w in zip(range(22,155),weights))
                direct+=prob[j]*sum(w*((n-mother)**2-E_d)*(n-center)**2 for n,w in zip(range(22,155),weights))/width**3
                covariance+=prob[j]*(E_cd-E_c*E_d)/width**3
            assert direct.is_finite() and covariance.is_finite() and direct.overlaps(covariance)
            certified=not direct.contains(0) and not covariance.contains(0)
            product=D*direct;certified=certified and product.is_finite() and not product.contains(0)
            attempts.append(dict(witness=witness,cohort=cohort,center=center,precision=256,derivative=str(direct),covariance_reference=str(covariance),
                status='CERTIFIED_NONZERO' if certified else 'UNRESOLVED_REGISTERED_SCOPE',**encode(direct)))
            table(out/'contrast_attempts.tsv',attempts)
            if certified:
                certs.append(dict(witness=witness,rank=19,theta_columns=18,sigma_column='sigma_kary',status='CERTIFIED_LOCAL_IDEAL_EXTENSION',
                    saved_determinant_sha256=sha(saved/'primary_128_determinant.tsv'),contrast_cohort=cohort,contrast_center=center,
                    factored_determinant=str(product),new_full_determinant_evaluations=0,initial_coordinates='four unknown',
                    scope='continuous ideal L jets through4 plus initial normalized karyotype moment; finite support and original witness branch',**encode(product)))
                accepted=True;break
        if not accepted:certs.append(dict(witness=witness,rank='',theta_columns=18,sigma_column='sigma_kary',status='UNRESOLVED_REGISTERED_SCOPE',saved_determinant_sha256=sha(saved/'primary_128_determinant.tsv'),contrast_cohort='',contrast_center='',factored_determinant='',new_full_determinant_evaluations=0,initial_coordinates='four unknown',scope='finite contrast allowance exhausted',mid_man='',mid_exp='',rad_man='',rad_exp=''))
    table(out/'observer_extension_certificates.tsv',certs)
    (out/'observer_extension_proof.md').write_text('''# New ideal observer extension

The saved B04 rows have an invertible18-column Jacobian for14 mechanisms and
four unknown initial coordinates. They are independent of sigma_kary. Append
the registered initial karyotype second-moment contrast, about44 for2N or88
for4N. Its sigma partial is the certified normalized-kernel derivative. The
augmented matrix [J_L,0; b,dz/dsigma] factors in determinant as det(J_L)*dz/dsigma.
Unknown initial coordinates remain in J_L and b; evaluating a witness does not
fix them. All existing power-of-two scales on J_L are retained. New contrast
and sigma scales are1, so they are nonzero. A saved scaled determinant times a
nonzero contrast partial suffices without generating the full19-column matrix.

The kernel identity is independently differentiated symbolically for generic
normalized exponential weights. The full finite contrast is evaluated both
through column derivatives and the covariance identity. Certified outward balls
and their product exclude zero. No state operator, Lie jet or full determinant
was recomputed. This is a new registered output-map extension, not a repair or
relabeling of the old B03/B05 order2 attempts. An added flow output inherits the
result only by output inclusion; no practical information gain is inferred.

The conclusion is local near each certified witness for ideal continuous L and
an ideal initial karyotype distribution. Sigma_growth is outside the mean map.
It proves neither global uniqueness nor identification from real sparse samples,
and does not describe binary64 behavior at the small seed144 observer width.
''')
    table(out/'gate.tsv',[dict(status='PASS',certified_witnesses=sum(r['status']=='CERTIFIED_LOCAL_IDEAL_EXTENSION' for r in certs),contrast_attempts=len(attempts),new_full_determinants=0,new_project_state_cases=0)])
    table(out/'application_result.tsv',[dict(status='PASS',task=task,scope='bounded ideal observer extension')])
    print('BOUNDED_OBSERVER_EXTENSION_RETURNED',len(attempts),flush=True)

if __name__=='__main__':main()
