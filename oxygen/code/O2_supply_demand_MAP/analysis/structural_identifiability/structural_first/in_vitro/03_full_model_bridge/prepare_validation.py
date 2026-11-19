"""Freeze added saved-reference checks before any perturbation propagation."""
from pathlib import Path
import ast
import shutil
import subprocess
import sys
from common import *

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();dest=study/'validation_release'
assert not dest.exists();dest.mkdir()
for name in ('common.py','policy.py','runtime.py'):shutil.copy2(study/'release'/name,dest/name)
r=(study/'release/sif_checks.R').read_text()
anchor="    b<-readRDS(file.path(study,'verification/M01/anchor_bundle.rds'));baseline<-readRDS(file.path(study,'verification/R01/result.rds'))"
assert r.count(anchor)==1
addition='''
    fit<-readRDS(file.path(study,'inputs/data/fit_result.rds'));reference<-list()
    for(cohort in c('2N','4N')) {
      actual<-baseline$runs[[cohort]]$segment_results;old<-fit$best_components[[paste0('run_',cohort)]]$segment_results
      check(paste0('production_segment_identity_',cohort),length(actual)==length(old))
      for(i in seq_along(actual)) {
        check(paste0('production_selection_',cohort,'_',i),isTRUE(all.equal(actual[[i]]$selection,old[[i]]$selection,tolerance=1e-12,check.attributes=TRUE)))
        for(name in intersect(names(actual[[i]]$sim),names(old[[i]]$sim))) {
          if(name%in%c('cache_g_build','cache_g_hit','cache_g_hysteresis'))next
          a<-actual[[i]]$sim[[name]];z<-old[[i]]$sim[[name]]
          if(!is.numeric(a)||!is.numeric(z))next
          stopifnot(identical(dim(a),dim(z)),identical(is.finite(a),is.finite(z)))
          finite<-is.finite(a);gap<-if(any(finite))max(abs(a[finite]-z[finite])/pmax(1,abs(z[finite])))else 0
          stopifnot(gap<=1e-10)
          reference[[length(reference)+1L]]<-data.frame(cohort=cohort,segment=i,quantity=name,max_scaled_gap=gap)
        }
      }
    }
    check('archived_production_objective',abs(baseline$score$objective-fit$best_components$objective)<=1e-8)
    sf_emit(do.call(rbind,reference),file.path(out,'production_reference_checks.tsv'))
    reset<-unlist(lapply(b$adapters,function(a)vapply(a$segments,function(seg){v<-as.numeric(seg$initial_cells);if(!is.finite(v)||v<=0)b$cfg$init_total_size else v},numeric(1))))
    domain<-list()
    domain_check<-function(p,id) {
      a<-b$cfg$DT*p$lam_max;positive<-1-b$cfg$DT*(p$lam_max+p$mu_hp*(b$cfg$N_MAX/44)^p$gamma_mu);stock<-1-b$cfg$DT*p$k_clear
      data.frame(eval_id=id,DT_lambda_max=a,live_positive_margin=positive,retained_positive_margin=stock,maximum_reset_over_K=max(reset)/b$cfg$K,
        certified=a<1&&positive>0&&stock>0&&max(reset)<b$cfg$K,scope='source real-arithmetic sufficient invariant region; not literal floating-bit proof')
    }
    domain[[1]]<-domain_check(b$p,'R01')
    registered<-sf_tab(file.path(study,'verification/M01/bridge_points.tsv'))
    for(i in seq_len(nrow(registered)))if(registered$status[i]=='registered') {p<-b$p;p[[registered$parameter[i]]]<-registered$new_value[i];domain[[length(domain)+1L]]<-domain_check(p,registered$eval_id[i])}
    sf_emit(do.call(rbind,domain),file.path(out,'clamp_invariance_domain.tsv'))
'''
r=r.replace(anchor,anchor+addition)
(dest/'sif_checks.R').write_text(r)
g=(study/'release/geometry.py').read_text()
g=g.replace("    names=[];matrix=[];diagnostics=[]", "    invariant={r['eval_id']:r['certified']=='TRUE' for r in read(study/'verification/O02/clamp_invariance_domain.tsv')}\n    y,a=sp.symbols('y a');assert sp.expand(1-(y+a*y*(1-y))-(1-y)*(1-a*y))==0\n    table(out/'clamp_polynomial_identity.tsv',[dict(residual=0,domain='0<=y<1;0<=a<1',identity='1-[y+a*y*(1-y)]=(1-y)*(1-a*y)',meaning='invariant region given nonnegative live steps and source total birth bound')])\n    names=[];matrix=[];diagnostics=[]")
g=g.replace("if branch['clamp_branch_certified']!='TRUE':tags.append('clamp_branch_not_certified')", "if branch['clamp_branch_certified']!='TRUE' and not invariant.get(r['eval_id'],False):tags.append('clamp_branch_not_certified')")
g=g.replace("status='NUMERICALLY_FLAT_AT_REGISTERED_PROBES' if signal<=10*noise", "status='UNRESOLVED_NO_FEASIBLE_PROBE' if not samples else 'NUMERICALLY_FLAT_AT_REGISTERED_PROBES' if signal<=10*noise")
g=g.replace("# Sufficient clamp-branch bounds", "# Sufficient clamp-branch bounds and invariant-region refinement")
g=g.replace("Between observations, total live growth", "An independent tighter source bound uses y=L/K and a=DT*lam_max<1. The\nidentity1-[y+a*y*(1-y)]=(1-y)*(1-a*y)>=0 preserves L<K at every\nnonnegative real-arithmetic step, provided every root/reseed total is below K.\nThis uses sum(G_column)<=lam_max and mu>=0. Its margins were frozen in O02\nbefore perturbation outcomes were inspected, with an exact polynomial check\nhere. It refines an inconclusive upper bound, not an observed clamp switch.\n\nBetween observations, total live growth")
(dest/'geometry.py').write_text(g)
s=(study/'release/runner.sh').read_text()
s=s.replace('"$study/release/verify_source.py"','"$study/validation_release/verify_source.py"')
s=s.replace('"$study/release/sif_checks.R"','"$study/validation_release/sif_checks.R"')
s=s.replace('"$study/release/$program"','"$study/validation_release/$program"')
s=s.replace('--bind "$study/release:$study/release:ro"','--bind "$study/release:$study/release:ro" --bind "$study/validation_release:$study/validation_release:ro"')
(dest/'runner.sh').write_text(s)
(dest/'verify_source.py').write_text('''from pathlib import Path
import sys
from common import verify,read,sha
s=Path(sys.argv[1]);verify(s)
for r in read(s/'validation_source_manifest.tsv'):
 p=s/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
print('BASE_AND_VALIDATION_SOURCE_PASS')
''')
ast.parse(g);ast.parse((dest/'verify_source.py').read_text())
subprocess.run(['bash','-n',str(dest/'runner.sh')],check=True)
subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(dest/'sif_checks.R'))+'))'],check=True,stdout=subprocess.DEVNULL)
table(study/'reports/additional_validation_tests.tsv',[dict(check=n,passed=True) for n in ('Python_syntax','R_parse','shell_parse','zero_additional_forward_calls','same_27_observer_call_scope','no_missing_probe_flat_classification')])
table(study/'validation_source_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.iterdir()) if p.is_file()])
(study/'reports/additional_validation.md').write_text('''# Pre-perturbation saved-reference validation

The new baseline is compared with original seed144 saved selections, full numeric
states and score, without any additional propagation. The SIF observer task also
records a source-based invariant-region check for all already registered vectors:
nonnegative steps, DT*lam_max<1 and all root/reseed totals below K imply K is
invariant because1-[y+a*y*(1-y)]=(1-y)*(1-a*y)>=0. G has nonnegative off-diagonals,
G_NN>=-lambda and column sums<=lam_max. Direct death only reduces the total.
This provides a sharper sufficient check when the original exponential gap bound
is inconclusive. Raw conservative branch receipts remain unchanged. The exact
polynomial identity is checked in G01. This is not a positive structural or
literal floating-point certificate. Missing feasible probes remain unresolved,
never classified numerically flat. No new parameter point, forward call,
contrast, precision, observation, fit or solver is introduced. The supplement
was frozen before perturbation outcomes; it is necessary validation, not a retry
of a failed allocation or an outcome-based change to a scientific threshold.
''')
print('SAVED_REFERENCE_AND_BRANCH_VALIDATION_FROZEN')
