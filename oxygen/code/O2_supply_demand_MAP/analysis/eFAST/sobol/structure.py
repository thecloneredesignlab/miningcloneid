"""Auditable parameter/rate/operator map for the separate fixed-O2 in-vivo model."""
import argparse
from pathlib import Path
import core as c

DIRECT = {
    'lam_max':['growth'], 'alpha_o2':['growth'], 'gamma_growth':['growth'],
    'mu_hp':['death','missegregation'], 'gamma_mu':['death','missegregation'],
    'O2_crit':['oxygen_stress','growth','death','missegregation'],
    'n_O':['oxygen_stress','growth','death','missegregation'],
    'p_mis_base':['missegregation'], 'p_misseg':['missegregation'], 'k_o_mis':['missegregation'],
    'buffer_smax':['buffering'], 'buffer_beta':['buffering'], 'buffer_n_exp':['buffering'],
    'p_wgd':['WGD'],
}
DETAIL = {
    'lam_max':('lambda_eff_soft_o2_only_cpp','lam_max / (1 + alpha_o2*h*(N/44)^gamma_growth)',270,'division-rate scale'),
    'alpha_o2':('lambda_eff_soft_o2_only_cpp','lam_max / (1 + alpha_o2*h*(N/44)^gamma_growth)',270,'hypoxia/ploidy growth penalty strength'),
    'gamma_growth':('lambda_eff_soft_o2_only_cpp','lam_max / (1 + alpha_o2*h*(N/44)^gamma_growth)',270,'ploidy dependence of growth penalty'),
    'mu_hp':('mu_eff_soft_cpp; resolve_pmis_for_death','mu_hp*h*(N/44)^gamma_mu; p0+a*mu/(mu+k)',378,'death-rate scale also enters missegregation saturation'),
    'gamma_mu':('mu_eff_soft_cpp; resolve_pmis_for_death','mu_hp*h*(N/44)^gamma_mu; p0+a*mu/(mu+k)',378,'ploidy dependence of death also changes missegregation'),
    'O2_crit':('hypoxia_weight_cpp','h=O2_crit^n_O / (O2_crit^n_O + O2^n_O)',184,'shared oxygen scale for growth, death and death-linked missegregation'),
    'n_O':('hypoxia_weight_cpp','h=O2_crit^n_O / (O2_crit^n_O + O2^n_O)',184,'shared oxygen response steepness'),
    'p_mis_base':('resolve_pmis_for_death','p_mis=clip01(p_mis_base+p_misseg*mu/(mu+k_o_mis))',951,'baseline per-copy missegregation probability'),
    'p_misseg':('resolve_pmis_for_death','p_mis=clip01(p_mis_base+p_misseg*mu/(mu+k_o_mis))',951,'death-linked missegregation amplification'),
    'k_o_mis':('resolve_pmis_for_death','p_mis=clip01(p_mis_base+p_misseg*mu/(mu+k_o_mis))',951,'half-saturation scale on effective death rate'),
    'buffer_smax':('buffering_survival_modifier','s_N=clip01(buffer_smax*exp(-buffer_beta*(44/N)^buffer_n_exp)); survival=s_N^m for m>0',533,'per-event survival ceiling'),
    'buffer_beta':('buffering_survival_modifier','s_N=clip01(buffer_smax*exp(-buffer_beta*(44/N)^buffer_n_exp)); survival=s_N^m for m>0',533,'ploidy-dependent buffering penalty strength'),
    'buffer_n_exp':('buffering_survival_modifier','s_N=clip01(buffer_smax*exp(-buffer_beta*(44/N)^buffer_n_exp)); survival=s_N^m for m>0',533,'ploidy dependence of buffering penalty'),
    'p_wgd':('constant_p_wgd_cpp; cpp_o2simps_build_G_for_o2_triplet','scale_mis=lambda_N*(1-p_wgd); scale_wgd=lambda_N*p_wgd',1059,'division channel split between missegregation kernel and WGD'),
}


def pair_link(pi,pj):
    a,b=set(DIRECT[pi]),set(DIRECT[pj]); shared=a&b
    if shared: return 'shared_rate_or_upstream_pathway', ';'.join(sorted(shared))
    if ('buffering' in a and 'missegregation' in b) or ('buffering' in b and 'missegregation' in a):
        return 'missegregation_kernel_survival_weighting','binomial event weights multiplied by buffering survival'
    if ('growth' in a and ('buffering' in b or 'missegregation' in b or 'WGD' in b)) or ('growth' in b and ('buffering' in a or 'missegregation' in a or 'WGD' in a)):
        return 'division_generator_product','division rate scales offspring transition weights'
    if ('WGD' in a and ('buffering' in b or 'missegregation' in b)) or ('WGD' in b and ('buffering' in a or 'missegregation' in a)):
        return 'division_channel_weighting','1-p_wgd scales the missegregation/buffering channel'
    return 'coupled_eigenproblem','parameters enter the same live-state operator; output interaction must be measured'


def write_structure(out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    cpp=c.HERE.parents[2]/'model/model_O2_supply_demand_MAP.cpp'
    cpp_text=cpp.read_text()
    rows=[]
    for p in c.ACTIVE:
        function,formula,line,meaning=DETAIL[p]
        if any(name.strip() not in cpp_text for name in function.split(';')):
            raise ValueError('Mapped model function absent: '+function)
        rows.append(dict(parameter=p,primary_group=c.efast.GROUP[p],pathways=';'.join(DIRECT[p]),
            function=function,formula=formula,interpretation=meaning,present_fixed_O2='TRUE',
            code_file='oxygen/code/O2_supply_demand_MAP/model/model_O2_supply_demand_MAP.cpp',code_line=line,
            model_cpp_sha256=c.sha(cpp)))
    for p in c.efast.STRUCTURAL:
        rows.append(dict(parameter=p,primary_group='oxygen_dynamics_or_dead_clearance',pathways='',function='',formula='',
            interpretation='absent from the externally fixed-oxygen live-state operator; outside the 14-factor design',
            present_fixed_O2='FALSE',code_file='',code_line='',model_cpp_sha256=c.sha(cpp)))
    c.table(out/'parameter_equation_map.tsv',rows)
    pair_rows=[]
    for i,j in c.PAIRS:
        pi,pj=c.ACTIVE[i],c.ACTIVE[j]; category,link=pair_link(pi,pj)
        pair_rows.append(dict(parameter_i=pi,parameter_j=pj,process_i=c.efast.GROUP[pi],process_j=c.efast.GROUP[pj],
            structural_link_category=category,structural_link=link,
            evidence_status='mechanistic candidate; neither large S2 nor nonidentifiability is established'))
    c.table(out/'parameter_pair_mechanism_map.tsv',pair_rows)
    (out/'model_equations.md').write_text('''# Separate in-vivo fixed-oxygen model: parameter and equation map

This is a mechanistic equation audit, not a completed structural-identifiability proof.
The canonical model is `model_O2_supply_demand_MAP.cpp`, with N=22,...,154,
N_UNIT=22 and diploid chromosome-copy reference 44. The in-vivo configuration
enables O2_growth and uses ploidy_O2_death="ploidy_related".

## Rates and shared pathways

For oxygen O (percent), h(O)=C^n/(C^n+O^n), with C=O2_crit and n=n_O.
Growth: lambda_N=lam_max/[1+alpha_o2*h(O)*(N/44)^gamma_growth].
Death: mu_N=mu_hp*h(O)*(N/44)^gamma_mu.
Missegregation: p_N=clip01[p_mis_base+p_misseg*mu_N/(mu_N+k_o_mis)].
Thus the death module explicitly feeds the missegregation module. At O=0,
h=1 for all allowed C,n, so both oxygen-shape parameters are inactive at that
specific oxygen condition; their interaction estimates should be near zero
within numerical/sampling error.

For m missegregated copies and parent N, s_N=clip01[buffer_smax*
exp(-buffer_beta*(44/N)^buffer_n_exp)]. Survival is s_N^m for m>0 and 1 for m=0.
Binomial weights Binomial(m;N,p_N) multiply this survival and populate the
gain/loss daughter shifts N+m and N-m symmetrically. Kernel tail truncation
and the configured state-boundary rule remain exactly those of the fitted model.
The kernel weights represent expected viable daughter numbers, not a
renormalized conditional probability distribution.

WGD weight is lambda_N*p_wgd; the ordinary division-kernel weight is
lambda_N*(1-p_wgd). The code deposits the WGD contribution at 2N and removes
the mother at rate lambda_N. In schematic matrix form,
G[q,N]=lambda_N*((1-p_wgd)*K[q,N]+p_wgd*B[q,N])-lambda_N*delta[q,N],
with K and B including the original boundary handling.
M=G-diag(mu_N); dx/dt=Mx at fixed oxygen.

## Outputs

```mermaid
flowchart LR
    O[Oxygen, O2_crit, n_O] --> H[Hypoxia weight]
    H --> L[Division rate]
    H --> D[Death rate]
    D --> P[Missegregation probability]
    P --> K[Viable daughter kernel]
    B[Buffering parameters] --> K
    L --> G[Division generator G]
    K --> G
    W[p_wgd] --> G
    G --> M[Live-state operator M]
    D --> M
    M --> E[Leading eigenvalue: growth]
    M --> V[Leading eigenvector: mean ploidy]
```

Growth r=max Re(eigenvalues(M)), in day^-1.
Dominant mean ploidy=sum_N (N/22)*v_N, where Mv=rv and sum(v)=1.
The Sobol evaluator uses this unchanged canonical matrix, checks the eigenpair
residual and nonnegativity, and records 50/100-digit recovery if required.

Parameters interact through rate products, saturation, offspring survival,
channel weighting and the eigenvalue/eigenvector map. Co-occurrence in an
equation is a mechanistic link; it does not prove a large second-order Sobol
component. The 91-pair table records candidate links before numerical analysis.

## Scope and identifiability

o2_S0, kappa_O, eta_o2 and k_clear do not enter this fixed-O2 live-state operator.
The parameter-equation map is not a claim that the 14 active parameters are
individually identifiable. COMBOS/StructuralIdentifiability analyses would
require a specified observation model and a separate audit of the model's
non-rational parameter functions.
''')
    c.save_json(out/'source_manifest.json',dict(model_cpp_sha256=c.sha(cpp),
        fixed_O2_script_sha256=c.sha(c.HERE.parents[2]/'simulation/o2/fixed_o2/run_fixed_o2_simulation.R'),
        n_active_parameters=14,n_absent_parameters=4,n_pairs=91,
        structural_identifiability_tested=False))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--out-dir',required=True)
    write_structure(p.parse_args().out_dir)
