"""Prepare an isolated Stage A-Alternative audit; no scientific application calls."""
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

def table(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def read(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def git(repo, *args):
    return subprocess.check_output(['git', '--git-dir='+str(repo), *args])

def main(project, review, study):
    if study.exists() and (study/'inputs/source_manifest.tsv').exists():
        raise SystemExit('Refusing to overwrite an initialized evidence namespace')
    for name in ['inputs', 'reports', 'sources', 'controls', 'runtime', 'accounting', 'releases/r01']:
        (study/name).mkdir(parents=True, exist_ok=True)
    old = project/'oxygen/results/structural_identifiability/20261116_culture_step04_stageA_tools_01'
    verified = []
    for row in read(old/'reports/frozen_source_verification.tsv'):
        actual = digest(project/row['path'])
        assert actual == row['expected_sha256'], row['path']
        verified.append(dict(path=row['path'], sha256=actual, status='matched'))
    table(study/'reports/frozen_source_verification.tsv', verified)
    for source in [project/'oxygen/results/analysis_documentation/in_vitro/structural_certificates/culture_structural_step04_public_nonrational_alternative_plan_20261116.md',
                   project/'oxygen/results/analysis_documentation/frozen_inputs/in_vitro/design_and_contracts/meeting_aligned_synthetic_identifiability_master_plan_20261112.md',
                   old/'inputs/parameter_inventory.tsv', old/'reports/parameter_and_nuisance_audit.tsv',
                   old/'reports/per_case_backend_status.tsv']:
        shutil.copy2(source, study/'inputs'/source.name)
    roles = read(study/'inputs/parameter_inventory.tsv')
    assert len(roles) == 20 and sum(r['scientific_role']=='mechanism' for r in roles)==14
    table(study/'inputs/authorization.tsv', [dict(key=k, value=v) for k,v in {
        'scope':'Step 4 Stage A-Alternative only', 'confirmation_date':'2026-11-16 America/New_York',
        'independent_allocated_core_hour_cap':'24', 'qos':'xxlarge', 'concurrency_cap':'none',
        'author_contact':'not_authorized', 'Java_JavaScript_JSON':'forbidden',
        'Stage_B_Step_5':'not_authorized', 'synthetic_recovery_profiles_tumor':'not_authorized',
        'old_controllers_tasks':'preserve_unchanged'}.items()])
    packages = {
        'FISPO': ('strike-goldd.git','c3caa2eb3af0c7508470e3bbf3e034e3334dfa65','https://github.com/afvillaverde/strike-goldd','GPL-3.0'),
        'GenSSI': ('GenSSI.git','9b1bba281f782a59a851e726fc83197836be4115','https://github.com/genssi-developer/GenSSI','BSD-2-Clause text in README'),
        'StrikePy': ('StrikePy.git','26386200d116077158d37f5ade7bea1c219d8692','https://github.com/afvillaverde/StrikePy','GPL-3.0')}
    manifest, access, features = [], [], []
    terms = ('json', 'java', 'jvm', 'webread', 'webwrite', 'system(', 'subprocess', 'figure(', 'eval(')
    for tool, (dirname, commit, url, license_name) in packages.items():
        repo = review/dirname
        assert git(repo,'rev-parse','HEAD').decode().strip()==commit
        paths=git(repo,'ls-tree','-r','--name-only','HEAD').decode().splitlines()
        for name in paths:
            p = Path(name)
            select = name in ('LICENSE','README.md')
            if tool=='FISPO':
                select |= name in ('STRIKE-GOLDD/STRIKE_GOLDD.m','STRIKE-GOLDD/options.m')
                select |= name.startswith('STRIKE-GOLDD/functions/') and p.suffix=='.m' and 'aux_APP/' not in name
            elif tool=='GenSSI':
                select |= p.suffix=='.m' and (name.startswith('Auxiliary/') or len(p.parts)==1)
            else:
                select |= p.suffix=='.py' and ('/functions/' in name or name in ('StrikePy/strike_goldd.py','StrikePy/options.py'))
            if not select:
                continue
            data=git(repo,'show',commit+':'+name)
            dest=study/'sources'/tool/name
            dest.parent.mkdir(parents=True,exist_ok=True)
            dest.write_bytes(data)
            manifest.append(dict(tool=tool,commit=commit,upstream_path=name,path=str(dest.relative_to(study)),bytes=len(data),sha256=digest(dest)))
            if p.suffix in ('.py','.m'):
                for number,line in enumerate(data.decode('utf-8', errors='replace').splitlines(),1):
                    found=[term for term in terms if term in line.lower()]
                    if found:
                        features.append(dict(tool=tool,path=name,line=number,features=','.join(found),text=line.strip(),interpretation='static flag; reachability requires audit'))
        assert (study/'sources'/tool/'README.md').is_file()
        license_file=study/'sources'/tool/('README.md' if tool=='GenSSI' else 'LICENSE')
        assert 'Redistribution and use' in license_file.read_text() if tool=='GenSSI' else 'GNU GENERAL PUBLIC LICENSE' in license_file.read_text()
        access.append(dict(tool=tool,origin=url,commit=commit,license=license_name,source_obtained=True,license_obtained=True,upstream_execution=False))
    table(study/'inputs/source_manifest.tsv',manifest)
    table(study/'reports/source_access.tsv',access)
    table(study/'reports/static_source_flags.tsv',features)
    table(study/'accounting/budget_ledger.tsv',[dict(phase=k,cap_core_hours=v,actual_core_seconds=0,reserved_core_seconds=0) for k,v in [('runtime',4),('method_controls',8),('nuisance_controls',6),('crosscheck',2),('report',2),('one_technical_retry',2)]])
    table(study/'accounting/jobs.tsv',[],['task','job_id','qos','state','exit_code','allocated_cpus','elapsed_seconds','actual_core_seconds','reserved_core_seconds'])
    tracked=subprocess.check_output(['git','diff','--name-only'],cwd=project,text=True).splitlines()
    table(study/'inputs/protected_worktree_hashes.tsv',[dict(path=p,sha256=digest(project/p)) for p in tracked])
    shutil.copy2(__file__,study/'releases/r01/prepare.py')
    print('Prepared',len(manifest),'licensed source files;',len(verified),'frozen source hashes;',len(roles),'original roles. No upstream application executed.')

if __name__=='__main__':
    main(*(Path(arg).resolve() for arg in sys.argv[1:]))
