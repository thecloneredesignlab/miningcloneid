#!/usr/bin/env python3
"""Current Figure7/8 entry plumbing and static-only checks (stdlib only)."""
import argparse
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def table(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def current_path(path, workspace, required=True):
    p = Path(path).expanduser()
    if not p.is_absolute():
        p = workspace / p
    for value in (str(p), str(p.resolve())):
        if re.search(r'(^|/)(arxives|archives|frozen_sources)(/|$)|/figure[78]_[^/]*_v[0-9]+(/|$)', value):
            raise ValueError('Old-version/archive roots are not runtime inputs: ' + value)
    p = p.resolve()
    if required and not p.exists():
        raise ValueError('Missing current input: ' + str(p))
    return p


def model_sources(workspace):
    root = current_path('Code/vendor/figure78_model', workspace)
    override = os.environ.get('FIGURE_MODEL_CODE_ROOT')
    if override and current_path(override, workspace) != root:
        raise ValueError('External/shared Figure7/8 model override is forbidden')
    rows = table(workspace/'Code/config/manifests/figure78_model_sources.tsv')
    if len(rows) != 12 or len({r['relative_path'] for r in rows}) != 12:
        raise ValueError('Invalid model lock')
    for r in rows:
        p = (root/r['relative_path']).resolve()
        if root not in p.parents or sha(p) != r['sha256']:
            raise ValueError('Private model lock mismatch: ' + str(p))
    return root, rows


def r_command(workspace):
    explicit = os.environ.get('FIGURE78_RSCRIPT')
    if explicit:
        return shlex.split(explicit)
    return ['bash', str(workspace/'Code/hpc/figure78_Rscript.sh')]


def check_prepared(workspace, figure, data_root, index=None):
    model_sources(workspace)
    root = current_path(data_root, workspace)
    if figure == 7 and root == (workspace/'data/Figures').resolve():
        raise ValueError('Published Figure7 is read-only for generation; prepare an isolated campaign root')
    run = root/('Figure'+str(figure))
    required = [run/'prepared/analysis_bundle.rds', run/'prepared/task_manifest.tsv',
                run/'provenance/preparation.tsv']
    if not all(p.is_file() for p in required):
        raise ValueError('A current prepared campaign is required; retained plotting data are not prepared tasks')
    receipt = table(required[2])
    tasks = table(required[1])
    if len(receipt) != 1 or receipt[0]['profile'] != 'full' or receipt[0]['status'] != 'prepared':
        raise ValueError('Only an explicitly prepared full campaign can be used on HPC')
    if not tasks or [int(r['task_id']) for r in tasks] != list(range(1, len(tasks)+1)):
        raise ValueError('Current task IDs must be complete and consecutive')
    if int(receipt[0]['n_tasks']) != len(tasks):
        raise ValueError('Preparation and task manifest counts disagree')
    if figure == 7:
        if receipt[0].get('engine') != 'Figure7_native_fit_last5_current':
            raise ValueError('Figure7 HPC requires the native-fit last-five current campaign; legacy caches are forbidden')
        for row in tasks:
            expected = 'native_fit_last5' if row['model_context'] == 'in vitro' else 'unchanged_continuous'
            if row.get('engine') != expected:
                raise ValueError('Figure7 context engine mismatch')
            if row['model_context'] == 'in vitro' and int(row['maximum_ploidy']) not in (137,139,141):
                raise ValueError('Figure7 full grid certificate requires 137/139/141N')
    if index is not None and not 1 <= int(index) <= len(tasks):
        raise ValueError('Task index outside current manifest')
    return len(tasks)


def static_check(workspace, parse_r=True, source_only=False):
    checks = []
    def check(name, condition, detail=''):
        checks.append({'check': name, 'passed': bool(condition), 'detail': detail})
        if not condition:
            raise ValueError(name+': '+str(detail))
    root, lock = model_sources(workspace)
    check('private_model_12_locked_files', len(lock) == 12)
    cfg = json.loads((workspace/'Code/config/figure78_inputs.json').read_text())
    check('model_default_is_private', cfg['model_root'] == 'Code/vendor/figure78_model')
    fit = current_path(os.environ.get('FIGURE_JOINT_RESULT_ROOT', cfg['fit_root']), workspace)
    inputs = [fit/f['pair_id']/('seed'+str(seed))/name for f in cfg['families'].values()
              for seed in range(1, 501) for name in ('joint_best_params_long.tsv','fit_config.rds','fit_result.rds')]
    check('all_3000_original_fit_inputs_exist', all(p.is_file() for p in inputs), len(inputs))
    manifest = [] if source_only else table(workspace/'Code/config/manifests/figure78_current_plot_inputs.tsv')
    current_figures = set()
    # Newly generated data authenticate against their own completed current
    # bundle, just as the renderer does. The old consolidation map is only for
    # figures that still retain that older provenance layout (e.g. Figure8).
    for number in (() if source_only else (7, 8)):
        tag = 'Figure'+str(number)
        run = workspace/'data/Figures'/tag
        bundle = run/'prepared/analysis_bundle.rds'
        if not bundle.is_file():
            continue
        done = table(run/'provenance/completion.tsv')
        check('current_full_completion_'+tag, len(done)==1 and
              done[0]['status']=='validated_complete' and done[0]['profile']=='full' and
              done[0]['production_ensemble']=='TRUE')
        record = done[0]
        for rel, field in [('prepared/analysis_bundle.rds','bundle_sha256'),
                           ('provenance/numerical_validation.tsv','validation_sha256'),
                           ('provenance/plot_inputs.tsv','plot_manifest_sha256')]:
            check('current_receipt_'+tag+'_'+field, sha(run/rel)==record[field])
        validation = table(run/'provenance/numerical_validation.tsv')
        check('current_numerical_gates_'+tag, bool(validation) and all(r['passed']=='TRUE' for r in validation))
        rows = table(run/'provenance/plot_inputs.tsv')
        for r in rows:
            r['figure']=tag
        historical = [r for r in manifest if tag=='Figure7' and r['figure']==tag and
                      r.get('group')=='lineage' and
                      r['relative_path'] not in {x['relative_path'] for x in rows} and
                      (workspace/'data/Figures'/r['relative_path']).is_file()]
        if tag=='Figure7':
            wanted={'history_time_summary.tsv','history_source_daughter_flux_binned.tsv',
                    'history_low_population_budget.tsv'}
            check('preserved_original_historical_summary_bindings_'+tag,
                  wanted <= {r['logical_name'] for r in historical+rows})
        manifest = [r for r in manifest if r['figure']!=tag]+historical+rows
        for r in table(run/'provenance/current_code_manifest.tsv'):
            slot = re.sub(r'^.*/(Code/)',r'\1',r['path'])
            check('current_code_portable_slot_'+tag+'_'+Path(slot).name, slot.startswith('Code/'))
            p = current_path(slot, workspace)
            check('current_code_hash_'+tag+'_'+Path(slot).name,sha(p)==r['sha256'])
        current_figures.add(tag)
    for row in manifest:
        p = workspace/'data/Figures'/row['relative_path']
        check('retained_input_'+row['relative_path'], p.is_file() and sha(p) == row['sha256'])
        check('retained_input_not_archive_'+row['relative_path'],
              not any(x in p.resolve().parts for x in ('arxives','archives')))
    if not source_only:
        check('current_and_retained_input_inventory', bool(manifest) and
              len({r['relative_path'] for r in manifest})==len(manifest))
    for rel in ('data/Figure3_PCAWG_included_samples.csv', 'data/Figure3_DepMap_included_lines.csv', 'originals/organ_lut.csv'):
        check('empirical_'+rel, (workspace/'data/share_bundle_figures_3_5'/rel).is_file())
    bindings = [] if source_only else table(workspace/'Code/config/manifests/figure78_bundle_bindings.tsv')
    by_name = {Path(r['relative_path']).name:r for r in lock}
    for row in bindings:
        if row['relative_path'].split('/')[0] in current_figures:
            continue
        p = workspace/'data/Figures'/row['relative_path']
        check('retained_binding_'+row['binding']+'_'+row['order'], sha(p) == row['sha256'])
        if row['binding'] in ('engine','initializer'):
            r = by_name[p.name]
            check('private_model_matches_original_'+p.name, r['sha256'] == row['sha256'])
    entries = table(workspace/'Code/config/manifests/figure78_entrypoints.tsv')
    check('entry_names_portable_across_case_sensitive_filesystems',
          len({r['path'].lower() for r in entries}) == len(entries))
    global_rows=[r for r in table(workspace/'Code/config/manifests/figure_entrypoints.tsv')
                 if re.match(r'^(Figure[78](?:_ver2)?|Supp_Figure[78](?:_|$))',r['figure'])]
    check('global_registry_all_27_publication_products',len(global_rows)==27)
    for row in global_rows:
        figure=7 if re.match(r'^(Figure7|Supp_Figure7)',row['figure']) else 8
        check('global_registry_current_'+row['figure'],
              row['analysis_entry']=='Code/Figures/data_Figure%d.R'%figure and
              row['data_directory']=='data/Figures/Figure%d'%figure and
              (workspace/row['drawing_entry']).is_file())
    for row in entries:
        p = workspace/row['path']
        check('entry_exists_'+row['path'], p.is_file())
        if row['status'] == 'compatibility_guard':
            check('legacy_entry_routed_'+row['path'],
                  any(s in p.read_text() for s in ('figure78_compat.R','legacy_main','figure78_legacy_hpc.sh','Legacy Slurm spool')))
    current = set(row['path'] for row in entries)
    for p in (workspace/'Code/Figures').glob('*'):
        if p.suffix in ('.R','.py') and re.search(r'figure[78]|fig8',p.name,re.I):
            check('all_figure_entries_registered_'+p.name, str(p.relative_to(workspace)) in current)
    drawings = [r for r in entries if r['role'] == 'drawing']
    check('32_drawing_entries', len(drawings) == 32)
    for r in drawings:
        check('unified_draw_'+r['path'], 'figure78_render_only.R' in (workspace/r['path']).read_text())
    for p in (workspace/'Code/hpc').iterdir():
        if p.suffix in ('.sh','.sub') and re.search(r'figure[78]|fig8',p.name,re.I):
            check('all_hpc_entries_registered_'+p.name, str(p.relative_to(workspace)) in current)
    retired = json.loads((workspace/'Code/config/manifests/figure78_retired_paths.json').read_text())
    for rel in retired:
        check('retired_implementation_not_active_'+rel, not (workspace/rel).exists())
    analysis = workspace/'Code/Figures/util/analysis'
    modules = [p for p in analysis.iterdir() if p.suffix in ('.R','.cpp') and re.search(r'figure[78]',p.name,re.I)]
    for p in modules:
        s = p.read_text()
        check('no_version_directory_dependency_'+p.name,
              not re.search(r'figure[78]_[A-Za-z0-9_]+_20[0-9]{6}_v[0-9]+',s))
    python_files = [workspace/r['path'] for r in entries if r['path'].endswith('.py')]
    for p in python_files:
        ast.parse(p.read_text(),filename=str(p))
    check('python_ast_only', True, len(python_files))
    shells = [workspace/r['path'] for r in entries if r['path'].endswith(('.sh','.sub','.env'))]
    for p in shells:
        result = subprocess.run(['bash','-n',str(p)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        check('shell_syntax_'+str(p.relative_to(workspace)), result.returncode == 0, result.stderr)
    rfiles = sorted(set([workspace/r['path'] for r in entries if r['path'].endswith('.R')]
                       + [p for p in modules if p.suffix == '.R']
                       + [root/r['relative_path'] for r in lock if r['relative_path'].endswith('.R')]))
    if parse_r:
        expression = 'for(p in commandArgs(TRUE)) invisible(parse(file=p,keep.source=FALSE));cat("PARSE_ONLY_OK\\n")'
        result = subprocess.run(r_command(workspace)+['--vanilla','-e',expression,*map(str,rfiles)],
                                stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
        check('r_parse_only_no_source_or_model_execution', result.returncode == 0,
              {'files':len(rfiles),'stdout':result.stdout,'stderr':result.stderr})
    else:
        check('r_parse_not_requested', True, 'not tested')
    worker = (workspace/'Code/hpc/run_Figure78_task.sh').read_text()
    check('worker_auto_loads_private_paths', 'figure78_hpc_paths.env' in worker and 'check-prepared' in worker)
    submit = (workspace/'Code/hpc/submit_Figure78.sh').read_text()
    check('submit_explicit_opt_in', '--submit=yes' in submit)
    check('xxlarge_without_array_concurrency_cap', '--qos=xxlarge' in submit and '--array="1-${count}"' in submit)
    check('renderer_uses_private_initializer_not_snapshot',
          'sys.source(private,envir=env)' in (analysis/'figure78_provenance.R').read_text())
    negative_roots = ('data/Figures/arxives/not_read','data/Figures/archives/not_read',
                      'data/Figures/Figure7/figure7_instant_history_20261014_v8')
    for value in negative_roots:
        rejected = False
        try:
            current_path(value,workspace,required=False)
        except ValueError:
            rejected = True
        check('reject_old_root_without_read_'+value,rejected)
    for figure in (7,8):
        run = workspace/'data/Figures'/('Figure'+str(figure))
        if figure == 7 or not (run/'prepared/analysis_bundle.rds').exists():
            rejected = False
            try:
                check_prepared(workspace,figure,workspace/'data/Figures')
            except ValueError:
                rejected = True
            check('retained_Figure'+str(figure)+'_cannot_be_submitted_as_prepared',rejected)
    source_paths = sorted(set([workspace/r['path'] for r in entries]+modules+
                            [workspace/'Code/config/figure78_inputs.json',
                             workspace/'Code/config/manifests/figure78_entrypoints.tsv',
                             workspace/'Code/config/manifests/figure78_model_sources.tsv']))
    return {'status':'source_static_passed' if source_only else 'static_passed',
            'retained_dataset_validation':not source_only,'checks':checks,'registered_entries':len(entries),
            'locked_model_files':len(lock),'retained_inputs':len(manifest),
            'r_syntax_files':len(rfiles),'python_syntax_files':len(python_files),'shell_syntax_files':len(shells),
            'source_sha256':{str(p.relative_to(workspace)):sha(p) for p in source_paths},
            'global_figure78_entry_rows':global_rows,
            'numerical_or_render_entry_execution':False,'jobs':0,'simulations':0,'redraws':0}


def legacy_main(name):
    # Legacy validators now expose only the explicit unified interfaces.
    workspace = Path(__file__).resolve().parents[2]
    args = sys.argv[1:]
    if '--mode=static' in args:
        args.remove('--mode=static')
        return main(['static','--workspace='+str(workspace),*args])
    if '--mode=retained' in args:
        args.remove('--mode=retained')
        return subprocess.call(r_command(workspace)+[str(workspace/'Code/Figures/validate_Figure78.R'),*args])
    raise SystemExit(name+': legacy results/indices are not accepted. Use '+
                     'figure78_entry_tools.py static or validate_Figure78.R; aliases require --mode=static|retained.')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=('static','check-prepared'))
    ap.add_argument('--workspace',type=Path,default=Path(__file__).resolve().parents[2])
    ap.add_argument('--figure',type=int,choices=(7,8))
    ap.add_argument('--data-root',type=Path)
    ap.add_argument('--task-index',type=int)
    ap.add_argument('--report',type=Path)
    ap.add_argument('--parse-r',choices=('yes','no'),default='yes')
    ap.add_argument('--source-only',choices=('yes','no'),default='no',
                    help='Syntax/entry/model/raw-input checks only; explicitly does not certify retained numerical data')
    a = ap.parse_args(argv)
    workspace = current_path(a.workspace.resolve(),Path.cwd())
    if a.action == 'check-prepared':
        if a.figure is None or a.data_root is None:
            ap.error('check-prepared requires --figure and --data-root')
        print(check_prepared(workspace,a.figure,a.data_root,a.task_index))
        return 0
    result = static_check(workspace,a.parse_r == 'yes',a.source_only == 'yes')
    if a.report:
        report = a.report.resolve()
        if workspace/'audit' not in report.parents:
            raise ValueError('Static reports must be written under workspace/audit')
        report.parent.mkdir(parents=True,exist_ok=True)
        report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('checks','source_sha256','global_figure78_entry_rows')},sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
