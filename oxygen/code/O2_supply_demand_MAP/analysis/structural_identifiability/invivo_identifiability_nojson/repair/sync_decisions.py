#!/usr/bin/env python3
"""Mirror explicit metadata decisions and their local proofs; never control state."""
import argparse
import shlex
import subprocess
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    a = p.parse_args()
    study = a.study.resolve()
    if study.name != '20261111_invivo_baseline_repair_045318Z':
        raise ValueError('Unexpected study identity')
    names = [
        'baseline_gate*.tsv', 'baseline_policy_gates.tsv', 'baseline_chain_summary.tsv',
        'baseline_*.rds', 'extended_*_gate*.tsv', 'engineering_feasible_points.tsv',
        'cross_policy_feasible_floors.tsv', 'root_*_selected.rds', 'root_*_selected.tsv',
        'reserve_activation.tsv', 'reserve_decision.tsv',
        'density_*_selected.rds', 'profile_gate*.tsv', 'profile_point_review.tsv',
        'profile_counterexamples.rds', 'profile_reserve_activation.tsv',
        'profile_reserve_decision.tsv', 'profile_*_verified.rds',
    ]
    paths = set()
    for pattern in names:
        paths.update(study.joinpath('tables').glob(pattern))
    for directory in ('record_audits', 'gate_history'):
        paths.update(study.joinpath('tables', directory).glob('*.tsv'))
    files = sorted(str(path.relative_to(study)) for path in paths if path.is_file())
    if not files:
        raise ValueError('No metadata decisions found')
    listing = study / 'tables/decision_sync_files.txt'
    listing.write_text('\n'.join(files) + '\n')
    ssh = ['ssh', '-o', 'BatchMode=yes', '-o', 'IdentitiesOnly=yes', '-i',
           '/Users/4482173/.ssh/id_ed25519', '-o', 'ControlMaster=auto', '-o',
           'ControlPath=/private/tmp/ivi_repair_ssh_%C', '-o', 'ControlPersist=600']
    remote = '/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/' + study.name
    subprocess.run(['rsync', '-az', '--files-from=' + str(listing), '-e',
                    ' '.join(shlex.quote(x) for x in ssh), str(study) + '/',
                    '4482173@red.moffitt.org:' + remote + '/'], check=True)
    print('Mirrored', len(files), 'explicit metadata decision/proof files; no ledger, source or release copy.')


if __name__ == '__main__':
    main()
