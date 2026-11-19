"""Assemble and verify a new release without changing closed archives."""
from pathlib import Path
import shutil
import sys

from common import sha, table, verify

ROOT_NAME = '20261118_culture_D19_independent_01'
UPSTREAM = '20261117_culture_D19_recovery_restart_01'

def main(parent, release, starts):
    parent = Path(parent).resolve()
    previous = parent / UPSTREAM
    root = parent / ROOT_NAME
    root.mkdir()
    for subdir in ('verification', 'controller', 'accounting', 'logs', 'reports',
                   'inputs/public', 'inputs/private', 'inputs/canaries'):
        (root / subdir).mkdir(parents=True)
    shutil.copytree(release, root / 'release', symlinks=True)
    shutil.copytree(previous / 'runtime', root / 'runtime', symlinks=True)
    shutil.copytree(previous / 'inputs/production', root / 'inputs/production', symlinks=True)
    for item in ('D19.rds', 'contract.rds', 'parameter_roles.tsv'):
        shutil.copy2(previous / 'inputs/public' / item, root / 'inputs/public')
    for name, digest in (
        ('D19.rds','0b1787b85af14c024c420c281c59c50ef92ac5a4f10ada926f1cc614ced1cbfe'),
        ('contract.rds','db3a2d7aed3356facf901e31d887548649920c4c8db9689706d2018cd408a618')):
        if sha(root / 'inputs/public' / name) != digest:
            raise RuntimeError('Unmatched frozen input: ' + name)
    for item in ('truth.rds', 'reference.rds'):
        shutil.copy2(previous / 'inputs/private' / item, root / 'inputs/private')
    for item in ('banks', 'rng'):
        shutil.copytree(Path(starts) / item, root / 'inputs' / item)
    shutil.copy2(Path(starts) / 'start_registry.tsv', root / 'inputs')
    (root / 'inputs/canaries/alias_truth.rds').symlink_to('../private/truth.rds')
    (root / 'inputs/canaries/alias_B.rds').symlink_to('../banks/D_B.rds')
    canary = root / 'verification/D_B_fixture'
    canary.mkdir()
    (canary / 'state.rds').write_text('Foreign origin read canary.\n')
    rows = []
    for subdir in ('release', 'runtime', 'inputs'):
        for path in sorted((root / subdir).rglob('*')):
            if path.is_file() and not path.is_symlink():
                rows.append(dict(path=str(path.relative_to(root)), bytes=path.stat().st_size, sha256=sha(path)))
    table(root / 'source_manifest.tsv', rows)
    verify(root)
    table(root / 'accounting/authority.tsv', [dict(
        scope='D19 independent blind validation A/B/D with conditional F',
        qos='xxlarge', core_seconds=432000, max_full_forwards=336400,
        mechanisms=14, nuisances=6, technical_retries=1, reserve_rounds=1)])
    print(root)

if __name__ == '__main__':
    main(*sys.argv[1:])
