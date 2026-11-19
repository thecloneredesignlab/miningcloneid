"""Re-read candidate intervals and verify repeated certificates and observer factorization."""
from pathlib import Path

from certify import arb, arb_mat, ctx, exact, nonzero, independent_minor_check, read, table
from algebra import decode, encode

OUT = Path(__file__).parent
TAGS = ('growth_only_18', 'identity_kernel_18', 'real_kernel_19')


def load(directory, tag):
    path = OUT / directory
    rows = read(path / f'{tag}_rows.tsv')
    cols = read(path / f'{tag}_cols.tsv')
    entries = read(path / f'{tag}_matrix.tsv')
    assert len(rows) == len(cols) and len(entries) == len(rows)**2
    matrix = arb_mat(len(rows), len(cols))
    for i, row in enumerate(entries):
        assert (int(row['row']), int(row['column'])) == divmod(i, len(cols))
        matrix[int(row['row']), int(row['column'])] = decode(row)
    check, norm = independent_minor_check(matrix)
    det = matrix.det()
    assert check and nonzero(det)
    assert det.overlaps(decode(read(path / f'{tag}_determinant.tsv')[0]))
    return matrix, rows, cols, det, norm


def main():
    ctx.prec = 640
    summaries = []
    anchors = {}
    for tag in TAGS:
        a = load('anchor_384', tag)
        b = load('anchor_512', tag)
        c = load('interior_384', tag)
        assert a[1:3] == b[1:3]
        count = a[0].nrows()
        assert all(a[0][i,j].overlaps(b[0][i,j]) for i in range(count) for j in range(count))
        assert a[3].overlaps(b[3])
        anchors[tag] = a
        summaries.append(dict(map=tag, rank=count, anchor_determinant=str(a[3]),
                              interior_determinant=str(c[3]), anchor_neumann_bound=str(a[4]),
                              entries_overlapping_across_anchor_precisions=count**2,
                              rereview_bits=ctx.prec))
    table(OUT / 'verified_rank.tsv', summaries)

    for filename, keys in [
        ('operator_branches.tsv', ('oxygen_pct', 'mask_sha256')),
        ('passage_branches.tsv', ('cohort', 'segment', 'selected_day', 'selection_rule',
                                 'parent_segment', 'parent_day', 'reset_rule')),
    ]:
        reference = None
        for run in ('anchor_384','anchor_512','interior_384'):
            signature = [tuple(r[k] for k in keys) for r in read(OUT / run / filename)]
            if reference is None:
                reference = signature
            assert reference == signature

    M, rows, cols, _, _ = anchors['real_kernel_19']
    assert cols[-1]['parameter'] == 'sigma_kary'
    observed = {(r['passage_id'], r['observed_kary_N'])
                for r in read(Path(__file__).resolve().parents[4] / 'data/InVitroData_LTEE/invitro_kary_cells.tsv')}
    possible = [i for i, row in enumerate(rows) if (row['passage'], row['chromosome']) in observed]
    chosen = max(possible, key=lambda i: abs(M[i,18]).mid())
    entry = M[chosen,18] * arb((1, -int(cols[-1]['rescale_power_of_2'])))
    assert nonzero(entry)
    # Growth rows have exactly zero width derivative. Appending this observed
    # bin as the last row makes the corresponding minor block triangular.
    factorized = anchors['growth_only_18'][3] * entry
    assert nonzero(factorized)
    table(OUT / 'observer_factorization.tsv', [dict(**rows[chosen], **encode(entry),
          raw_width_partial=str(entry), determinant_with_growth_column_scaling=str(factorized))])
    print('Verified 9 saved matrices, all repeated anchor entries, and the observed-bin width factorization.')
    print('Observed bin', rows[chosen]['passage'], rows[chosen]['chromosome'], 'width partial', entry)
    print('Minors and complete branch signatures match the declared checks.')


if __name__ == '__main__':
    main()
