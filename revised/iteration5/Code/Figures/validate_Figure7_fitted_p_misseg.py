#!/usr/bin/env python3
"""Independent table/figure QA; never imports a simulator or submits a job."""
import argparse
import hashlib
import json
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--campaign-root', required=True, type=Path)
    parser.add_argument('--audit-root', required=True, type=Path)
    parser.add_argument('--output-root', type=Path)
    parser.add_argument('--protected-before', type=Path)
    parser.add_argument('--workspace', type=Path)
    args = parser.parse_args()
    root = args.campaign_root
    args.audit_root.mkdir(parents=True, exist_ok=False)
    checks = []

    def check(name, passed, detail=None):
        checks.append(dict(name=name, passed=bool(passed), detail=detail))
        (args.audit_root / 'checks.json').write_text(json.dumps(checks, indent=2))
        if not passed:
            raise AssertionError(name)

    def table(relative):
        return pd.read_csv(root / relative, sep='\t')

    done = json.loads((root / 'provenance/completion.json').read_text())
    check('validated_complete', done['status'] == 'validated_complete')
    check('bundle_authenticated', sha(root / 'prepared/analysis_bundle.rds') == done['bundle_sha256'])
    check('manifest_authenticated', sha(root / 'provenance/validated_files.tsv') == done['manifest_sha256'])
    mf = table('provenance/validated_files.tsv')
    check('all_registered_files_authenticated', all(sha(root / r.relative_path) == r.sha256 for r in mf.itertuples()))
    expected = table('provenance/fitted_p_misseg.tsv').set_index(['pair_label', 'model_context']).p_misseg
    z = table('tables/instant_time_summary.tsv')
    check('fitted_mode_only', z.condition_mode.eq('fitted_p_misseg').all())
    check('native_p_per_environment', all(np.isclose(r.p_misseg, expected.loc[(r.pair_label, r.model_context)], atol=1e-14, rtol=0)
          for r in z[['pair_label', 'model_context', 'p_misseg']].drop_duplicates().itertuples()))
    key = ['pair_label', 'model_context', 'O2_pct', 'initial_ploidy', 'q', 'day']
    check('no_duplicate_daily_metrics', not z.duplicated(key + ['metric']).any())
    a = z[z.metric.eq('flux_percent')].set_index(key).sort_index()
    b = z[z.metric.eq('growth_percent')].set_index(key).sort_index()
    check('matched_source_keys', a.index.equals(b.index))
    check('matched_source_NA', a['median'].isna().equals(b['median'].isna()))
    err = (a['median'] + b['median'] - 100).abs().max()
    check('A_B_complement', err < 1e-10, float(err))
    check('reverse_IQR_complement', max((a.q25 + b.q75 - 100).abs().max(),
          (a.q75 + b.q25 - 100).abs().max()) < 1e-10)
    check('day_zero_NA', z.loc[z.day.eq(0), 'median'].isna().all())
    check('source_shares_in_range', z['median'].dropna().between(-1e-10, 100 + 1e-10).all())
    check('48_complete_conditions', len(z[key[:-1]].drop_duplicates()) == 48)
    counts = table('tables/daily_source_counts_summary.tsv')
    cp = counts.pivot(index=key, columns='counter', values='median')
    reconstructed = 100 * cp.HL_per_start_live_cell / (cp.HL_per_start_live_cell + cp.LL_per_start_live_cell)
    # Protocol medians of counts are not generally the median of ratios. The
    # deterministic continuous route is the appropriate exact independent cut.
    vivo = a.index.get_level_values('model_context') == 'in vivo'
    check('independent_continuous_daily_ratio', (reconstructed.reindex(a.index)[vivo] - a['median'][vivo]).abs().max() < 1e-9)
    rt = table('postprocess/route_validity.tsv.gz')
    means = ['mean_low_birth_rate', 'mean_division_rate', 'mean_ploidy']
    source = rt.group.isin(['direct_first_generation', 'historical_lineage'])
    structural = (source & (rt.q.eq(1) | rt.day.eq(0))) | (
        rt.initial_ploidy.eq(4) & (rt.group.eq('never_HL_lineage') | rt.q.eq(1) | rt.day.eq(0)))
    finite = np.isfinite(rt[means].to_numpy()).all(axis=1)
    valid = (~structural) & rt.ancestry_identifiable.fillna(False) & finite & np.isfinite(rt.mass_fraction) & rt.mass_fraction.ge(1e-8)
    check('independent_relative_only_mask', np.array_equal(valid.to_numpy(), rt.valid.to_numpy()))
    check('theoretical_empty_mask', np.array_equal(structural.to_numpy(), rt.structurally_empty.to_numpy()))
    check('unique_route_points', not rt.duplicated(key + ['group', 'route']).any())
    groupkey = key + ['group']
    eligibility = rt.groupby(groupkey, observed=True).valid.agg(['all', 'sum', 'size'])
    check('route_count_by_environment', all(r.size == (3 if r.Index[1] == 'in vitro' else 1)
          for r in eligibility.itertuples()))
    raw = table('tables/low_population_time_summary.tsv').set_index(groupkey + ['metric']).sort_index()
    shown = table('postprocess/low_population_time_summary.tsv').set_index(groupkey + ['metric']).sort_index()
    check('raw_derived_same_keys', raw.index.equals(shown.index))
    ploidy = shown.xs('mean_ploidy', level='metric')['median'].dropna()
    check('conditional_mean_ploidy_within_low_grid', ploidy.between(1 - 1e-10, 65 / 22 + 1e-10).all())
    for metric in ['mean_division_rate', 'mean_low_birth_rate']:
        check('nonnegative_' + metric, shown.xs(metric, level='metric')['median'].dropna().ge(-1e-12).all())
    for metric in means:
        rr = raw.xs(metric, level='metric')
        dd = shown.xs(metric, level='metric')
        support = eligibility.reindex(dd.index)['all']
        check('masked_' + metric, dd.loc[~support, ['median', 'q25', 'q75']].isna().all().all())
        check('retained_' + metric, np.allclose(dd.loc[support, ['median', 'q25', 'q75']],
              rr.loc[support, ['median', 'q25', 'q75']], atol=1e-12, rtol=0, equal_nan=True))
        check('eligible_route_count_' + metric, np.array_equal(dd.n_valid.to_numpy(), eligibility.reindex(dd.index)['sum'].to_numpy()))
        rows = rt[groupkey + [metric]].copy()
        rows.loc[~np.isfinite(rows[metric]), metric] = np.nan
        for p, column in [(.25, 'q25'), (.5, 'median'), (.75, 'q75')]:
            independent = rows.groupby(groupkey, observed=True)[metric].quantile(p).reindex(rr.index)
            check('independent_route_quantile_' + metric + '_' + column,
                  np.allclose(independent, rr[column], atol=1e-10, rtol=0, equal_nan=True))
    for metric in ['cell_share_percent', 'low_birth_share_percent']:
        check('nonmean_unchanged_' + metric, np.allclose(raw.xs(metric, level='metric')[['median', 'q25', 'q75']],
              shown.xs(metric, level='metric')[['median', 'q25', 'q75']], atol=1e-12, rtol=0, equal_nan=True))
    sens = table('postprocess/sensitivity_summary.tsv')
    check('relative_sensitivity_only', set(sens.relative_mass_min) == {1e-9, 1e-8, 1e-7} and not sens.absolute_count_filter.any())
    grid = table('provenance/grid_comparisons.tsv')
    check('all_grid_comparisons_pass', len(grid) > 0 and grid.passed.all())
    ev = table('provenance/passage_events.tsv')
    check('native_passage_no_fallback', not ev.uniform_fallback.any())
    check('1000_days_each_passage_route', ev.groupby(['task_id', 'q', 'route']).retained_days.sum().eq(1000).all())
    if args.output_root:
        import pypdfium2 as pdfium
        from PIL import Image, ImageDraw
        from pypdf import PdfReader
        renders = []
        specs = [('Figure7_fitted_p_misseg', 1, 14, 11), ('supp_fig7_8_fitted_p_misseg', 24, 24, 13)]
        for stem, pages, width, height in specs:
            pdf = args.output_root / 'Figures' / (stem + '.pdf')
            reader = PdfReader(pdf)
            check(stem + '_pages', len(reader.pages) == pages)
            doc = pdfium.PdfDocument(pdf)
            for index, page in enumerate(reader.pages):
                # Cairo/Helvetica encodes ordinary "fi" as its standard ligature.
                text = unicodedata.normalize('NFKC', page.extract_text() or '')
                check(stem + '_page_%02d_geometry' % (index + 1),
                      abs(float(page.mediabox.width) - width * 72) < 1 and abs(float(page.mediabox.height) - height * 72) < 1)
                check(stem + '_page_%02d_labels' % (index + 1), 'fitted p_misseg' in text and 'C01' in text and 'C02' in text)
                image = doc[index].render(scale=1).to_pil().convert('RGB')
                check(stem + '_page_%02d_nonblank' % (index + 1), np.asarray(image).std() > 8)
                image.save(args.audit_root / (stem + '_%02d.png' % (index + 1)))
                renders.append((stem + ' p%d' % (index + 1), image.copy()))
            for ext in ['pdf', 'png']:
                p = args.output_root / 'Figures' / (stem + '.' + ext)
                check(stem + '_' + ext + '_mirror', sha(p) == sha(args.output_root / 'manuscript/Figures' / p.name))
        for sheet in range((len(renders) + 7) // 8):
            canvas = Image.new('RGB', (1600, 1200), 'white')
            draw = ImageDraw.Draw(canvas)
            for i, (label, img) in enumerate(renders[sheet * 8:(sheet + 1) * 8]):
                img.thumbnail((790, 270));x = (i % 2) * 800;y = (i // 2) * 300
                draw.text((x + 8, y + 5), label, fill='black');canvas.paste(img, (x + 5, y + 25))
            canvas.save(args.audit_root / ('overview_%d.png' % (sheet + 1)))
    if args.protected_before:
        assert args.workspace
        rows = json.loads(args.protected_before.read_text())
        check('all_original_files_unchanged', all(sha(args.workspace / r['path']) == r['sha256'] for r in rows), len(rows))
    summary = dict(status='passed', checks=len(checks), fingerprint=done['fingerprint'],
                   daily_pairs=int(a['median'].notna().sum()), route_group_days=len(rt))
    (args.audit_root / 'result.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
