#!/usr/bin/env python3
"""Independent saved-value/48-page checks for the conditional-mean filter."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import pypdfium2 as pdfium
from PIL import Image, ImageChops, ImageDraw


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--output-root', type=Path, required=True)
    p.add_argument('--overlay-root', type=Path, required=True)
    p.add_argument('--audit-root', type=Path, required=True)
    p.add_argument('--before-pdf', type=Path, required=True)
    a = p.parse_args()
    a.audit_root.mkdir(parents=True, exist_ok=False)
    checks = []

    def check(name, passed, detail=None):
        checks.append(dict(test=name, passed=bool(passed), detail=detail))
        (a.audit_root/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
        if not passed:
            raise RuntimeError(name)

    key_columns = ['pair_label', 'model_context', 'O2_pct', 'p_misseg',
                   'threshold_N', 'initial_ploidy', 'q', 'group', 'day']

    def key(row):
        return tuple(row[k] for k in key_columns)

    means = ['mean_low_birth_rate', 'mean_division_rate', 'mean_ploidy']
    actual = {}
    route_rows = 0
    zero_rates = 0
    valid_subcell_rows = 0
    with gzip.open(a.overlay_root/'route_validity.tsv.gz', 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            mass = float(row['mass_fraction']) if row['mass_fraction'] != 'NA' else math.nan
            logn = float(row['log_total_cells']) if row['log_total_cells'] != 'NA' else math.nan
            rates = [float(row[m]) if row[m] != 'NA' else math.nan for m in means]
            source = row['group'] in {'direct_first_generation', 'historical_lineage'}
            beta, day, initial = float(row['q']), int(row['day']), float(row['initial_ploidy'])
            structural = (source and (beta == 1 or day == 0)) or (initial == 4 and
                (row['group'] == 'never_HL_lineage' or beta == 1 or day == 0))
            valid = (not structural and row['ancestry_identifiable'] == 'TRUE'
                     and all(math.isfinite(x) for x in rates)
                     and math.isfinite(mass) and mass >= 1e-8)
            if valid != (row['valid'] == 'TRUE'):
                raise RuntimeError('independent route validity differs: '+str(key(row)))
            k = key(row)
            n, good = actual.get(k, (0, 0))
            actual[k] = (n+1, good+int(valid))
            route_rows += 1
            zero_rates += int(valid and (rates[0] == 0 or rates[1] == 0))
            log_expected = (logn + math.log(mass) if math.isfinite(logn) and math.isfinite(mass) and mass > 0
                            else (-math.inf if math.isfinite(logn) and mass == 0 else math.nan))
            lt1 = None if math.isnan(log_expected) else log_expected < 0
            stored_lt1 = None if row['expected_cells_lt1'] == 'NA' else row['expected_cells_lt1'] == 'TRUE'
            if lt1 != stored_lt1:
                raise RuntimeError('independent expected-count flag differs: '+str(key(row)))
            valid_subcell_rows += int(valid and lt1 is True)
    check('independent_relative_only_thresholds_all_route_points', route_rows > 0, route_rows)
    check('independent_lt1_flags_do_not_mask_valid_means', valid_subcell_rows > 0, valid_subcell_rows)
    summary = {}
    with gzip.open(a.overlay_root/'summary_validity.tsv.gz', 'rt') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            k = key(row)
            n, good = actual[k]
            required = 3 if row['model_context'] == 'in vitro' else 1
            valid = n == required and good == required
            if (n != int(row['n_routes']) or good != int(row['n_eligible'])
                    or valid != (row['summary_valid'] == 'TRUE')):
                raise RuntimeError('route aggregation differs: '+str(k))
            summary[k] = valid
    check('independent_all_route_aggregation', summary.keys() == actual.keys(), len(summary))
    mean_rows = 0
    masked_rows = 0
    with (a.overlay_root/'low_population_time_summary.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if row['metric'] not in means:
                continue
            valid = summary[key(row)]
            missing = all(row[n] == 'NA' for n in ('median', 'q25', 'q75'))
            if valid == missing:
                raise RuntimeError('derived mean mask differs: '+str(key(row)))
            mean_rows += 1
            masked_rows += int(missing)
    check('independent_derived_mean_masks', mean_rows == 3*len(summary),
          dict(mean_rows=mean_rows, masked_rows=masked_rows, valid_zero_rate_route_points=zero_rates))
    names = {'supp_fig7_8.pdf', 'supp_fig7_8.png'}
    for folder in ('Figures', 'manuscript/Figures'):
        found = {f.name for f in (a.output_root/folder).iterdir() if f.is_file()}
        check('exact_export_scope_'+folder, found == names, sorted(found))
    check('publication_mirrors_identical', all(sha(a.output_root/'Figures'/name) ==
          sha(a.output_root/'manuscript/Figures'/name) for name in names))
    pdf = a.output_root/'Figures/supp_fig7_8.pdf'
    page_folder = a.audit_root/'pages'
    page_folder.mkdir()
    images = []
    with pdfium.PdfDocument(pdf) as doc, pdfium.PdfDocument(a.before_pdf) as old:
        check('page_count_48', len(doc) == len(old) == 48)
        for i in range(48):
            page = doc[i]
            check('page_geometry_'+str(i+1), page.get_size() == old[i].get_size())
            text = page.get_textpage().get_text_range()
            check('readable_labels_'+str(i+1), len(text)>100 and '\ufffd' not in text and '\x00' not in text)
            im = page.render(scale=1).to_pil().convert('RGB')
            check('nonblank_page_'+str(i+1), ImageChops.difference(im, Image.new('RGB', im.size, 'white')).getbbox() is not None)
            im.save(page_folder/f'page_{i+1:02d}.png')
            if i % 6 >= 3:
                previous = old[i].render(scale=1).to_pil().convert('RGB')
                check('unchanged_CDF_pixels_'+str(i+1), im.size == previous.size and
                      ImageChops.difference(im, previous).getbbox() is None)
            images.append(im)
    for offset in range(0, 48, 12):
        canvas = Image.new('RGB', (1728, 4*332), 'white')
        draw = ImageDraw.Draw(canvas)
        for j, im in enumerate(images[offset:offset+12]):
            thumb = im.copy()
            thumb.thumbnail((576, 312))
            x, y = j%3*576, j//3*332
            draw.text((x+5, y+3), f'Page {offset+j+1}', fill='black')
            canvas.paste(thumb, (x, y+20))
        canvas.save(a.audit_root/f'overview_{offset//12+1}.png')
    preview = a.output_root/'Figures/supp_fig7_8.png'
    with Image.open(preview) as im, Image.open(a.before_pdf.with_suffix('.png')) as old:
        check('preview_dimensions_unchanged', im.size == old.size, im.size)
    (a.audit_root/'result.json').write_text(json.dumps(dict(status='automated_checks_passed',
        checks=len(checks), route_points=route_rows, group_days=len(summary), mean_rows=mean_rows,
        masked_mean_rows=masked_rows, PDF_pages=48, unchanged_CDF_pages=24,
        valid_expected_cells_lt1_route_points=valid_subcell_rows, absolute_count_filter=False,
        visual_review_required=True), indent=2)+'\n')
    print('PASS:', len(checks), 'checks; 48 pages; 24 CDF pages unchanged')


if __name__ == '__main__':
    main()
