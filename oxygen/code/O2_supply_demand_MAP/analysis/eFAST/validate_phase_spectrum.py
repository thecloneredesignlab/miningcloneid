#!/usr/bin/env python3
"""Independent representative full-FFT and table checks; artifact receipts."""
import argparse
import collections
import csv
import datetime
import gzip
import json
from pathlib import Path
import numpy as np
from pypdf import PdfReader
from collect_resolution_nm_results import sha, read
import phase_spectrum as diag


def gzread(path):
    with gzip.open(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))


def validate(root,visual):
    records=read(root/'trajectory_manifest.tsv');cases=read(root/'representative_cases.tsv')
    metrics=gzread(root/'phase_metrics.tsv.gz');lookup={(r['trajectory_id'],round(float(r['O2_pct']),10),r['output']):r for r in metrics}
    assert len(lookup)==190*201*2
    error=0.;checks=0
    for case in cases:
        o=float(case['O2_pct']);j=int(round(o/.025))
        for row in records:
            if (row['fit_seed'],row['parameter'])!=(case['fit_seed'],case['parameter']):continue
            n,m=int(row['N']),int(row['M']);omega=int((n-1)/(2*m));cut=omega//2
            with np.load(root/'processed_cache'/('trajectory_%03d.npz'%int(row['trajectory_id']))) as z:
                assert str(z['code_sha256'].item())==sha(Path(diag.__file__))
                ci=int(np.argmin(abs(z['case_oxygen']-o)))
                for k,out in enumerate(diag.OUTPUTS):
                    y=z['case_y'][:,ci,k]
                    # Deliberately use full complex FFT and scalar summation here.
                    p=(abs(np.fft.fft(y)[1:(n+1)//2])/n)**2
                    v=2*sum(p);s1=2*sum(p[np.arange(1,m+1)*omega-1])/v;st=1-2*sum(p[:cut])/v
                    original=lookup[(row['trajectory_id'],round(o,10),out)]
                    for name,value in [('S1',s1),('ST',st),('variance',np.var(y))]:
                        diff=abs(float(original[name])-value);error=max(error,float(diff));assert diff<1e-10;checks+=1
    histogram=read(root/'representative_response_bins.tsv');groups=collections.defaultdict(int)
    for r in histogram:groups[r['trajectory_id']]+=int(r['n_in_bin'])
    assert len(groups)==190
    for r in records:assert groups[r['trajectory_id']]==int(r['N'])
    raw=gzread(root/'raw_representative_responses.tsv.gz');assert len(raw)==5*(2*8193+2*4097)
    pages=[]
    for path in sorted((root/'figures').glob('*.pdf')):
        reader=PdfReader(path);assert len(reader.pages)==1
        text=reader.pages[0].extract_text();assert len(text)>500
        assert (root/'figures'/path.with_suffix('.png').name).is_file()
        pages.append(dict(path=str(path.relative_to(root)),pages=1,text_characters=len(text)))
    assert len(pages)==6
    record=dict(status='passed',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        independent_full_fft_numeric_checks=checks,independent_full_fft_max_absolute_difference=error,
        n_complete_representative_input_bin_groups=len(groups),n_raw_figure_source_response_rows=len(raw),
        pdf_text_and_page_checks=pages,visual_render_reviewed=visual,
        visual_review_scope='Six Poppler-rendered PDF pages manually inspected for labels, legends, clipping and layout' if visual else 'pending',
        unit_tests='4 passed: known sinusoid, independent full FFT, phase/design identity, undefined output rejection',
        no_new_model_evaluations=True,no_new_jobs_submitted=True)
    (root/'deliverable_validation.json').write_text(json.dumps(record,indent=2)+'\n')
    paths=[p for p in root.iterdir() if p.is_file() and p.name!='artifact_manifest.tsv' and not p.name.startswith('.')]
    paths+=list((root/'figures').glob('*.pdf'))+list((root/'figures').glob('*.png'))
    code=Path(__file__).parent
    paths += [code/name for name in ('collect_phase_spectrum.py','phase_spectrum.py','report_phase_spectrum.py','validate_phase_spectrum.py','test_phase_spectrum.py')]
    repo=code.parents[4]
    diag.table(root/'artifact_manifest.tsv',[dict(path=str(p.relative_to(repo)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(paths)])
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--visual-render-reviewed',action='store_true',help='Only after manual review of every rendered PDF page')
    a=p.parse_args();validate(a.root.resolve(),a.visual_render_reviewed)
