#!/usr/bin/env python3
"""Exercise real CLI entry points, not mocked publishers. Writes only an audit tree."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import csv
import hashlib
import json
import subprocess
import time
from figure78_entry_tools import r_command

from PIL import Image, ImageChops, ImageStat
from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--audit-root', type=Path, required=True)
    parser.add_argument('--existing-full-root', type=Path)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--compare-pdf-raster', action='store_true',
                        help='Compare every canonical PDF page to production pixels using pypdfium2')
    args = parser.parse_args()
    workspace = args.workspace.resolve()
    audit = args.audit_root.resolve()
    audit.mkdir(parents=True, exist_ok=True)
    code = workspace / 'Code/Figures'
    rscript = r_command(workspace)
    results = []
    def snapshot():
        roots = ('Code/Figures','Code/config','Figures','manuscript/Figures',
                 'manuscript/figure_legends','data/Figures/Figure7','data/Figures/Figure8')
        paths = sorted({p for root in roots for p in (workspace/root).rglob('*')
                        if p.is_file() and p.name!='.DS_Store' and '__pycache__' not in p.parts})
        return {str(p):sha(p) for p in paths}
    def check(name, condition, detail=''):
        results.append({'check':name,'passed':bool(condition),'detail':str(detail)})
        print(('PASS' if condition else 'FAIL'),name,detail,flush=True)
        if not condition:
            save()
            raise AssertionError(name + ': ' + str(detail))
    def save():
        (audit/'checks.json').write_text(json.dumps(results,indent=2)+'\n')
    def sha(path):
        h = hashlib.sha256()
        with Path(path).open('rb') as handle:
            for chunk in iter(lambda:handle.read(8*1024*1024),b''):
                h.update(chunk)
        return h.hexdigest()
    def run(label, entry, options):
        command = rscript+['--vanilla',str(code/entry),*options]
        start = time.monotonic()
        completed = subprocess.run(command,cwd=workspace,capture_output=True,text=True)
        logs = audit/'logs';logs.mkdir(exist_ok=True)
        (logs/f'{label}.log').write_text('COMMAND '+json.dumps(command)+'\n'+completed.stdout+completed.stderr)
        return completed,time.monotonic()-start
    def components(stem):
        if stem not in matrices:return []
        return [f'{stem}_{context}_{family}_{initial}N' for context in ('invivo','invitro')
                for family in ('c01','c02') for initial in (4,2)]
    baseline = snapshot()
    (audit/'baseline.json').write_text(json.dumps(baseline,indent=2)+'\n')
    matrices = {'supp_fig7_3','supp_fig7_4','supp_fig7_3_ver2','supp_fig7_4_ver2'}
    stems = ['assembled_fig7']+[f'supp_fig7_{i}' for i in range(1,8)]+['figure7_ver2']+[
        f'supp_fig7_{i}_ver2' for i in range(1,8)]+['assembled_fig8']+[f'supp_fig8_{i}' for i in range(1,10)]
    full_root = args.existing_full_root.resolve() if args.existing_full_root else audit/'full'
    if not args.existing_full_root:
        completed,_ = run('full','draw_Figure78_all.R',[f'--output-root={full_root}'])
        check('full_native_cli_exit_zero',completed.returncode==0,completed.stderr[-500:])
    def one(stem):
        root = audit/'individual'/stem
        if stem=='assembled_fig7': entry='draw_Figure7.R'
        elif stem=='figure7_ver2': entry='draw_Figure7_ver2.R'
        elif stem=='assembled_fig8': entry='draw_Figure8.R'
        else:entry='draw_Supp_Figure'+stem[len('supp_fig'):]+'.R'
        options = [f'--output-root={root}','--formats=pdf']
        if stem in ('assembled_fig7','figure7_ver2','assembled_fig8'):options.append(f'--products={stem}')
        completed,elapsed = run(stem,entry,options)
        expected = {x+'.pdf' for x in [stem,*components(stem)]}
        found = {p.name for p in (root/'Figures').glob('*.pdf')}
        pages = len(PdfReader(root/'Figures'/f'{stem}.pdf').pages) if completed.returncode==0 else 0
        return stem,completed.returncode,found==expected,pages,elapsed,completed.stderr[-600:]
    with ThreadPoolExecutor(max_workers=max(1,args.workers)) as pool:
        futures = [pool.submit(one,stem) for stem in stems]
        for future in as_completed(futures):
            stem,exit_code,coverage,pages,elapsed,tail = future.result()
            check('individual_cli_'+stem,exit_code==0 and coverage and pages==(8 if stem in matrices else 1),
                  f'exit={exit_code}; pages={pages}; seconds={elapsed:.1f}; '+('' if exit_code==0 else tail))
    expected_all = {x+'.'+ext for s in stems for x in [s,*components(s)] for ext in ('pdf','png')}
    found_all = {p.name for p in (full_root/'Figures').iterdir() if p.suffix in ('.pdf','.png')}
    check('full_exact_export_coverage',found_all==expected_all,f'{len(found_all)} files per mirror')
    check('full_mirrors_byte_identical',all(sha(full_root/'Figures'/n)==sha(full_root/'manuscript/Figures'/n) for n in expected_all))
    pages_total = 0
    image_comparisons = []
    for stem in stems:
        reader = PdfReader(full_root/'Figures'/f'{stem}.pdf')
        expected_pages = 8 if stem in matrices else 1
        check('pdf_pages_'+stem,len(reader.pages)==expected_pages)
        pages_total += len(reader.pages)
        check('pdf_nonblank_'+stem,all(len(p.extract_text() or '')>100 and len(p.get_contents().get_data())>1000 for p in reader.pages))
        with Image.open(full_root/'Figures'/f'{stem}.png') as im:
            current = im.convert('RGB')
            check('png_nonblank_'+stem,ImageChops.difference(current,Image.new('RGB',current.size,'white')).getbbox() is not None)
            original_path = workspace/'Figures'/f'{stem}.png'
            with Image.open(original_path) as old:
                check('png_dimensions_unchanged_'+stem,current.size==old.size,f'{current.size}')
                current.thumbnail((1000,1000));reference=old.convert('RGB');reference.thumbnail((1000,1000))
                delta=ImageChops.difference(current,reference)
                mean=sum(ImageStat.Stat(delta).mean)/3
                image_comparisons.append({'stem':stem,'mean_absolute_channel_difference_0_255':mean})
                check('png_byte_identity_'+stem,sha(full_root/'Figures'/f'{stem}.png')==sha(original_path))
    check('all_component_png_bytes_identical',all(sha(full_root/'Figures'/n)==sha(workspace/'Figures'/n)
          for n in expected_all if n.endswith('.png')))
    check('canonical_pages_total' ,pages_total==54,pages_total)
    (audit/'image_comparisons.json').write_text(json.dumps(image_comparisons,indent=2)+'\n')
    if args.compare_pdf_raster:
        import pypdfium2 as pdfium
        pdf_comparisons = []
        for stem in stems:
            with pdfium.PdfDocument(full_root/'Figures'/f'{stem}.pdf') as current_pdf, \
                 pdfium.PdfDocument(workspace/'Figures'/f'{stem}.pdf') as original_pdf, \
                 pdfium.PdfDocument(audit/'individual'/stem/'Figures'/f'{stem}.pdf') as individual_pdf:
                check('pdf_reference_page_count_'+stem,len(current_pdf)==len(original_pdf))
                check('individual_reference_page_count_'+stem,len(current_pdf)==len(individual_pdf))
                for index in range(len(current_pdf)):
                    current=current_pdf[index].render(scale=96/72).to_pil().convert('RGB')
                    reference=original_pdf[index].render(scale=96/72).to_pil().convert('RGB')
                    identical=current.size==reference.size and ImageChops.difference(current,reference).getbbox() is None
                    pdf_comparisons.append({'stem':stem,'page':index+1,'pixels_identical':identical})
                    check(f'pdf_pixels_identical_{stem}_page{index+1}',identical)
                    single=individual_pdf[index].render(scale=96/72).to_pil().convert('RGB')
                    check(f'individual_pdf_pixels_identical_{stem}_page{index+1}',
                          current.size==single.size and ImageChops.difference(current,single).getbbox() is None)
        (audit/'pdf_visual_fidelity.json').write_text(json.dumps(pdf_comparisons,indent=2)+'\n')
    # Re-render just one PDF amid a full existing output set: all other bytes stay fixed.
    before = {str(p):sha(p) for folder in ('Figures','manuscript/Figures') for p in (full_root/folder).iterdir() if p.is_file()}
    completed,_ = run('single_in_full','draw_Supp_Figure7_1.R',[f'--output-root={full_root}','--formats=pdf'])
    check('single_in_full_exit_zero',completed.returncode==0)
    check('single_leaves_all_other_images_identical',all(sha(Path(p))==h for p,h in before.items() if Path(p).name!='supp_fig7_1.pdf'))
    # Missing/corrupt inputs are private fixtures, never changes to production data.
    negatives = audit/'negative_fixtures';negatives.mkdir(exist_ok=True)
    empty = negatives/'empty';empty.mkdir(exist_ok=True)
    bad = negatives/'corrupt';target=bad/'Figure7/tables/source_daughter_flux_binned.tsv'
    target.parent.mkdir(parents=True,exist_ok=True);target.write_text('corrupt fixture\n')
    forbidden = negatives/'arxives';forbidden.mkdir(exist_ok=True)
    cases = [
        ('missing_retained',[f'--data-root={empty}'],'Missing retained plotting input'),
        ('corrupt_retained',[f'--data-root={bad}'],'SHA-256 mismatch'),
        ('archive_root',[f'--data-root={forbidden}'],
         ('Archive access is forbidden','Old-version/archive roots are not Figure7/8 runtime inputs')),
        ('unknown_product',['--products=not_a_plot'],'Unknown, empty, or duplicated'),
        ('forbidden_output',[f'--output-root={workspace / "data/Figures"}'],'Output root must not'),
    ]
    for name,extra,message in cases:
        root = negatives/name
        completed,_ = run(name,'draw_Supp_Figure7_1.R',[f'--output-root={root}',*extra,'--formats=pdf'])
        messages = (message,) if isinstance(message,str) else message
        check('negative_'+name,completed.returncode!=0 and any(m in completed.stdout+completed.stderr for m in messages),
              (completed.stdout+completed.stderr)[-300:])
        check('negative_'+name+'_no_export',not (root/'Figures').exists())
    # Mutate only a validator-owned non-target fixture after a real renderer snapshots it.
    guard_root = negatives/'guard';(guard_root/'Figures').mkdir(parents=True,exist_ok=True)
    fixture = guard_root/'Figures/non_target_fixture.png';fixture.write_bytes(b'before')
    process = subprocess.Popen([*rscript,'--vanilla',str(code/'draw_Supp_Figure7_1.R'),
        f'--output-root={guard_root}','--formats=pdf'],cwd=workspace,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    deadline = time.monotonic()+120;changed=False
    while time.monotonic()<deadline and process.poll() is None:
        if list((guard_root/'audit').glob('**/protected_before.tsv')):
            fixture.write_bytes(b'changed by validator');changed=True;break
        time.sleep(.05)
    stdout,stderr = process.communicate(timeout=180)
    (audit/'logs/protection_negative.log').write_text(stdout+stderr)
    check('non_target_guard_detects_real_change',changed and process.returncode!=0 and 'non-target image' in stdout+stderr)
    check('guard_failure_does_not_publish_targets',not (guard_root/'Figures/supp_fig7_1.pdf').exists())
    # Current default group selections, without spending a second full raster pass.
    groups = [('draw_Figure7.R',8),('draw_Figure7_ver2.R',8),('draw_Figure8.R',10),
              ('draw_Figure7_versions.R',16),('draw_Figure8_time_snapshots.R',2),('draw_Figure8_family_time_snapshots.R',2)]
    for index,(entry,count) in enumerate(groups):
        root=audit/'default_groups'/str(index)
        completed,_=run('default_'+entry,entry,[f'--output-root={root}','--preflight-only=yes'])
        receipts=list(root.glob('audit/**/requested_products.tsv'))
        actual=len(list(csv.DictReader(receipts[0].open(),delimiter='\t'))) if len(receipts)==1 else 0
        check('default_group_'+entry,completed.returncode==0 and actual==count,f'{actual} products')
    # The repair-named compatibility entry is drawing-only too; exercise it
    # for each environment family without calling a repair/generation worker.
    for figure,stem in [('7','assembled_fig7'),('8','assembled_fig8')]:
        root=audit/'repair_alias'/figure
        completed,_=run('repair_alias_'+figure,'draw_Figure78_invitro_repair.R',
                       [f'--output-root={root}',f'--figure={figure}',f'--products={stem}','--formats=pdf'])
        check('repair_alias_'+figure,completed.returncode==0 and (root/'Figures'/f'{stem}.pdf').exists())
    after = snapshot()
    diffs=[p for p in set(baseline)|set(after) if baseline.get(p)!=after.get(p)]
    check('code_data_documents_and_published_figures_unchanged',not diffs,
          f'{len(baseline)} files; differences={diffs}')
    expected_inputs=list(csv.DictReader((workspace/'Code/config/manifests/figure78_current_plot_inputs.tsv').open(),delimiter='\t'))
    existing={str(p.relative_to(workspace/'data/Figures')) for fig in ('Figure7','Figure8')
              for p in (workspace/'data/Figures'/fig).rglob('*') if p.is_file() and p.name!='.DS_Store'}
    check('canonical_data_inventory_matches_manifest',existing=={r['relative_path'] for r in expected_inputs},f'{len(existing)} files')
    check('no_versioned_data_directories',all(not __import__('re').search(r'_v\d+(?:/|$)',p)
          for p in existing))
    # Data entries now perform real forward computation. A drawing-only
    # validator must explicitly use the no-simulation drawing preflight.
    for entry,count,extra in [('draw_Figure7.R',16,['--analysis=all']),('draw_Figure8.R',10,[])]:
        root=audit/'drawing_preflight'/entry
        completed,_=run('verify_'+entry,entry,[f'--output-root={root}','--preflight-only=yes',*extra])
        receipts=list(root.glob('audit/**/requested_products.tsv'))
        actual=len(list(csv.DictReader(receipts[0].open(),delimiter='\t'))) if len(receipts)==1 else 0
        check('canonical_drawing_preflight_'+entry,completed.returncode==0 and actual==count and not list((root/'Figures').glob('*')),
              f'exit={completed.returncode}; products={actual}')
    for selection,count in [('instantaneous',8),('history',8),('all',16)]:
        root=audit/'analysis_options'/selection
        completed,_=run('analysis_'+selection,'draw_Figure7.R',
          [f'--output-root={root}',f'--analysis={selection}','--preflight-only=yes'])
        receipts=list(root.glob('audit/**/requested_products.tsv'))
        actual=len(list(csv.DictReader(receipts[0].open(),delimiter='\t'))) if len(receipts)==1 else 0
        check('analysis_option_'+selection,completed.returncode==0 and actual==count)
    completions=list(full_root.glob('audit/**/completion.tsv'))
    receipts=[r for p in completions for r in csv.DictReader(p.open(),delimiter='\t')]
    check('receipts_record_zero_archive_reads_and_simulations',bool(receipts) and all(r['archive_reads']=='0' and r['new_simulations']=='0' for r in receipts))
    save()
    summary={'status':'passed','checks':len(results),'canonical_figures':26,'canonical_pdf_pages':54,
             'full_export_files_per_mirror':116,'archive_restorations':0,'new_simulations':0,
             'retained_inputs_unchanged':len(expected_inputs),'production_files_content_verified':len(baseline),
             'full_output_root':str(full_root),'audit_root':str(audit)}
    (audit/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    main()
