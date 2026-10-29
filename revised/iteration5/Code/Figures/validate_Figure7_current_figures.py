#!/usr/bin/env python3
"""Inspect current Figure7 exports only; never simulates or redraws any figure."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageStat
from pypdf import PdfReader
import pypdfium2 as pdfium


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--audit-root', type=Path, required=True)
    args = parser.parse_args()
    workspace, root, audit = (p.resolve() for p in (args.workspace,args.output_root,args.audit_root))
    if workspace/'audit' not in audit.parents or workspace/'audit' not in root.parents:
        raise ValueError('Inspect staged exports under workspace/audit; do not overwrite production')
    audit.mkdir(parents=True,exist_ok=True)
    stems = ['assembled_fig7']+['supp_fig7_%d'%i for i in range(1,9)]
    matrices = {'supp_fig7_3','supp_fig7_4'}
    components = [s+'_'+ctx+'_'+family+'_'+initial+'N' for s in sorted(matrices)
                  for ctx in ('invivo','invitro') for family in ('c01','c02') for initial in ('4','2')]
    expected = {s+'.'+ext for s in stems+components for ext in ('pdf','png')}
    checks=[];records=[];raster_pages={}

    def check(name, passed, detail=''):
        checks.append(dict(test=name,passed=bool(passed),detail=detail))
        if not passed:
            (audit/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
            raise AssertionError(name+': '+str(detail))

    for folder in ('Figures','manuscript/Figures'):
        found={p.name for p in (root/folder).iterdir() if p.suffix in ('.pdf','.png')}
        check('exact_current_export_scope_'+folder,found==expected,dict(found=len(found),expected=50))
    check('publication_mirrors_identical',all(sha(root/'Figures'/n)==sha(root/'manuscript/Figures'/n) for n in expected))
    total=0
    for stem in stems+components:
        path=root/'Figures'/(stem+'.pdf')
        reader=PdfReader(path)
        wanted=48 if stem=='supp_fig7_8' else 8 if stem in matrices else 1
        check('page_count_'+stem,len(reader.pages)==wanted,len(reader.pages))
        total+=len(reader.pages)
        baseline=workspace/'Figures'/(stem+'.pdf')
        if baseline.is_file():
            previous=PdfReader(baseline)
            check('existing_page_count_preserved_'+stem,len(previous.pages)==len(reader.pages))
            check('existing_page_geometry_preserved_'+stem,
                  all(tuple(a.mediabox)==tuple(b.mediabox) for a,b in zip(previous.pages,reader.pages)))
        page_root=audit/'pages'/stem;page_root.mkdir(parents=True,exist_ok=True)
        images=[];texts=[];vector_marks=[]
        with pdfium.PdfDocument(path) as doc:
            for i in range(len(doc)):
                page=doc[i]
                # Dense heatmaps can exceed pypdf's stream-decompression cap.
                # Use the already-required native renderer for Unicode text
                # and path presence; never relax a global decompression limit.
                textpage=page.get_textpage()
                texts.append(textpage.get_text_bounded(errors='strict'))
                textpage.close()
                paths=page.get_objects(filter=[pdfium.raw.FPDF_PAGEOBJ_PATH])
                vector_marks.append(next(paths,None) is not None)
                del paths
                bitmap=page.render(scale=1)
                im=bitmap.to_pil().convert('RGB')
                check('raster_nonblank_'+stem+'_'+str(i+1),
                      max(ImageStat.Stat(im).stddev)>2)
                raster=page_root/('page_%02d.png'%(i+1));im.save(raster)
                images.append(raster)
                bitmap.close();page.close()
        check('pdf_text_and_marks_present_'+stem,all(len(t)>100 and marks
                                                    for t,marks in zip(texts,vector_marks)))
        check('no_replacement_glyph_'+stem,all('\ufffd' not in t and '\x00' not in t for t in texts))
        if stem in ['assembled_fig7','supp_fig7_3','supp_fig7_4','supp_fig7_5','supp_fig7_6','supp_fig7_7','supp_fig7_8']:
            check('blocking_label_beta_not_q_'+stem,all('q =' not in t for t in texts))
        raster_pages[stem]=[str(p.relative_to(audit)) for p in images]
        png=root/'Figures'/(stem+'.png')
        with Image.open(png) as image:
            rgb=image.convert('RGB')
            check('png_nonblank_'+stem,ImageChops.difference(rgb,Image.new('RGB',rgb.size,'white')).getbbox() is not None)
            baseline_png=workspace/'Figures'/(stem+'.png')
            if baseline_png.is_file():
                with Image.open(baseline_png) as old:
                    check('existing_png_dimensions_preserved_'+stem,image.size==old.size,[image.size,old.size])
        for ext in ('pdf','png'):
            p=root/'Figures'/(stem+'.'+ext)
            records.append(dict(path=str(p),sha256=sha(p),size=p.stat().st_size,product=stem,format=ext))
        # Overview sheets aid human inspection; individual full-page rasters
        # above remain the authority for readable labels and clipping checks.
        for start in range(0,len(images),12):
            subset=images[start:start+12]
            sheet=Image.new('RGB',(1200,((len(subset)+2)//3)*260),'#eeeeee')
            draw=ImageDraw.Draw(sheet)
            for j,p in enumerate(subset):
                with Image.open(p) as image:
                    image.thumbnail((390,235));x=(j%3)*400;y=(j//3)*260
                    sheet.paste(image,(x,y+22));draw.text((x+5,y+4),stem+' page '+str(start+j+1),fill='black')
            sheet.save(audit/(stem+'_overview_%02d.png'%(start//12+1)))
    check('all_current_PDF_pages',total==86,total)
    (audit/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    (audit/'raster_pages.json').write_text(json.dumps(raster_pages,indent=2)+'\n')
    with (audit/'artifact_checksums.tsv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,records[0].keys(),delimiter='\t');writer.writeheader();writer.writerows(records)
    result=dict(status='automated_current_figure_checks_passed',checks=len(checks),PDF_pages=total,
                artifacts_per_mirror=50,visual_review_of_page_rasters_required=True,
                old_pixels_expected_identical=False,simulations=0,redraws=0,Figure8_or_history_writes=False)
    (audit/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':
    main()
