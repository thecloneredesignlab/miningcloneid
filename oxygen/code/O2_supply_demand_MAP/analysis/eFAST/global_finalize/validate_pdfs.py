#!/usr/bin/env python3
"""Parse and render every exported PDF, reusing only hash-matched prior renders."""
import argparse
from pathlib import Path
import shutil
import subprocess

from PIL import Image,ImageStat,ImageDraw
from pypdf import PdfReader
from common import load,save,sha,now


def validate(root,pdftoppm):
    root=Path(root);out=root/'validation';out.mkdir(exist_ok=True)
    receipt=out/'pdf_render_validation.json'
    previous=load(receipt) if receipt.exists() else {'files':[]}
    prior={r['path']:r for r in previous['files']}
    rows=[];thumbnails=[];publication_thumbnails=[]
    paths=sorted(p for p in root.rglob('*.pdf') if 'validation' not in p.relative_to(root).parts)
    for path in paths:
        rel=str(path.relative_to(root));digest=sha(path)
        reader=PdfReader(str(path))
        if reader.is_encrypted or not reader.pages:raise ValueError('Unreadable/empty PDF '+rel)
        old=prior.get(rel)
        if old and old['sha256']==digest and len(old['pages'])==len(reader.pages):
            good=all((root/r['render_path']).exists() and sha(root/r['render_path'])==r['render_sha256'] for r in old['pages'])
        else:good=False
        pages=[]
        for page in range(len(reader.pages)):
            if good:
                r=old['pages'][page];render=root/r['render_path']
            else:
                stem=out/'pdf_renders'/path.relative_to(root).with_suffix('')
                stem=stem.parent/(stem.name+'_page%02d'%(page+1));stem.parent.mkdir(parents=True,exist_ok=True)
                subprocess.run([pdftoppm,'-f',str(page+1),'-l',str(page+1),'-singlefile','-scale-to','1400','-png',str(path),str(stem)],check=True,
                               stdout=subprocess.DEVNULL,stderr=subprocess.PIPE)
                render=Path(str(stem)+'.png')
                with Image.open(render) as image:
                    image.load();std=ImageStat.Stat(image.convert('L')).stddev[0]
                    if std<.5:raise ValueError('Blank PDF render '+rel)
                    size=image.size
                r=dict(page=page+1,render_path=str(render.relative_to(root)),render_sha256=sha(render),
                       image_size=size,grayscale_stddev=std)
            pages.append(r)
            with Image.open(render) as image:
                image=image.convert('RGB');image.thumbnail((430,540))
                thumb=Image.new('RGB',(450,575),'white')
                thumb.paste(image,((450-image.width)//2,25))
                ImageDraw.Draw(thumb).text((5,5),rel[-68:]+(' p%d'%(page+1)),fill='black')
                thumbnails.append(thumb)
                if rel.startswith('publication/'):publication_thumbnails.append(thumb.copy())
        rows.append(dict(path=rel,sha256=digest,n_pages=len(reader.pages),pages=pages,render_reused=good))
    sheet=Image.new('RGB',(450*3,575*((len(thumbnails)+2)//3)),'white')
    for i,thumb in enumerate(thumbnails):sheet.paste(thumb,((i%3)*450,(i//3)*575))
    sheet.save(out/'pdf_contact_sheet.png')
    for start in range(0,len(publication_thumbnails),12):
        selected=publication_thumbnails[start:start+12]
        preview=Image.new('RGB',(450*3,575*((len(selected)+2)//3)),'white')
        for i,thumb in enumerate(selected):preview.paste(thumb,((i%3)*450,(i//3)*575))
        preview.save(out/('publication_contact_sheet_%02d.png'%(start//12+1)))
    save(receipt,dict(verified_utc=now(),pdf_files=len(rows),rendered_pages=sum(r['n_pages'] for r in rows),
                     all_pdf_structures_and_nonblank_renders_passed=True,visual_review='separate manual review',files=rows))
    print('Validated',len(rows),'PDFs /',sum(r['n_pages'] for r in rows),'pages')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True);parser.add_argument('--pdftoppm',default=shutil.which('pdftoppm'))
    args=parser.parse_args()
    if not args.pdftoppm:parser.error('Provide a Poppler pdftoppm executable')
    validate(args.root,args.pdftoppm)
