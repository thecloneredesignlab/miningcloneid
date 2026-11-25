"""Render Panel A from sealed three-readout evidence and verify the local HTML deliverable."""

from datetime import datetime
from fractions import Fraction
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import re
from urllib.parse import unquote, urlsplit
from zoneinfo import ZoneInfo

from scan import ARCHIVES, read, write
from v6_sides import RESOLUTION

WITNESSES = ARCHIVES/"20261122_in_vitro_joint_structural_steps04_07_01"
HISTORICAL = ARCHIVES.parent/"invitro_passage_structural_identifiability_v5.html"
DEST = HISTORICAL.with_name("invitro_passage_structural_identifiability_v6.html")


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = []
        self.hrefs = []
    def handle_starttag(self,tag,attrs):
        a = dict(attrs)
        if tag=="script" or any(k.startswith("on") for k in a):
            raise ValueError("Executable browser content found")
        if "id" in a:
            self.ids.append(a["id"])
        if "href" in a:
            self.hrefs.append(a["href"])


def natural(x):
    return format(float(x),".7g")


def side_label(s):
    st = s["status"]
    if st=="Statistical endpoint bracket":
        a,b=sorted((float(s["supported"]),float(s["excluded"])))
        return f'{natural(a)}–{natural(b)} <small>bracket</small>'
    if st.endswith("boundary"):
        return f'{natural(s["domain_bound"])} <small>{escape(st)}</small>'
    outer = "" if s["excluded"]=="" else f' · excluded {natural(s["excluded"])}'
    return f'<strong>Unresolved</strong><small>Supported {natural(s["supported"])}{outer}</small>'


def plot(rows,points,conf,table_path):
    width,height=980,615
    left,right=150,930
    def xx(x):
        return left+(right-left)*float(x)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img" aria-labelledby="v6-plot-title v6-plot-description" style="min-width:720px;width:100%;height:auto">',
           '<title id="v6-plot-title">Joint profile support and structural witnesses for 14 mechanisms</title>',
           '<desc id="v6-plot-description">Rows show the fitted estimate and sampled joint support on each registered parameter scale. Exact endpoint brackets and unresolved sides are tabulated below.</desc>',
           '<g font-family="system-ui,sans-serif" fill="#203844">']
    def line(x1,y1,x2,y2,color,weight=1,extra=""):
        svg.append(f'<line x1="{x1:.2f}" x2="{x2:.2f}" y1="{y1}" y2="{y2}" stroke="{color}" stroke-width="{weight}" {extra}/>')
    for tick in (0,.25,.5,.75,1):
        x=xx(tick)
        line(x,33,x,550,"#d7e0e2")
        svg.append(f'<text x="{x}" y="572" text-anchor="middle" font-size="13">{tick:g}</text>')
    records = []
    table = []
    witness = [{r["parameter"]:float(Fraction(r["value"])) for r in read(WITNESSES/f/"witness.tsv")}
               for f in ("inherit_512","inherit_384_second")]
    names = list(dict.fromkeys(s["parameter"] for s in rows))
    assert len(names)==14 and len(rows)==28
    for i,name in enumerate(names):
        y = 47+38*i
        pair = {s["side"]:s for s in rows if s["parameter"]==name}
        L,R = pair["Lower"],pair["Upper"]
        estimate = float(L["estimate"])
        def xvalue(v):
            return xx(conf.target_coordinate(name,float(v)))
        svg.append(f'<text x="{left-14}" y="{y+4}" text-anchor="end" font-size="14">{escape(name)}</text>')
        line(left,y,right,y,"#e1e7e8")
        for s in (L,R):
            sx = xvalue(s["supported"])
            cx = xvalue(estimate)
            pending = [float(v) for v in s["unresolved_points"].split()]
            evidence = [p for p in points if p["parameter"]==name and
                        min(cx,sx)-1e-6<=xvalue(p["value"])<=max(cx,sx)+1e-6]
            # Only join neighboring supported observations. An intermediate undecided
            # value breaks the line; resolved support can resume beyond it.
            marks = sorted([(cx,"Retained"),(sx,"Retained")]+[(xvalue(p["value"]),p["decision"]) for p in evidence])
            for (a,da),(b,db) in zip(marks,marks[1:]):
                if da==db=="Retained":
                    line(a,y,b,y,"#17636b",5)
            if s["status"].endswith("boundary"):
                line(sx,y-9,sx,y+9,"#134e55",3)
            elif s["status"]=="Statistical endpoint bracket":
                ex = xvalue(s["excluded"])
                line(sx,y,ex,y,"#a56519",6)
                for cap in (sx,ex):
                    line(cap,y-5,cap,y+5,"#74470e",1.5)
            elif s["status"]=="Unresolved":
                # Arrow tip points outward only as a status marker, never as an endpoint.
                k = -1 if s["side"]=="Lower" else 1
                svg.append(f'<path d="M {sx+k*4:.2f} {y-5} l {k*6} 5 l {-k*6} 5" fill="none" stroke="#7b4c68" stroke-width="2.2"/>')
            for v in pending:
                px=xvalue(v)
                svg.append(f'<circle cx="{px:.2f}" cy="{y}" r="2.5" fill="white" stroke="#7b4c68"/>')
        for p in (p for p in points if p["parameter"]==name):
            px=xvalue(p["value"])
            if p["decision"]=="Excluded":
                line(px-2.4,y+7,px+2.4,y+12,"#86512d",1.4)
                line(px+2.4,y+7,px-2.4,y+12,"#86512d",1.4)
            elif p["decision"]=="Retained":
                svg.append(f'<circle cx="{px:.2f}" cy="{y}" r="2.7" fill="#17636b"/>')
        for k,w in enumerate(witness):
            wx=xvalue(w[name])
            if k==0:
                svg.append(f'<circle cx="{wx:.2f}" cy="{y-11}" r="4" fill="white" stroke="#525c70" stroke-width="1.8"/>')
            else:
                svg.append(f'<rect x="{wx-2.3:.2f}" y="{y-13.3}" width="4.6" height="4.6" fill="#525c70"/>')
        fx=xvalue(estimate)
        svg.append(f'<path d="M {fx:.2f} {y-6} l 5.5 6 l -5.5 6 l -5.5 -6 Z" fill="#172f43" stroke="white" stroke-width="1"/>')
        table.append(f'<tr><th scope="row"><code>{escape(name)}</code></th><td>{natural(estimate)}</td><td>{side_label(L)}</td><td>{side_label(R)}</td></tr>')
        records.extend(dict(s,witness_1=witness[0][name],witness_2=witness[1][name],
                            registered_estimate=conf.target_coordinate(name,estimate)) for s in (L,R))
    svg+=['<text x="540" y="600" text-anchor="middle" font-size="14">Position on each parameter’s registered fitting scale (0–1)</text></g></svg>']
    write(table_path,records)
    return "".join(svg),"".join(table)


def render(out,rows,points,conf,final=False):
    stamp = datetime.now(ZoneInfo("America/New_York"))
    chart,body=plot(rows,points,conf,out/"panel_a_data.tsv")
    complete = sum(all(r["status"]!="Unresolved" for r in rows if r["parameter"]==n)
                   for n in set(r["parameter"] for r in rows))
    unresolved = sum(r["status"]=="Unresolved" for r in rows)
    stage = "Final review" if final else "Sealed evidence to date"
    version = f'{stage}. Both sides are characterized for {complete} of 14 parameters. {unresolved} of 28 sides remain unresolved.'
    prefix=escape(str(out.relative_to(ARCHIVES.parent)))
    panel = f'''<div id="panel-a-support"><h4>A. Three-readout joint profile support for all 14 mechanisms</h4>
<p><strong>{version}</strong> Each row uses corrected visible counts, measured death classifications and verified chromosome histograms under strict donor inheritance. The retained range describes the sampled component around the fit after refitting the other active parameters. Individual predictions can change within this range.</p>
<p class="small">Blue-green line, sampled support · dark diamond, fitted estimate · vertical cap, declared domain boundary · amber segment, endpoint bracket · open chevron, unresolved side. Filled dots mark retained values, crosses mark excluded values and open dots mark undecided values. The circle and square above each row mark structural witnesses 1 and 2.</p>
<div class="plot-scroll" role="region" aria-label="V6 joint profile intervals" tabindex="0">{chart}</div>
<p>Each fixed-value test uses the five-percent pointwise bootstrap decision and includes unresolved pairs in the tail bounds. Each statistical bracket spans at most {RESOLUTION:g} on the registered coordinate. The table reports original units. Boundaries mark limits of the declared model domain. Unresolved gaps remain open; coverage is parameter-wise.</p>
<div class="table-scroll wide"><table class="panel-support"><thead><tr><th>Mechanism</th><th>Joint estimate</th><th>Lower side</th><th>Upper side</th></tr></thead><tbody>{body}</tbody></table></div>
<p class="sources"><a href="{prefix}/panel_a_data.tsv">All 28 sides and structural witnesses</a> · <a href="{prefix}/tested_values.tsv">Tested values and sources</a> · <a href="{prefix}/endpoint_status.tsv">Endpoint search and stopping records</a></p></div>
'''
    source=HISTORICAL.read_text()
    start=source.index("<h4",source.index('<div class="joint-figure-body">',source.index('id="joint-witness-profile-figure"')))
    stop=source.index('<h4 id="joint-j-profile-panel">',start)
    doc=source[:start]+panel+source[stop:]
    doc=doc.replace("identifiability in culture | V5","identifiability in culture | V6")
    doc=doc.replace("Structural identifiability · V5","Structural identifiability · V6")
    doc=doc.replace("<footer>V5 ·","<footer>V6 ·")
    doc=re.sub(r'<time datetime="[^"]+">Evidence through [^<]+</time>',
               f'<time datetime="{stamp:%Y-%m-%d}">Evidence through {stamp:%d %B %Y}</time>',doc,count=1)
    doc=re.sub(r'(<footer>V6 · [^<]+? · )Evidence through [^<]+',rf'\g<1>Evidence through {stamp:%d %B %Y}',doc,count=1)
    doc=re.sub(r'<p><strong>V5 extends the fitted-data evidence\.</strong>.*?</p>',
               f'<p><strong>V6 maps the joint profile endpoints.</strong> {version} <a href="#panel-a-support">See Panel A.</a></p>',doc,count=1)
    doc=doc.replace("Panel A gives the tested values for each mechanism","Panel A links each joint profile to its endpoint evidence")
    doc=doc.replace("Constructing parameter-specific support intervals requires calibration at candidate endpoints along each profile.",
                    'The subsequent endpoint calibrations and their decisions appear in <a href="#panel-a-support">Panel A</a>.')
    doc=doc.replace('<a href="invitro_passage_structural_identifiability_v1.html">Original V1</a>',
       '<a href="invitro_passage_structural_identifiability_v5.html">V5 joint profiles</a> · <a href="invitro_passage_structural_identifiability_v1.html">Original V1</a>',1)
    links=Links()
    links.feed(doc)
    assert len(links.ids)==len(set(links.ids)) and len(rows)==28 and "C1" not in panel
    for href in links.hrefs:
        loc=urlsplit(href)
        if loc.scheme or loc.netloc:
            continue
        if loc.path:
            assert (DEST.parent/unquote(loc.path)).is_file(),href
        elif loc.fragment:
            assert unquote(loc.fragment) in links.ids,href
    tmp=DEST.with_suffix(".tmp")
    tmp.write_text(doc)
    tmp.replace(DEST)
