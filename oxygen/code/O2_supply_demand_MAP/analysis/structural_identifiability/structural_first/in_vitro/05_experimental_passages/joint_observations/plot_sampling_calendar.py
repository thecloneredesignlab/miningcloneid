"""Draw the admitted culture sampling calendar as a static, self-contained SVG.

Run with Python 3 and redirect stdout to a file to reproduce the V5 diagram.
Inputs are the sealed joint design, the A7K crosswalk and the recorded passages.
"""
from collections import Counter, defaultdict
import csv
from datetime import date
from hashlib import sha256
from pathlib import Path
from xml.etree import ElementTree as ET

from paths import FIT, FROZEN_DESIGN, SYNTHETIC


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def calendar_day(stamp):
    return date.fromisoformat(stamp[:10])


def read_records():
    sources = rows(SYNTHETIC / "20261122_experimental_joint_observation_freeze_01" / "source_manifest.tsv")
    source, = [r for r in sources if r["source_id"] == "cloneid_passaging_snapshot"]
    raw_path = Path(source["locator"])
    assert sha256(raw_path.read_bytes()).hexdigest() == source["sha256"]
    raw = {r["id"]: r for r in rows(raw_path)}
    graph = {r["passage_id"]: r for r in rows(FROZEN_DESIGN / "passage_graph.tsv")}
    aliases = {}
    for r in rows(FIT / "a7k_event_contract.tsv"):
        observed = raw[r["harvest"]]
        assert observed["date"] == r["harvest_timestamp"] and observed["passaged_from_id1"] == r["seed"]
        aliases[r["harvest"].replace("SUM159_", "SUM-159_", 1)] = r["harvest"]
    events = {}
    for r in rows(FROZEN_DESIGN / "event_mask.tsv"):
        event, passage = r["event"], r["passage"]
        original = raw[aliases.get(event, event)]
        assert original["date"] == r["timestamp"] and original["passaged_from_id1"] == passage
        if event not in events:
            events[event] = (passage, calendar_day(r["timestamp"]), set())
        p, day, modalities = events[event]
        assert p == passage and day == calendar_day(r["timestamp"]) and r["modality"] not in modalities
        modalities.add(r["modality"])
    ends = {}
    for passage, day, _ in events.values():
        ends[passage] = max(ends.get(passage, day), day)
    for p, r in graph.items():
        donor = r["donor_passage"]
        if donor:
            draw = raw[r["donor_record"]]
            assert draw["date"] == r["donor_timestamp"] and draw["passaged_from_id1"] == donor
            assert raw[p]["passaged_from_id1"] == r["donor_record"]
            day = calendar_day(draw["date"])
            assert day <= calendar_day(raw[p]["date"])
            ends[donor] = max(ends.get(donor, day), day)
    assert set(ends) == set(graph)
    assert Counter(len(v[2]) for v in events.values()) == {1: 30, 2: 189, 3: 6}
    assert sum(not r["donor_passage"] for r in graph.values()) == 2
    assert sum(r["mixture_rule_free_mixtures"] == "free_simplex" for r in graph.values()) == 14
    return graph, raw, events, ends


def add(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {k.replace("_", "-"): str(v) for k, v in attrs.items()})


def label(parent, x, y, value, **attrs):
    add(parent, "text", x=x, y=y, **attrs).text = value


def line(parent, x1, y1, x2, y2, **attrs):
    add(parent, "line", x1=x1, y1=y1, x2=x2, y2=y2, **attrs)


def symbol(parent, shape, x, y, size, **attrs):
    if shape == "circle":
        return add(parent, "circle", cx=x, cy=y, r=size, **attrs)
    if shape == "square":
        return add(parent, "rect", x=x-size, y=y-size, width=2*size, height=2*size, **attrs)
    pts = ((x,y-size),(x+size,y),(x,y+size),(x-size,y)) if shape == "diamond" else (
        (x,y-size),(x+size,y+size*.65),(x-size,y+size*.65))
    return add(parent, "polygon", points=" ".join(f"{a:.2f},{b:.2f}" for a,b in pts), **attrs)


INK = "#253c50"
COHORT = {"2N": "#08785f", "4N": "#a94d0a"}
OXYGEN = {"20.5":"#154b86", "2":"#286ca4", "1":"#3e87b2", "0.5":"#589cc0",
          "0.3":"#71b0ce", "0.2":"#87bfd7", "0.1":"#9dcddb", "0":"#b5dce5"}
CHANNELS = {"corrected_count": "corrected count", "death_fraction": "death fraction", "karyotype": "karyotype"}


def draw():
    graph, raw, events, ends = read_records()
    terminal_passages = set(graph).difference(r["donor_passage"] for r in graph.values())
    last_harvest = {p: event for event,(p,day,_) in events.items()
                    if p in terminal_passages and day == ends[p]}
    assert set(last_harvest) == terminal_passages
    start = date(2024,11,12)
    def x(day):
        return round(151 + (day-start).days * 3.71, 2)
    def lane(p):
        cohort, branch, *_ = p.split("_")[2:]
        base = 0 if cohort == "2N" else 270
        if branch == "dp":
            return base + 138
        return base + {"C":103, "O1":170, "O2":248}[branch] + (35 if "_A7K_" in p else 0)

    svg = ET.Element("svg", {"xmlns":"http://www.w3.org/2000/svg","viewBox":"0 0 1400 592",
        "class":"sampling-calendar", "role":"img", "aria-labelledby":"sampling-title sampling-desc",
        "font-family":"Arial, Helvetica, sans-serif", "fill":INK})
    add(svg,"title",id="sampling-title").text = "Recorded joint-analysis harvests by date, lineage and assay"
    add(svg,"desc",id="sampling-desc").text = ("Each cohort starts at its dp seed and branches to control, O1 and O2. "
        "Passage strips follow calendar dates. Each harvest is marked by its admitted assays. "
        "The two squares mark the model starts. A vertical cap marks the last fitted harvest on each "
        "branch, including the A7K stock branches. Triangles mark the additional free-mixture restart roots.")
    # Keep both backgrounds on a single, unbroken calendar axis, including time in storage.
    ticks=[date(2024,12,1)]+[date(2025,m,1) for m in range(1,11)]
    for d in ticks:
        pos=x(d)
        line(svg,pos,44,pos,572,stroke="#dbe3e8",stroke_width=.9)
        label(svg,pos,35,d.strftime("%b"),font_size=14,text_anchor="middle",fill="#526578")
    label(svg,x(start),35,"12 Nov ’24",font_size=14,text_anchor="middle",fill="#526578")
    label(svg,x(date(2025,1,1)),16,"2025",font_size=13,text_anchor="middle",font_weight="bold")
    line(svg,12,306,1384,306,stroke="#ccd9e1")
    for c,offset in (("2N",0),("4N",270)):
        label(svg,17,79+offset,f"{c} lineage",font_size=20,font_weight="bold",fill=COHORT[c])
        for name,y in (("Control",103),("O1",170),("O2",248)):
            label(svg,113,y+5+offset,name,font_size=15,text_anchor="end")
        label(svg,133,143+offset,"Start · dp",font_size=14,text_anchor="end",fill=COHORT[c],font_weight="bold")
    # Retain physical stock-donor edges in both starting-composition designs.
    links=add(svg,"g",fill="none")
    bands=add(svg,"g")
    nodes=add(svg,"g")
    samples=add(svg,"g")
    for p,r in graph.items():
        day=calendar_day(raw[p]["date"])
        left,right=x(day),x(ends[p])
        y=lane(p)
        assert left<right
        donor=r["donor_passage"]
        if donor:
            dx,dy=x(calendar_day(r["donor_timestamp"])),lane(donor)
            if dy!=y:
                if "_A7K_" in p:
                    # The stock branch leaves the day-5 A7 donor, not the resumed main chain.
                    line(links,dx,dy,dx,y,stroke="#647a8a",stroke_width=1.4)
                    label(nodes,dx+5,dy+17,"A7 donor",font_size=11.5,fill="#546576")
                else:
                    line(links,dx,dy,dx,y,stroke=COHORT[r["cohort"]],stroke_width=1.6)
            if left>dx:
                line(links,dx,y,left,y,stroke="#758793",stroke_width=1.35,stroke_dasharray="2.8 3.6")
        band=add(bands,"rect",x=left,y=y-3.5,width=round(right-left,2),height=7,
                 fill=OXYGEN[r["oxygen_pct"]],stroke="white",stroke_width=.55,data_passage=p)
        add(band,"title").text=f"{p} · {day.isoformat()} to {ends[p].isoformat()} · nominal oxygen {r['oxygen_pct']}%"
        if r["mixture_rule_free_mixtures"]=="free_simplex":
            if donor:
                mark=symbol(nodes,"triangle",left,y+8.5,4.9,fill=COHORT[r["cohort"]],data_root=p)
                label(nodes,left,y+26,p.split("_")[4],font_size=12.5,text_anchor="middle",fill=COHORT[r["cohort"]])
            else:
                mark=symbol(nodes,"square",left,y,5.6,fill=COHORT[r["cohort"]],data_root=p,data_fit_start=p)
            add(mark,"title").text=f"{p} · {day.isoformat()} · "+("free-mixture restart" if donor else "model start (dp seed)")
        if p in last_harvest:
            final_day = ends[p]
            cap = add(nodes,"line",x1=right,y1=y-5,x2=right,y2=y+8,stroke=INK,
                      stroke_width=2.8,data_fit_end=last_harvest[p])
            add(cap,"title").text = f"Last fitted harvest · {last_harvest[p]} · {final_day.isoformat()}"
            passage_label = p.split("_")[4]
            # A7K already carries its date beside the three-readout diamond.
            if passage_label != "A7K":
                text = f"{passage_label} · {final_day.day} {final_day.strftime('%b')}"
                label(nodes,right,y-24,text,font_size=13,text_anchor="end",fill=INK)
    for event,(p,day,mods) in events.items():
        assert "corrected_count" in mods
        if len(mods)>1:
            assert "death_fraction" in mods
        pos,y=x(day),lane(p)-11.2
        if "karyotype" in mods:
            passage_label = p.split("_")[4]
            # Lift main-chain diamonds above neighboring harvests, keeping the exact date on the strip.
            diamond_y = y - (2 if passage_label == "A7K" else 15)
            line(samples,pos,diamond_y,pos,lane(p)-3.5,stroke="#97410c",stroke_width=1.1)
            marker=symbol(samples,"diamond",pos,diamond_y,12,fill="#b44909",stroke="white",stroke_width=1.4)
            sample_label = f"{passage_label} · {day.day} {day.strftime('%b')}"
            if passage_label == "A7K":
                label(samples,pos+17,diamond_y+4.5,sample_label,font_size=14,font_weight="bold",fill="#82400d")
            else:
                label(samples,pos,diamond_y-16,sample_label,font_size=14,font_weight="bold",text_anchor="middle",fill="#82400d")
        else:
            marker=symbol(samples,"circle",pos,y,3.2,fill=(INK if "death_fraction" in mods else "white"),
                          stroke=INK,stroke_width=1.4)
        marker.set("data-event",event)
        add(marker,"title").text=f"{event} · {day.isoformat()} · "+", ".join(CHANNELS[m] for m in CHANNELS if m in mods)
    return ET.tostring(svg,encoding="unicode")


if __name__ == "__main__":
    print(draw())
