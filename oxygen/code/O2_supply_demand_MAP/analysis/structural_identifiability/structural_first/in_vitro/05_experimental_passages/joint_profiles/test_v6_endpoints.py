"""Exercise release selection and HTML rendering without simulation or model fitting."""

from pathlib import Path
import re
from tempfile import TemporaryDirectory
from unittest.mock import patch
from xml.etree import ElementTree

import v6_panel as panel
from v6_sides import PRIORITY, evaluate_side, queue
from scan import write


class CoordinateFixture:
    idx = {n:i for i,n in enumerate(PRIORITY)}

    def natural(self,b,n):
        return b[self.idx[n]]

    def target_coordinate(self,n,v):
        return float(v)


def check():
    conf=CoordinateFixture()
    base=[.5]*14
    reg=dict(parameter="p_misseg",domain_lower=0.,domain_upper=1.,lower_search="0.2",upper_search="0.75")
    def point(v,d):
        return dict(parameter="p_misseg",value=v,coordinate=v,decision=d)
    def lower(p):
        return evaluate_side(conf,base,reg,"Lower",p,[])
    assert lower([])["candidate"]==0.
    pts=[point(.19,"Excluded"),point(.205,"Retained"),point(.05,"Retained"),point(.85,"Inconclusive")]
    assert lower(pts)["status"]=="Statistical endpoint bracket" and lower(pts)["possible_reentry"]
    assert lower([point(0,"Retained")])["status"]=="Supported to domain boundary"
    assert lower([point(x,"Inconclusive") for x in (.1,.2,.3)])["stop_reason"]
    base[0]=0.
    assert lower([])["status"]=="Fitted domain boundary"
    base[0]=.5
    rows=[evaluate_side(conf,base,dict(reg,parameter=n),side,pts,[]) for n in PRIORITY for side in ("Lower","Upper")]
    assert len(queue(rows))==4
    with TemporaryDirectory() as temp:
        oxygen=Path(temp)/"oxygen"
        root=oxygen/"results"
        archive=root/"synthetic_identifiability"
        out=archive/"fixture"
        out.mkdir(parents=True)
        for f in ("endpoint_status.tsv","tested_values.tsv"):
            (out/f).touch()
        for f in panel.ARCHIVES.iterdir():
            (archive/f.name).symlink_to(f)
        for f in panel.HISTORICAL.parent.iterdir():
            if not (root/f.name).exists():
                (root/f.name).symlink_to(f)
        for f in panel.HISTORICAL.parent.parent.iterdir():
            if not (oxygen/f.name).exists():
                (oxygen/f.name).symlink_to(f)
        witnesses=Path(temp)/"witnesses"
        for directory,v in (("inherit_512",.45),("inherit_384_second",.47)):
            write(witnesses/directory/"witness.tsv",[dict(parameter=n,value=v) for n in PRIORITY])
        dest=root/"fixture_report.html"
        with patch.object(panel,"ARCHIVES",archive),patch.object(panel,"DEST",dest),patch.object(panel,"WITNESSES",witnesses):
            panel.render(out,rows,pts,conf)
        doc=dest.read_text()
        fig=doc[doc.index('<div id="panel-a-support">'):doc.index('<h4 id="joint-j-profile-panel">')]
        assert fig.count("<tr>")==15 and "C1" not in fig and "V1" not in fig
        assert "culture | V6</title>" in doc
        ElementTree.fromstring(re.search(r"<svg\b.*?</svg>",fig,re.S).group())
    print("Passed 28-side inventory, endpoint logic, re-entry, domain contact, report links and 14-row SVG/table tests.")


if __name__=="__main__":
    check()
