from pathlib import Path
import sys
from base_contract import policy_findings,table
folder=Path(sys.argv[1]);bad=policy_findings(folder)
table(folder/'policy_gate.tsv',[dict(status='FAIL' if bad else 'PASS',finding_count=len(bad),detail='; '.join(bad))])
if bad:raise SystemExit(77)
