"""Current-chat baseline supervision; every allocation uses the finite driver.

This is a foreground host process, not an app automation or numerical worker.
It exits on a controller/transport defect or a finite baseline endpoint. A
successful baseline endpoint hands profile work back to this chat.
"""
from pathlib import Path
import argparse,subprocess,time,datetime

parser=argparse.ArgumentParser();parser.add_argument('study',type=Path)
s=parser.parse_args().study.resolve()
if s.name!='20261112_invivo_continuation_141206Z':raise ValueError('Wrong owned study')
remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+s.name
log=s/'reports/controller_passes.txt';previous=None
while True:
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ControlMaster=auto',
        '-o','ControlPath=/private/tmp/ivc_ssh_%C','-o','ControlPersist=600',
        '4482173@red.moffitt.org','python3',remote+'/policy_implementation/resolution01/driver.py',remote],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
    stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
    with log.open('a') as output:output.write(stamp+'\n'+result.stdout+result.stderr+'\n')
    if result.stdout!=previous:
        print(stamp+'\n'+result.stdout,flush=True);previous=result.stdout
    if result.returncode:
        print('CONTROLLER_HOLD exit='+str(result.returncode)+'\n'+result.stderr,flush=True)
        raise SystemExit(result.returncode)
    if 'BASELINE_TERMINAL' in result.stdout:
        print('FINITE_BASELINE_ENDPOINT; reconcile and close or implement conditional profiles in this chat.',flush=True)
        break
    time.sleep(45)
