"""Freeze native-comparison values of the preregistered mathematical witnesses."""
from pathlib import Path
import sys
from runtime import boot
from protocol import table

study=Path(sys.argv[1]);task=sys.argv[2]
sp,flint,np,out=boot(study,task)
from algebra import point,encode
flint.ctx.prec=256
native=[];balls=[]
for wid in ('W1','W2','W3'):
    p=point(study,wid)
    for name,value in p.items():
        native.append(dict(witness=wid,parameter=name,value=format(float(value.mid()),'.17g')))
        balls.append(dict(witness=wid,parameter=name,precision=256,**encode(value)))
target=study/'verification/B2';target.mkdir(exist_ok=True)
table(target/'native_witness_values.tsv',native)
table(target/'point_enclosures.tsv',balls)
table(out/'application_result.tsv',[dict(status='PASS',task=task,scope='mathematical certificate points; not culture observations or generating truths')])
print('REGISTERED_NATIVE_WITNESS_VALUES_READY',flush=True)
