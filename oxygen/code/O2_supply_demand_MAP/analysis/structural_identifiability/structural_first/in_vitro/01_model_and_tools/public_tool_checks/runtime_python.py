"""Dependency/import policy gate only; never calls a structural solver."""
import os
from pathlib import Path
import sys

study=Path(sys.argv[1])
log=study/'runtime/python/blocked_imports.tsv'
modules_log=study/'runtime/python/imported_modules.tsv'
forbidden={'json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer'}

def audit(event,args):
    if event=='import' and args[0].split('.')[0]=='ipydex' and os.environ.get('CULTURE_DISABLE_DEBUGGER')=='1':
        raise ImportError('Optional interactive debugger disabled in the private runtime')
    if event=='import' and args[0].split('.')[0] in forbidden:
        with log.open('a') as stream:
            stream.write('BLOCKED\timport\t'+args[0]+'\n')
        raise ImportError('Analysis policy excludes this module: '+args[0])

sys.addaudithook(audit)
sys.path[:0]=[
    str(study/('runtime/private_prefix' if os.environ.get('CULTURE_DISABLE_DEBUGGER')=='1' else 'sources/dependencies/symbtools/src')),
    *([] if os.environ.get('CULTURE_DISABLE_DEBUGGER')=='1' else [str(study/'sources/dependencies/ipydex/src')]),
    str(study/'sources/StrikePy/StrikePy'),
    '/home/4482173/.local/lib/python3.11/site-packages',
    *[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]
try:
    import numpy
    import sympy
    print('NUMPY',numpy.__version__,numpy.__file__,flush=True)
    print('SYMPY',sympy.__version__,sympy.__file__,flush=True)
    import symbtools
    import strike_goldd
    assert callable(strike_goldd.strike_goldd)
    print('SYMBTOOLS',symbtools.__version__,flush=True)
    x=sympy.Symbol('audit_x',positive=True)
    assert sympy.simplify(sympy.diff(sympy.exp(x),x)-sympy.exp(x))==0
    print('PYTHON_RUNTIME_GATE_COMPLETED',flush=True)
finally:
    with modules_log.open('w') as stream:
        stream.write('module\tpath\n')
        for name,module in sorted(sys.modules.items()):
            stream.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')
