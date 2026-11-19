"""Scoped scientific runtime guard, used only in allocated workers."""
import atexit
import os
from pathlib import Path
import sys
from protocol import table

def boot(study,task):
    study=Path(study);out=study/'verification'/task;out.mkdir(parents=True,exist_ok=True)
    assert not (out/'application_result.tsv').exists(),'No stale worker return may be reused'
    blocked={'json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','ipydex'}
    def guard(event,args):
        if event=='import' and args[0].split('.')[0] in blocked:
            with (out/'blocked_imports.tsv').open('a') as f:f.write('BLOCKED\t'+args[0]+'\n')
            raise ImportError('Excluded scientific runtime dependency')
    sys.addaudithook(guard)
    def inventory():
        table(out/'imported_modules.tsv',[dict(module=n,path=str(getattr(m,'__file__','builtin'))) for n,m in sorted(sys.modules.items())])
    atexit.register(inventory)
    os.environ['SYMPY_GROUND_TYPES']='python'
    sys.path[:0]=[str(study/'runtime/private_prefix'),'/home/4482173/.local/lib/python3.11/site-packages',
                 *[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]
    import sympy
    import flint
    import numpy
    assert sympy.__version__=='1.12' and flint.__version__=='0.9.0'
    flint.ctx.threads=1
    table(out/'runtime_versions.tsv',[dict(package='SymPy',version=sympy.__version__,path=sympy.__file__),
          dict(package='python-flint',version=flint.__version__,path=flint.__file__),
          dict(package='FLINT',version=flint.__FLINT_VERSION__,path='private bundled native dependency'),
          dict(package='NumPy',version=numpy.__version__,path=numpy.__file__)])
    return sympy,flint,numpy,out
