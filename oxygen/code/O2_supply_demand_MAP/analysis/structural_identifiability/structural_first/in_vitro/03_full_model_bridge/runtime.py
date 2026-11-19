"""Guarded initialization before every transitive mathematical import."""
import atexit
import os
from pathlib import Path
import sys
from common import table

def boot(study,task):
    study=Path(study).resolve();out=study/'verification'/task
    out.mkdir(parents=True,exist_ok=True)
    assert not (out/'application_result.tsv').exists(),'Do not reuse a stale worker return'
    table(out/'bootstrap_receipt.tsv',[dict(stage='before_optional_imports',task=task,mathematical_imports_started=False)])
    forbidden={'json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','ipydex'}
    def allowed(path):
        if path is None or isinstance(path,int):return True
        p=Path(os.fsdecode(path)).resolve()
        return p==Path('/dev/null') or p==out or out in p.parents
    def deny(path):
        if not allowed(path):raise PermissionError('Scientific write outside phase output: '+str(path))
    def guard(event,args):
        if event=='import' and args[0].split('.')[0] in forbidden:raise ImportError('Excluded runtime dependency')
        if event=='open':
            path,mode,flags=args
            if flags is not None and (flags&os.O_ACCMODE!=os.O_RDONLY or flags&(os.O_CREAT|os.O_TRUNC|os.O_APPEND)):deny(path)
        if event in ('os.remove','os.rmdir','os.mkdir','os.chmod','os.truncate','os.utime'):deny(args[0])
        if event in ('os.rename','os.link','os.symlink'):
            deny(args[0]);deny(args[1])
    sys.addaudithook(guard)
    def inventory():
        table(out/'imported_modules.tsv',[dict(module=n,path=str(getattr(m,'__file__','builtin'))) for n,m in sorted(sys.modules.items())])
    atexit.register(inventory)
    # The guard and inventory are installed before path resolution/imports below.
    sys.path[:0]=[str(study/'runtime/private_prefix'),'/home/4482173/.local/lib/python3.11/site-packages',
                 *[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]
    os.environ['SYMPY_GROUND_TYPES']='python'
    import sympy
    import flint
    import numpy
    assert sympy.__version__=='1.12' and flint.__version__=='0.9.0'
    assert Path(flint.__file__).resolve().is_relative_to(study/'runtime/private_prefix')
    flint.ctx.threads=1
    table(out/'runtime_versions.tsv',[dict(package='SymPy',version=sympy.__version__,path=sympy.__file__),
        dict(package='python-flint',version=flint.__version__,path=flint.__file__),
        dict(package='FLINT',version=flint.__FLINT_VERSION__,path='bundled pinned native library'),
        dict(package='NumPy',version=numpy.__version__,path=numpy.__file__)])
    table(out/'bootstrap_receipt.tsv',[dict(stage='guarded_runtime_ready',task=task,mathematical_imports_started=True)])
    return sympy,flint,numpy,out
