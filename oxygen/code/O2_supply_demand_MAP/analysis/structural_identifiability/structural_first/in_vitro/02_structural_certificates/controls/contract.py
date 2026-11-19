"""Finite control input/result contracts. Standard library only."""
import ast
import csv
import hashlib
from pathlib import Path
import re

CASES=('C02','C03','C04','C05')
F={'C02':'-a*x1','C03':'-(a+b)*x1','C04':'-a*x1**b','C05':'-sym.exp(a+b)*x1'}
OPTIONS="modelname='case_model'\ncheckObser=1\nmaxLietime=60\nnnzDerU=[]\nnnzDerW=[]\nprev_ident_pars=[]\n"

def model_text(case):
    if case not in CASES: raise ValueError('Out-of-scope case')
    params='[[a]]' if case=='C02' else '[[a],[b]]'
    return "import sympy as sym\nx1,a,b=sym.symbols('x1 a b',positive=True)\nx=[x1]\nh=[x1]\nf=["+F[case]+"]\np="+params+"\nu=[]\nw=[]\nvariables_locales=locals().copy()\n"

def table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields or list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def read(path):
    with Path(path).open() as stream: return list(csv.DictReader(stream,delimiter='\t'))

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def static_contract(case,model,options=OPTIONS):
    if case not in CASES: raise ValueError('Out-of-scope case')
    tree=ast.parse(model)
    fields={node.targets[0].id:node.value for node in tree.body if isinstance(node,ast.Assign) and isinstance(node.targets[0],ast.Name)}
    for name in ('x','h','f'):
        node=fields[name]
        if not isinstance(node,ast.List) or len(node.elts)!=1 or isinstance(node.elts[0],(ast.List,ast.Tuple)):
            raise ValueError('Invalid scalar-list field '+name)
    assert isinstance(fields['x'].elts[0],ast.Name) and fields['x'].elts[0].id=='x1'
    assert isinstance(fields['h'].elts[0],ast.Name) and fields['h'].elts[0].id=='x1'
    expected=['a'] if case=='C02' else ['a','b']
    actual=[]
    for element in fields['p'].elts:
        if not isinstance(element,ast.List) or len(element.elts)!=1 or not isinstance(element.elts[0],ast.Name):
            raise ValueError('Invalid parameter column')
        actual.append(element.elts[0].id)
    if actual!=expected or len(set(actual))!=len(actual): raise ValueError('Parameter role mismatch')
    expression=fields['f'].elts[0]
    if ast.dump(expression)!=ast.dump(ast.parse(F[case],mode='eval').body): raise ValueError('Equation changed')
    free={n.id for n in ast.walk(expression) if isinstance(n,ast.Name)}-{'sym'}
    if free!=set(['x1']+expected): raise ValueError('Free symbol mismatch')
    for field in ('u','w'):
        if ast.literal_eval(fields[field])!=[]: raise ValueError('Unexpected input')
    opts={n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(options).body if isinstance(n,ast.Assign)}
    if opts!={'modelname':'case_model','checkObser':1,'maxLietime':60,'nnzDerU':[],'nnzDerW':[],'prev_ident_pars':[]}:
        raise ValueError('Options mismatch')
    return True

def named_set(text,allowed):
    node=ast.parse(text.strip(),mode='eval').body
    if not isinstance(node,(ast.List,ast.Tuple)): raise ValueError('Expected named list')
    result=[]
    for item in node.elts:
        if not isinstance(item,ast.Name) or item.id not in allowed: raise ValueError('Unexpected parameter token')
        result.append(item.id)
    if len(set(result))!=len(result): raise ValueError('Duplicate result parameter')
    return set(result)

def receipt_valid(row,case,input_hash,run_id):
    return row.get('case')==case and row.get('input_sha256')==input_hash and row.get('run_id')==run_id and row.get('status')=='PASS'

def select_native_matrix(paths,final_order):
    """Select by authoritative native order; retain legitimate intermediate files."""
    if not isinstance(final_order,int) or not 1<=final_order<=4:raise ValueError('Invalid final order')
    expected='obs_ident_matrix_case_model_'+str(final_order)+'_Lie_deriv.txt'
    selected=[Path(p) for p in paths if Path(p).name==expected]
    if len(selected)!=1:raise ValueError('Missing or ambiguous final matrix')
    return selected[0]

def safe_matrix(text,sp,symbols):
    if not text.startswith('onx = '): raise ValueError('Unexpected native matrix prefix')
    root=ast.parse(text[6:],mode='eval').body
    def convert(node):
        if isinstance(node,ast.List): return [convert(e) for e in node.elts]
        if isinstance(node,ast.Constant) and isinstance(node.value,int): return sp.Integer(node.value)
        if isinstance(node,ast.Name) and node.id in symbols: return symbols[node.id]
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,(ast.USub,ast.UAdd)):
            return -convert(node.operand) if isinstance(node.op,ast.USub) else convert(node.operand)
        if isinstance(node,ast.BinOp):
            left,right=convert(node.left),convert(node.right)
            for typ,op in [(ast.Add,lambda a,b:a+b),(ast.Sub,lambda a,b:a-b),(ast.Mult,lambda a,b:a*b),(ast.Div,lambda a,b:a/b),(ast.Pow,lambda a,b:a**b)]:
                if isinstance(node.op,typ): return op(left,right)
        if isinstance(node,ast.Call) and isinstance(node.func,ast.Name) and node.func.id in ('exp','log') and len(node.args)==1 and not node.keywords:
            return getattr(sp,node.func.id)(convert(node.args[0]))
        raise ValueError('Unapproved symbolic token '+type(node).__name__)
    return sp.Matrix(convert(root))

def policy_findings(folder):
    folder=Path(folder);bad=[]
    for path in folder.glob('blocked*'):
        if path.stat().st_size: bad.append('policy denial: '+path.name)
    for path in folder.glob('trace*'):
        for line in path.read_text().splitlines():
            match=re.search(r'(?:openat\([^,]+,\s*"|open\("|execve\(")([^"\n]+)"',line)
            ret=re.search(r'=\s*(-?\d+)\s*$',line)
            if match and ret and int(ret.group(1))>=0:
                low=match.group(1).lower();base=Path(low).name
                if '.json' in low or '/sys/java/' in low or any(s in low for s in ('libjvm','libjava.','libv8','libmozjs')) or ('json' in base and '.so' in base) or base in ('java','node','nodejs','js') or low.endswith('.js'):
                    bad.append(line)
    module_file=folder/'imported_modules.tsv'
    if not module_file.exists(): bad.append('missing final module inventory')
    else:
        for row in read(module_file):
            if row['module'].split('.')[0] in ('json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer'):
                bad.append('forbidden loaded module '+row['module'])
    return bad
