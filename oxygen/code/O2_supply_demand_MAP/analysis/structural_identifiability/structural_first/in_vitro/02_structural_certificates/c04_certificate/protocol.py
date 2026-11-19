"""Metadata-only contracts for a saved-matrix C04 certificate."""
import ast
import csv
import hashlib
from pathlib import Path

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def matrix_tree(text):
    if not text.startswith('onx = '):raise ValueError('Unexpected matrix format')
    tree=ast.parse(text[6:],mode='eval').body
    if not isinstance(tree,ast.List) or len(tree.elts)!=3:raise ValueError('Need exactly three saved rows')
    def term(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return
        if isinstance(n,ast.Name) and n.id in ('x1','a','b'):return
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):term(n.operand);return
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):term(n.left);term(n.right);return
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='log' and len(n.args)==1 and not n.keywords:term(n.args[0]);return
        raise ValueError('Unapproved saved-matrix token')
    for row in tree.elts:
        if not isinstance(row,ast.List) or len(row.elts)!=3:raise ValueError('Saved column count')
        for entry in row.elts:term(entry)
    return tree

def parse_matrix(text,sp,symbols):
    tree=matrix_tree(text)
    def value(n):
        if isinstance(n,ast.List):return [value(v) for v in n.elts]
        if isinstance(n,ast.Constant):return sp.Integer(n.value)
        if isinstance(n,ast.Name):return symbols[n.id]
        if isinstance(n,ast.UnaryOp):return -value(n.operand) if isinstance(n.op,ast.USub) else value(n.operand)
        if isinstance(n,ast.Call):return sp.log(value(n.args[0]))
        a,b=value(n.left),value(n.right)
        if isinstance(n.op,ast.Add):return a+b
        if isinstance(n.op,ast.Sub):return a-b
        if isinstance(n.op,ast.Mult):return a*b
        if isinstance(n.op,ast.Div):return a/b
        if isinstance(n.op,ast.Pow):return a**b
        raise ValueError('Unknown operator')
    return sp.Matrix(value(tree))
