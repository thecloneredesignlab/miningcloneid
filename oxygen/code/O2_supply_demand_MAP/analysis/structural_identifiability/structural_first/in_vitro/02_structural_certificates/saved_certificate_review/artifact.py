"""Pure collection contracts and bounded verification; no scientific imports on import."""
import ast
import csv
import os
from pathlib import Path
import re
from base_contract import read,table,sha,named_set,safe_matrix,select_native_matrix,static_contract,policy_findings

def atomic_table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():raise ValueError('Immutable receipt already exists')
    temp=path.with_name(path.name+'.pending')
    with temp.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields or list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows);stream.flush();os.fsync(stream.fileno())
    os.replace(temp,path)

def start_call(path,case,run_id):
    if case not in ('C04','C05'):raise ValueError('Native case forbidden')
    path=Path(path)
    if (path.parent/'native_completion.tsv').exists():raise ValueError('Native completion already exists')
    with path.open('x',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t');writer.writerow(['case','run_id','native_calls']);writer.writerow([case,run_id,1]);stream.flush();os.fsync(stream.fileno())

def matrix_syntax(text,columns=3):
    if not text.startswith('onx = '):raise ValueError('Matrix prefix')
    node=ast.parse(text[6:],mode='eval').body
    if not isinstance(node,ast.List) or not node.elts:raise ValueError('Matrix rows')
    def term(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return
        if isinstance(n,ast.Name) and n.id in ('x1','a','b'):return
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):return term(n.operand)
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow)):
            term(n.left);term(n.right);return
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id in ('exp','log') and len(n.args)==1 and not n.keywords:return term(n.args[0])
        raise ValueError('Unapproved matrix token')
    for row in node.elts:
        if not isinstance(row,ast.List) or len(row.elts)!=columns:raise ValueError('Matrix column count')
        for cell in row.elts:term(cell)
    return len(node.elts),columns

def summary_sets(text):
    ident=re.search(r'These parameters are identifiable:\s*(\[[^\]]*\])',text)
    nonident=re.search(r'These parameters are unidentifiable:\s*(\[[^\]]*\])',text)
    if not ident or not nonident or 'Total execution time:' not in text:raise ValueError('Incomplete negative summary')
    return named_set(ident.group(1),{'a','b'}),named_set(nonident.group(1),{'a','b'})

def validate_metadata(m):
    required={'case','nd','rank','identifiable','nonidentifiable','measured','observable','native_returned','column_order','positive_real_symbols'}
    if not required.issubset(m):raise ValueError('Incomplete native metadata')
    if m['case'] not in ('C03','C04','C05') or not 1<=int(m['nd'])<=4 or int(m['rank']) not in (2,3):raise ValueError('Invalid metadata case/order/rank')
    if str(m['native_returned'])!='True' or str(m['positive_real_symbols'])!='True' or m['column_order']!='x1,a,b':raise ValueError('Invalid completion/domain/columns')
    for key in ('identifiable','nonidentifiable'):
        names=list(filter(None,str(m[key]).split(',')))
        if len(set(names))!=len(names) or not set(names)<={'a','b'}:raise ValueError('Invalid named metadata')
    if 'x1' not in (str(m['measured'])+','+str(m['observable'])).split(','):raise ValueError('State observation missing')
    return True

def select_matrix(folder,nd):
    folder=Path(folder)
    paths=list(folder.glob('obs_ident_matrix_*.txt'))
    if any(p.is_symlink() for p in paths):raise ValueError('Matrix symlink rejected')
    selected=select_native_matrix(paths,nd)
    if selected.resolve().parent!=folder.resolve():raise ValueError('Matrix outside artifact directory')
    return selected

def reconstruct_C03(folder,core_text):
    folder=Path(folder);log=(folder/'console.log').read_text()
    model_text=(folder/'native_app/models/case_model.py').read_text()
    static_contract('C03',model_text,(folder/'native_app/options.py').read_text())
    declarations=[n for n in ast.parse(model_text).body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Tuple)]
    if len(declarations)!=1 or not isinstance(declarations[0].value,ast.Call):raise ValueError('Saved symbol declaration')
    if not any(k.arg=='positive' and isinstance(k.value,ast.Constant) and k.value.value is True for k in declarations[0].value.keywords):raise ValueError('Saved positive-domain declaration')
    built=[int(v) for v in re.findall(r'matrix built with (\d+) Lie derivatives',log)]
    if built!=[2,3] or 'The model is structurally unidentifiable as a whole' not in log:raise ValueError('Ambiguous saved order evidence')
    if 'maxLietime has been reached' in log or 'More Lie derivatives would be needed' in log:raise ValueError('Saved derivative timeout')
    if 'Total execution time:' not in log or 'assert len(summaries)==len(matrices)==1' not in log:raise ValueError('Saved completion provenance')
    tree=ast.parse(core_text)
    branch=[n for n in ast.walk(tree) if isinstance(n,ast.If) and ast.dump(n.test)==ast.dump(ast.parse('nd >= len(xaug)',mode='eval').body)]
    if len(branch)!=1:raise ValueError('Unknown pinned stopping branch')
    for item in branch[0].body:
        for n in ast.walk(item):
            if isinstance(n,(ast.Assign,ast.AugAssign)):
                targets=n.targets if isinstance(n,ast.Assign) else [n.target]
                if any(isinstance(t,ast.Name) and t.id=='nd' for t in targets):raise ValueError('Stopping branch changes order')
    writer='"results/obs_ident_matrix_{}_{}_Lie_deriv.txt".format(options.modelname, nd)'
    if core_text.count(writer)<2:raise ValueError('Pinned final writer not found')
    summaries=list((folder/'native_app/results').glob('id_results_*.txt'))
    if len(summaries)!=1:raise ValueError('Ambiguous summary')
    ids,nonids=summary_sets(summaries[0].read_text())
    if ids or nonids!={'a','b'}:raise ValueError('Saved named result disagrees')
    matrices=list((folder/'native_app/results').glob('obs_ident_matrix_*.txt'))
    selected=select_matrix(folder/'native_app/results',3)
    if matrix_syntax(selected.read_text())!=(4,3):raise ValueError('Final matrix shape')
    prior=select_matrix(folder/'native_app/results',2)
    if matrix_syntax(prior.read_text())!=(3,3):raise ValueError('Intermediate matrix shape')
    calls=read(folder/'rank_calls.tsv')
    if len(calls)!=4 or any(r['rank']!='2' for r in calls):raise ValueError('Saved rank record inconsistency')
    return dict(case='C03',nd=3,rank=2,identifiable='',nonidentifiable='a,b',measured='x1',observable='',native_returned=True,column_order='x1,a,b',positive_real_symbols=True,order_provenance='reconstructed_from_saved_log_and_pinned_source',selected_matrix=str(selected),summary=str(summaries[0]))

def verify_matrix(case,metadata,matrix_path,sp,output,frame_path=None):
    validate_metadata(metadata)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    x,a,b=sp.symbols('x1 a b',positive=True);symbols=dict(x1=x,a=a,b=b)
    nd=int(metadata['nd']);assert 1<=nd<=4
    observed=safe_matrix(Path(matrix_path).read_text(),sp,symbols)
    assert observed.shape==(nd+1,3)
    if frame_path:
        assert observed==safe_matrix(Path(frame_path).read_text(),sp,symbols),'Native frame/file discrepancy'
    if case=='C04':
        jets=[x,-a*x**b,a*a*b*x**(2*b-1),-a**3*b*(2*b-1)*x**(3*b-2),a**4*b*(2*b-1)*(3*b-2)*x**(4*b-3)]
    else:
        k=a+b if case=='C03' else sp.exp(a+b)
        jets=[(-k)**i*x for i in range(5)]
    reference=sp.Matrix(jets[:nd+1]).jacobian(sp.Matrix([x,a,b]))
    rows=[]
    for i in range(nd+1):
        for j in range(3):
            residual=sp.simplify(observed[i,j]-reference[i,j]);assert residual==0,'Exact saved/native row mismatch'
            rows.append(dict(row=i,column=j,residual=str(residual)))
    table(output/'row_checks.tsv',rows)
    ids=set(filter(None,str(metadata['identifiable']).split(',')));nonids=set(filter(None,str(metadata['nonidentifiable']).split(',')))
    if case=='C04':
        assert ids=={'a','b'} and not nonids
        witness=-a*a*x**(3*b-1)
        assert sp.simplify(observed[:3,:].det()-witness)==0
        rank=3
    else:
        assert not ids and nonids=={'a','b'}
        assert all(sp.simplify(v)==0 for v in observed[:,1]-observed[:,2])
        witness=-x if case=='C03' else -sp.exp(a+b)*x
        assert sp.simplify(observed[:2,:2].det()-witness)==0
        c=sp.Symbol('c',real=True)
        f=-(a+b)*x if case=='C03' else -sp.exp(a+b)*x
        assert sp.simplify(f.subs({a:a+c,b:b-c},simultaneous=True)-f)==0
        rank=2
    assert int(metadata['rank'])==rank
    table(output/'rank_witness.tsv',[dict(rank=rank,witness=str(witness),domain='x1,a,b positive; symmetry shift -a<c<b',scope='toy control only')])
    return dict(case=case,status='PASS',rank=rank,derivative_order=nd,identifiable=','.join(sorted(ids)),nonidentifiable=','.join(sorted(nonids)),scope='saved_native_verified' if case=='C03' else 'native_toy_control_verified')
