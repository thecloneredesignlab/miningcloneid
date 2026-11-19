"""Additional local static rejection checks; no scientific imports."""
import ast
from pathlib import Path
import sys
from contract import CASES,OPTIONS,model_text,static_contract,table
study=Path(sys.argv[1]);checks=[]
def positive_declaration(text):
    assignments=[n for n in ast.parse(text).body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Tuple)]
    assert len(assignments)==1
    call=assignments[0].value
    assert isinstance(call,ast.Call) and isinstance(call.func,ast.Attribute) and call.func.attr=='symbols'
    assert ast.literal_eval(call.args[0])=='x1 a b'
    assert any(k.arg=='positive' and ast.literal_eval(k.value) is True for k in call.keywords)
for case in CASES:
    text=model_text(case);positive_declaration(text)
    checks.append(dict(check=case+':symbol_declaration',passed=True))
    mutations=[('missing_parameter',text.replace('p=[[a]]' if case=='C02' else 'p=[[a],[b]]','p=[]')),
               ('duplicate_parameter',text.replace('p=[[a]]' if case=='C02' else 'p=[[a],[b]]','p=[[a],[a]]')),
               ('numeric_parameter',text.replace('p=[[a]]' if case=='C02' else 'p=[[a],[b]]','p=[[1]]'))]
    if case=='C04':mutations.append(('fixed_exponent',text.replace('x1**b','x1**2')))
    for label,bad in mutations:
        try:static_contract(case,bad,OPTIONS)
        except (AssertionError,ValueError,KeyError):checks.append(dict(check=case+':reject_'+label,passed=True))
        else:raise AssertionError(label)
    try:positive_declaration(text.replace('positive=True','positive=False'))
    except AssertionError:checks.append(dict(check=case+':reject_lost_domain',passed=True))
    else:raise AssertionError('Lost domain admitted')
for path in Path(__file__).parent.glob('*.py'):ast.parse(path.read_text())
table(study/'reports/extra_static_checks.tsv',checks)
print('Additional static checks passed:',len(checks),'; all adapter Python sources parse.')
