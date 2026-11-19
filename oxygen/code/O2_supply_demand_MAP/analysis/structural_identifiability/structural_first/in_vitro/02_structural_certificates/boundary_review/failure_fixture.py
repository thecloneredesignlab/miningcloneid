"""One expected post-bootstrap missing-dependency receipt control."""
from pathlib import Path
import sys
from common import table
from runtime import boot

def main():
    study=Path(sys.argv[1]);task=sys.argv[2]
    sp,flint,np,out=boot(study,task)
    try:__import__('b6_intentionally_absent_fixture_dependency')
    except ModuleNotFoundError as exc:
        table(out/'expected_failure.tsv',[dict(expected=True,exception=type(exc).__name__,application_exit=1)])
        return 1
    raise AssertionError('Missing-dependency fixture unexpectedly resolved')

if __name__=='__main__':raise SystemExit(main())
