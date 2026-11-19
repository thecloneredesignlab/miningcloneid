"""Standard-library immutable-source verification; no mathematical imports."""
from pathlib import Path
import sys
from common import verify
verify(Path(sys.argv[1]));print('IMMUTABLE_B6_INPUTS_AND_RELEASE_PASS')
