from pathlib import Path
import sys
from common import verify
verify(Path(sys.argv[1]));print('IMMUTABLE_SOURCE_PASS')
