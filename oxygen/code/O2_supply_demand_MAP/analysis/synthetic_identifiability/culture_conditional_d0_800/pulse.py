"""Lightweight login-node orchestration; no scientific calculations here."""
import fcntl,sys,time,traceback
from pathlib import Path
import driver
import workflow as w
s=Path(sys.argv[1]).resolve()
with (s/'state/CONTROLLER_LOCK').open('a') as handle:
    fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
    w.typed(s/'state/controller.tsv',dict(pid=__import__('os').getpid(),started_utc=w.utc(),release=str(Path(__file__).parent),scope='Only this isolated study'))
    # Bounded ten-day orchestration lifetime, independent of a laptop session.
    for tick in range(14400):
        if (s/'COMPLETE').exists():break
        try:
            with (s/'state/DRIVER_LOCK').open('a') as turn:
                fcntl.flock(turn,fcntl.LOCK_EX|fcntl.LOCK_NB);driver.advance(s)
        except Exception as e:
            traceback.print_exc()
            if not (s/'tables/STOP.tsv').exists():w.typed(s/'tables/STOP.tsv',dict(reason='controller_hold',message=str(e),utc=w.utc()))
        sys.stdout.flush();time.sleep(60)
    else:w.typed(s/'tables/STOP.tsv',dict(reason='controller_lifetime_exhausted',utc=w.utc()))
