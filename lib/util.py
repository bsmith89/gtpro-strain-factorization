import sys
from datetime import datetime


def info(*msg):
    now = datetime.now()
    print(f"[{now}]", *msg, file=sys.stderr, flush=True)
