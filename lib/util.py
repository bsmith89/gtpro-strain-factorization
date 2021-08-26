import sys
from datetime import datetime
from functools import reduce


def info(*msg):
    now = datetime.now()
    print(f"[{now}]", *msg, file=sys.stderr, flush=True)


def _mul(x, y):
    return x * y


def is_empty(x):
    return reduce(_mul, x.shape, 1) == 0
