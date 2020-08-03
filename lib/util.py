import sys
from datetime import datetime


def info(*msg):
    now = datetime.now()
    print(f'[{now}]', *msg, file=sys.stderr)


def idxwhere(condition):
    return list(condition[condition].index)


def normalize_rows(df):
    return df.divide(df.sum(1), axis=0)
