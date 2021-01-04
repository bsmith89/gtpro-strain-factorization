#!/usr/bin/env python3

from Bio.AlignIO import read
from tqdm import tqdm
import sys
import numpy as np

from tqdm import TqdmSynchronisationWarning
import warnings

warnings.simplefilter("ignore", TqdmSynchronisationWarning)


def contains_only(col, values):
    return set(col) - set(values) == set()


def main():
    alignment = read(sys.stdin, "fasta")
    remove_chars = np.asarray(list(sys.argv[1]))
    length = alignment.get_alignment_length()
    print("Num sequences: %d" % len(alignment), file=sys.stderr)
    print("Alignment length: %d" % length, file=sys.stderr)

    keep_cols = []
    for i in tqdm(range(length)):
        if not contains_only(alignment[:, i], remove_chars):
            keep_cols.append(i)

    out = np.empty((len(alignment), len(keep_cols)), dtype="<U1")
    for i, j in tqdm(enumerate(keep_cols), total=len(keep_cols)):
        out[:, i] = np.array(list(alignment[:, j]))

    print("Remaining columns: %d" % out.shape[1], file=sys.stderr)
    for i, seq in enumerate(alignment):
        print(">" + seq.id)
        print("".join(out[i, :]))


if __name__ == "__main__":
    main()
