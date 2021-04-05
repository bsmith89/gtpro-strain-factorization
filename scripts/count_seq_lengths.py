#!/usr/bin/env python3

import sys
from Bio.SeqIO import parse

if __name__ == "__main__":
    # print('seq_id', 'length', sep='\t')
    for rec in parse(sys.argv[1], "fasta"):
        print(rec.id, len(rec.seq), sep="\t")
