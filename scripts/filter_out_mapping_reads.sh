#!/usr/bin/env bash
#
# USAGE:
#     filter_out_bowtie2_mapped_reads.sh <index_name> <fastq_in_r1> <fastq_in_r2> <fastq_out_r1> <fastq_out_r2>
#
# <index_name> is not the full path to the bowtie index;
# instead it's the short name, which usually excludes extensions
# like .1.bt2, .2.bt2, .rev.1.bt2, etc.
#
# -f 12 refers to flags read unmapped (0x4) + mate unmapped (0x8)

set -euo pipefail

bowtie2 --threads $1 --mm -x "$2" -1 "$3" -2 "$4" \
    | samtools view -f 12 -O BAM \
    | samtools fastq -1 "$5" -2 "$6" -
