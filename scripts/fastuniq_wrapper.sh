#!/usr/bin/env bash
# USAGE: fastuniq_wrapper.sh <IN1> <IN2> <OUT1> <OUT2>
#
# Deduplicate paired end reads, taking fastq.gz input and outputting the same.

tmpdir=$(mktemp -d /tmp/tmp.XXXXXX)
pipe_r1=$tmpdir/$(basename $1).pipe
pipe_r2=$tmpdir/$(basename $2).pipe
mkfifo $pipe_r1
mkfifo $pipe_r2
gunzip < $1 > $pipe_r1 &
gunzip < $2 > $pipe_r2 &
fastuniq -i <(echo $pipe_r1 $pipe_r2 | tr ' ' '\n') -o >(gzip > $3) -p >(gzip > $4)
wait
