#!/usr/bin/env bash
remote='pollard-bueno'
project='Projects/ucfmt'

for path in $@
do
    dirname=$(dirname $path)
    basename=$(basename $path)
    mkdir -p "$dirname"

    # Use rsync --partial for restartable transfers.
    partialdir=$TMPDIR/rsync-tmpdir/$project/$dirname
    mkdir -p "$partialdir"

    rsync -hravz -P --partial-dir="$partialdir" $remote:$project/$path "$dirname"
done
