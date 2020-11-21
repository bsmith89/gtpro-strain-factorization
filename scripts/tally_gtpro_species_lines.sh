#!/usr/bin/env bash

library=$1
bzcat "$2" \
    | cut -f1 \
    | uniq -c \
    | awk -v OFS='\t' -v library=$library '{print library,$2,$1}'
