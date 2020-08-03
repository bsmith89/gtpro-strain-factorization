#!/bin/sh
# properties = {properties}

## Pre-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2

source ./env

{exec_job}
status=$?

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2

exit $status
