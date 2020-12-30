#!/bin/sh
# properties = {properties}

## Pre-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2

# # This *should* prevent CUDA from
# # using GPUs that were not assigned.
# # FIXME: How to be sure?
# # FIXME: Why is SGE_GPU not set??
# export CUDA_VISIBLE_DEVICES=$SGE_GPU
# echo CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES SGE_GPU=$SGE_GPU

source ./env

{exec_job}
status=$?

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" >&2

env

exit $status
