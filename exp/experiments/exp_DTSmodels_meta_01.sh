#!/bin/bash
#
# Job-submitting experiment for the first 10D GP model testing
# It is expected to be located in exp/experiments/

export EXPID='exp_DTSmodels_meta_01'
export DATASET="DTS_meta_004"

# Enable this option for using Matlab MCR compilated binaries:
export useMCR=1

# !!! THIS HAS TO BE HERE !!!
#     vvvvvvvvvvvvvvvvvvv
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
. $CWD/../metacentrum_testmodels_common.sh
#     ^^^^^^^^^^^^^^^^^^^
# !!! THIS HAS TO BE HERE !!!

# job submittion can be done via one of the two bash functions:
#
# (a) subtask [JOBNAME_SUFFIX]
#
#     submits a job with current-set $ID, $DIM, $FUNC, $INST and $OPTS
#
#
# (b) submit_sequence LOW_IDX STEP HIGH_IDX
#
#     call subtask() with 'modelOptionsIndices' set to a small
#     numbers of model-settings;
#     it submits a job for every $2 consecutive modelOptions' indices
#
#     e.g.: submit 11 4 22
#     submits jobs for the following modelOption indices:
#       - 11:14
#       - 15:18
#       - 19:22
#
#     Important: the HIGH_IDX should be the equal to (k*LOW_IDX) - 1
#                for some k


# critical characters has to be replaced in $OPTS:
# '|' with ',' and "%" with "'"
OPTS=""

INST="[11 12 13 14 15]"
INPUT_IDS="1 2 3 4 5 6 8 9" # 7 is ADD kernel

ID=1
QUEUE="04:00:00"
for DIM in 2; do
  for FUNC in `seq 1 24`; do
    for IDS in $INPUT_IDS; do
      subtask
      exit
    done
  done
done

#ID=201
#QUEUE="24:00:00"
#for DIM in 3; do
#  for FUNC in `seq 1 24`; do
#    for IDS in $INPUT_IDS; do
#      subtask
#    done
#  done
#done
#
#ID=301
#QUEUE="48:00:00"
#for DIM in 5 10 20; do
#  for FUNC in `seq 1 24`; do
#    for IDS in $INPUT_IDS; do
#      subtask
#    done
#  done
#done
