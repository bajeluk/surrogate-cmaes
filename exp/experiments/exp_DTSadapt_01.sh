#!/bin/bash
#
# Job-submitting experiment for the first 10D GP model testing
# It is expected to be located in exp/experiments/

export EXPID='exp_DTSadapt_01'
export DATASET="DTS_005_25_models"

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

INST="[1 2 3 4 5]"

ID=2
DIM=2
QUEUE="24:00:00"
FUNC="[`seq -s' ' 1 24`]"
subtask

ID=3
DIM=3
QUEUE="24:00:00"
FUNC="[`seq -s' ' 1 24`]"
subtask

ID=5
DIM=5
QUEUE="24:00:00"
FUNC="[`seq -s' ' 1 24`]"
subtask

ID=10
DIM=10
QUEUE="24:00:00"
FUNC="[`seq -s' ' 1 12`]"
subtask
FUNC="[`seq -s' ' 13 24`]"
subtask

ID=20
DIM=20
QUEUE="24:00:00"
FUNC="[`seq -s' ' 1 8`]"
subtask
FUNC="[`seq -s' ' 9 16`]"
subtask
FUNC="[`seq -s' ' 17 24`]"
subtask

exit 0


# Example of using submit_sequence
#
# # f1
# FUNC=1
# QUEUE="4:00:00"
# ID=490
# submit_sequence 146 1 146
# submit_sequence 147 2 170
# submit_sequence 217 2 256
