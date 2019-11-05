#!/bin/bash
#
# Job-submitting experiment example for generation
# of snapshots from multiple DTS experiments.
# The resulting experiment .sh file is expected to be located 
# in exp/experiments.
#
# usage: $ exp/experiments/exp_modeltestsets_example.sh
#
# options:
#   -n | --dry-run     no submition will happen, only print what should be done
#
# see also:
#   exp/bash_settings.sh
#   exp/metacentrum_modeltestsets_common.sh

export EXPID='exp_modeltestsets_example'
export EXPLOG='exp_doubleEC_28_log_nonadapt'
export DATASET="exp/experiments/dataset/DTS_meta_005.mat"

# Enable this option for using Matlab MCR compilated binaries:
export useMCR=1

# !!! THIS HAS TO BE HERE !!!
#     vvvvvvvvvvvvvvvvvvv
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
. $CWD/../metacentrum_modeltestsets_common.sh
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

#for i in {1..6} {8..9}; do
#  OPTS[${#OPTS[@]}]="struct(%modelOptionsIndices%| $i)"
#done
#MODEL_IDS="1"

#INST="[11 12 13 14 15]"
#IDS="[1 2 3 4 5 6 8 9]" # 7 is ADD kernel

ID=1

for DIM in 2 3; do
  case $DIM in
    [235])
      QUEUE="04:00:00"
      ;;
    10)
      QUEUE="24:00:00"
      ;;
    *)
      QUEUE="48:00:00"
      ;;
  esac

  for FUNC in `seq 6 7`; do
    for IN in `seq 11 12`; do
      INST="[${IN}]"
      subtask
    done
  done
done
