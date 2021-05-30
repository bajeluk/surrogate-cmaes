#!/bin/bash
#
# Job-submitting experiment for metafeature calculation
# of TSS knn data selected from dataset DTS_meta_005_train.mat.
#
# It is expected to be located in exp/experiments/

export EXPID='exp_DTSmodels_meta_03_mfts_train_knn'
export DATASET="exp/experiments/dataset/DTS_meta_005_train.mat"

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
MODEL_IDS=1
DIMS=(2 3 5 10 20)
FUNCS=($(seq 1 24))
# queues
Q2D="4:00:00"
Q3D="4:00:00"
Q5D="24:00:00"
Q10D="24:00:00"
QND="24:00:00" # other dimensions (e.g., 20D)

INST="[11 12 13 14 15]"
IDS="[1 2 3 4 5 6 7]" # IDS in DTS_meta_005_train were selected at random from DTS_meta_005

ID=1

for DIM in ${DIMS[@]}; do
  case $DIM in
    2)
      QUEUE=$Q2D
      ;;
    3)
      QUEUE=$Q3D
      ;;
    5)
      QUEUE=$Q5D
      ;;
    10)
      QUEUE=$Q10D
      ;;
    *)
      QUEUE=$QND
      ;;
  esac

  for FUNC in ${FUNCS[@]}; do
    for MODEL_ID in ${MODEL_IDS[@]}; do
      OPTS="struct(%modelOptionsIndices%|$MODEL_ID)"

      #echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
      #echo "${MISSING[@]}" | grep "\bf${FUNC}_${DIM}D\b" > /dev/null
      #if [[ $? -eq 0 ]]; then
        subtask
      #else
      #  ID=$((ID+1))
      #fi
    done
  done
done
