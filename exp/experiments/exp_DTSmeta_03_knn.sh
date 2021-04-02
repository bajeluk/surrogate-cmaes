#!/bin/bash
#
# Job-submitting experiment for generation
# of snapshots from multiple DTS experiments, TSS knn.

export EXPID='exp_DTSmeta_03_knn'
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

DIMS=(2 3 5 10 20)
FUNCS=({1..24})
INSTS=(11 12 13 14 15)
# queues
Q2D="168:00:00"
Q3D="168:00:00"
Q5D="168:00:00"
Q10D="336:00:00"
QND="336:00:00"

IDS="$*"

# no IDS as input -> run default: all settings
if [ -z "${IDS}" ]; then

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
      for IN in ${INSTS[@]}; do
        INST="[${IN}]"
        # echo subtask ID=$ID DIM=$DIM FUNC=$FUNC INST=$INST QUEUE=$QUEUE
        subtask 
      done
    done
  done

# IDS given as input
else

  for ID in ${IDS[@]}; do
  
    # dimension
    if [ $[$ID%120] -eq 0 ]; then
      DIM=${DIMS[ $[$ID/120 - 1] ]}
    else
      DIM=${DIMS[ $[$ID/120] ]}
    fi
    # function
    if [ $[$ID%120%5] -eq 0 ]; then
      FUNC=${FUNCS[ $[$ID%120/5 - 1] ]}
    else
      FUNC=${FUNCS[ $[$ID%120/5] ]}
    fi
    # instance number
    INST="[${INSTS[ $[$ID%5 - 1] ]}]"

    # queue according to dimension
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

    # submit task
    echo subtask ID=$ID:  DIM=$DIM FUNC=$FUNC INST=$INST QUEUE=$QUEUE
    subtask

  done

fi
