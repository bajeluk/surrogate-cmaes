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
#FAILED=(4 5 6 10 13 15 17 18 21 22 23 24 25 152 290 315 320 321 384 385 501 502 506 509 510 514 518 522 528 529 533 558 574 578 581 582 586 589 590 594 1003 1004 1005 2100 2103 2111 2175 2181 2182)
#FAILED=(4 6 25 528 529 533)
#FAILED=(4 6 510 528 533 1001 2100 2181 2182)
#FAILED=(4 528)
FAILED_LINE=`exp/metacentrum_list.sh -e $EXPID {1..192} {201..392} {501..692} {1001..1192} {2001..2192} | grep -A1 "Errored out" | tail -n1`
read -a FAILED <<< $FAILED_LINE

if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "No failed ids found."
  exit
fi

echo -n "Submit ${#FAILED[@]} failed ids (${FAILED[@]})? [y/N] "
read ANS

if [[ $ANS = [yY] ]]; then
  echo "Submitting..."
  sleep 1
else
  exit
fi

ID=1
QUEUE="04:00:00"
DIM=2
for FUNC in `seq 1 24`; do
  for IDS in $INPUT_IDS; do
    echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
    if [[ $? -eq 0 ]]; then
      subtask
    else
      ID=$((ID+1))
    fi
  done
done

ID=201
DIM=3
for FUNC in `seq 1 24`; do
  for IDS in $INPUT_IDS; do
    echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
    if [[ $? -eq 0 ]]; then
      subtask
    else
      ID=$((ID+1))
    fi
  done
done

ID=501
DIM=5
for FUNC in `seq 1 24`; do
  for IDS in $INPUT_IDS; do
    echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
    if [[ $? -eq 0 ]]; then
      subtask
    else
      ID=$((ID+1))
    fi
  done
done

ID=1001
QUEUE="24:00:00"
DIM=10
for FUNC in `seq 1 24`; do
  for IDS in $INPUT_IDS; do
    echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
    if [[ $? -eq 0 ]]; then
      subtask
    else
      ID=$((ID+1))
    fi
  done
done

ID=2001
QUEUE="48:00:00"
DIM=20
for FUNC in `seq 1 24`; do
  for IDS in $INPUT_IDS; do
    echo "${FAILED[@]}" | grep "\b$ID\b" > /dev/null
    if [[ $? -eq 0 ]]; then
      subtask
    else
      ID=$((ID+1))
    fi
  done
done
