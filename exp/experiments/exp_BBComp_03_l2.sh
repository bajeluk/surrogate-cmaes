#!/bin/bash
#
# Job-submitting experiment for the first 10D GP model testing
# It is expected to be located in exp/experiments/

export EXPID='exp_BBComp_03'

# Enable this option for using Matlab MCR compilated binaries:
export useMCR=1

# !!! THIS HAS TO BE HERE !!!
#     vvvvvvvvvvvvvvvvvvv
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
. $CWD/../metacentrum_bbcomp_common.sh
#     ^^^^^^^^^^^^^^^^^^^
# !!! THIS HAS TO BE HERE !!!

# job submittion can be done via the bash function:
#
# submit [JOBNAME_SUFFIX]
#
#     submits a job with current-set $INST and $OPTS

# critical characters has to be replaced in $OPTS:
# '|' with ',' and "%" with "'"

minutes_to_hours() {
  H=$(($1/60)); M=$(($1%60));
  TIME=`printf '%02d%02d' $((HOURS+H)) $((MINS+M))`
}

# 1-week jobs:

OPTS=""
MEMORY=""

QUEUE="168:00:00"
for INST in `seq 566 600` `seq 666 700` `seq 766 800`; do
  submit
done

MEMORY="4048mb"
QUEUE="168:00:00"
for INST in `seq 866 900` `seq 966 1000`; do
  submit
done

# 1-day jobs with increasing start time:

HOURS=1
MINS=00
PLUSMINS=0

MEMORY=""
QUEUE="24:00:00"
for INST in `seq 66 100` `seq 166 200` `seq 266 300`; do
  minutes_to_hours $PLUSMINS
  SUBMIT_PBS_PARAMS="-a $TIME"
  submit
  PLUSMINS=$((PLUSMINS+1))
done

QUEUE="48:00:00"
for INST in `seq 366 400` `seq 466 500`; do
  minutes_to_hours $PLUSMINS
  SUBMIT_PBS_PARAMS="-a $TIME"
  submit
  PLUSMINS=$((PLUSMINS+1))
done

exit 0
