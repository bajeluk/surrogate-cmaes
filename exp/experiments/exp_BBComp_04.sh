#!/bin/bash
#
# Job-submitting experiment for the first 10D GP model testing
# It is expected to be located in exp/experiments/

export EXPID='exp_BBComp_04'

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

minutes_to_hours() {
  H=$(($1/60)); M=$(($1%60));
  TIME=`printf '%02d%02d' $((HOURS+H)) $((MINS+M))`
}

# critical characters has to be replaced in $OPTS:
# '|' with ',' and "%" with "'"
OPTS=""

# 1-day jobs:

HOURS=18
MINS=20
PLUSMINS=0

QUEUE="24:00:00"
for INST in `seq 1 100 500`; do
  minutes_to_hours $PLUSMINS
  SUBMIT_PBS_PARAMS="-a $TIME"
  submit
  PLUSMINS=$((PLUSMINS+1))
done

# 2-day jobs:
SUBMIT_PBS_PARAMS=""
QUEUE="48:00:00"
for INST in `seq 501 100 700`; do
  submit
done

# 4-day jobs:

MEMORY="4048mb"
QUEUE="96:00:00"
for INST in `seq 701 100 1000`; do
  submit
done

exit 0
