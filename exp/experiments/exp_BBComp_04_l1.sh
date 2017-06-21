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
  H=$(((MINS+$1)/60)); M=$(((MINS+$1)%60));
  TIME=`printf '%02d%02d' $(((HOURS+H)%24)) $M`
}

# critical characters has to be replaced in $OPTS:
# '|' with ',' and "%" with "'"
OPTS=""

# 1-day jobs:

HOURS=20
MINS=00
PLUSMINS=0

QUEUE="24:00:00"
for INST in `seq 9 50` `seq 102 150` `seq 202 250` `seq 302 350` `seq 402 450`; do
  minutes_to_hours $PLUSMINS
  SUBMIT_PBS_PARAMS="-a $TIME"
  submit
  PLUSMINS=$((PLUSMINS+1))
done

# 2-day jobs:

HOURS=04
MINS=00
PLUSMINS=0

QUEUE="48:00:00"
for INST in `seq 502 550` `seq 602 650` `seq 702 750`; do
  minutes_to_hours $PLUSMINS
  SUBMIT_PBS_PARAMS="-a $TIME"
  submit
  PLUSMINS=$((PLUSMINS+1))
done

# 4-day jobs:

SUBMIT_PBS_PARAMS=""
MEMORY="4048mb"
QUEUE="96:00:00"
for INST in `seq 802 850` `seq 902 950`; do
  submit
done

exit 0
