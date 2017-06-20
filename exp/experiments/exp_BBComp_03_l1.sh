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
OPTS=""

QUEUE="24:00:00"
for INST in `seq 11 20` `seq 111 120` `seq 211 220`; do
  submit
done

QUEUE="48:00:00"
for INST in `seq 311 320` `seq 411 420`; do
  submit
done

QUEUE="168:00:00"
for INST in `seq 511 520` `seq 611 620` `seq 711 720` `seq 811 820` `seq 911 920`; do
  submit
done

exit 0
