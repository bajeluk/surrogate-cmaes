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
MEMORY=""

QUEUE="24:00:00"
for INST in `seq 1 10` `seq 101 110` `seq 201 210`; do
  submit
done

QUEUE="48:00:00"
for INST in `seq 301 310` `seq 401 410`; do
  submit
done

QUEUE="168:00:00"
for INST in `seq 501 510` `seq 601 610` `seq 701 710`; do
  submit
done

MEMORY="4048mb"
QUEUE="168:00:00"
for INST in `seq 801 810` `seq 901 910`; do
  submit
done

exit 0
