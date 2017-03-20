#!/bin/bash
#
# Job-submitting script for GP model testing

# QUEUE = Metacentrum walltime (2h/4h/1d/2d/1w) -- queue will be decided accordingly
export EXPID='exp_GPtest_01'

# Enable this option for using Matlab MCR compilated binaries:
# export useMCR=1

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
. $CWD/../bash_settings.sh
export EXPPATH_SHORT="$CWD"

export ID
export DIM
export FUNC
export INST
export OPTS
export MATLAB_FCN

subtask() {
  # make sure that these variables will be exported
  export EXPPATH_SHORT
  export ID
  export DIM
  export INST
  export OPTS
  export FUNC=`echo $FUNC | tr '\n' ' '`
  if [ "$useMCR" == 1 ]; then
    echo "MCR binary submit: ID=$ID : DIM=$DIM : FUNC=$FUNC : INST=$INST"
    qsub -N "${EXPID}__${ID}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_binary_metajob.sh
  else
    echo ID=$ID : DIM=$DIM : FUNC=$FUNC : INST=$INST : MATLAB_FCN=$MATLAB_FCN
    # submission on Metacentrum
    # TODO: convert into the new PBS-Pro task scheduler
    qsub -N "${EXPID}__$1" -l "walltime=$QUEUE" -v FUNC,DIM,INST,MATLAB_FCN,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_metajob.sh && echo "submitted ok."
  fi
  ID=$((ID+1))
}

OPTS=""

QUEUE=4h
ID=1;

MATLAB_FCN="exp_GPtest_01"
INST="[1 2 3 4 5 41 42 43 44 45 46 47 48 49 50]"

DIM=2
# for i in `seq 1 24`; do
for i in 1 2 3; do
  FUNC=$i; subtask $ID
done

exit 0

DIM=10
for i in `seq 1 24`; do
  FUNC=$i; subtask $ID
done

