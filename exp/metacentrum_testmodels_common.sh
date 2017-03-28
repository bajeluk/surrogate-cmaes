# exp/metacentrum_testmodels_common.sh
#
# support script for GP model testing experiment submitter
#
# Usage: to be called in experiment submitter via
#
# . ../metacentrum_testmodels_common.sh
#
# see also:
#   exp/experiments/exp_testmodels_example.sh
#   exp/metacentrum_master_template.sh

# CWD = Directory of this particular file
export EXPPATH_SHORT="$CWD"
METACENTRUM_MODEL_BINARY=$EXPPATH_SHORT/../metacentrum_testmodels

# Parse command line options
#
getopt -o nh --long dry-run,help -- "$@"
DRY_RUN=0
HELP=0
while true; do
  case "$1" in
    -h | --help )    HELP=1; shift ;;
    -n | --dry-run ) DRY_RUN=1; echo "---- dry-run ----"; shift ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

if [ "$HELP" = 1 ]; then
  echo "Job-submitting experiment example (submitting the all considered models"
  echo "in 2D)."
  echo "The resulting experiment .sh file is expected to be located in exp/experiments/"
  echo ""
  echo "usage: $ exp/exp_TESTMODELS_EXPERIMENT.sh"
  echo ""
  echo "options:"
  echo "  -n | --dry-run      no submittion will happen, only print what should be done"
  echo ""
  echo "see also:"
  echo "  exp/bash_settings.sh"
  echo "  exp/metacentrum_testmodels_common.sh "
  exit 0
fi

if [ "$useMCR" = 1 ]  &&  ! make -q model; then
  echo "Warning: the binary $METACENTRUM_MODEL_BINARY"
  echo "         is out of date."
  echo "         Running MCR compilation (make model) in 5 sec..."
  sleep 5
  make model
else
  echo "###############################################################"
  echo "Warning: the binary $METACENTRUM_MODEL_BINARY"
  echo "         seems to be up-to-date."
  echo "         No compilation will happen now."
  echo "###############################################################"
  sleep 3
fi

# make sure that these variables will be exported
export ID
export DIM
export FUNC
export INST
export OPTS
export MATLAB_FCN

subtask() {
  if [ -n "$1" ]; then
    JOBNAME_SUFFIX="__$1"
  else
    JOBNAME_SUFFIX=""
  fi

  if [ "$useMCR" = 1 ]; then
    echo "MCR binary submit: ID=$ID : DIM=$DIM : FUNC=$FUNC : INST=$INST : OPTS=$OPTS"
    if [ "$DRY_RUN" = 0 ]; then
      qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,EXPID,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_binary_metajob.sh
    else
      echo qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,EXPID,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_binary_metajob.sh
    fi

  else
    echo "ID=$ID : DIM=$DIM : FUNC=$FUNC : INST=$INST : MATLAB_FCN=$MATLAB_FCN : OPTS=$OPTS"
    if [ "$DRY_RUN" = 0 ]; then
      qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,MATLAB_FCN,EXPID,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_metajob.sh && echo "submitted ok."
    else
      echo qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,MATLAB_FCN,EXPID,EXPPATH_SHORT $EXPPATH_SHORT/../modelTesting_metajob.sh && echo "submitted ok."
    fi
  fi
  ID=$((ID+1))
}

function submit_sequence()
{
  # submit LOW_IDX STEP HIGH_IDX
  #
  # call subtask() with 'modelOptionsIndices' set to a small
  # numbers of model-settings
  #
  # it submits a job for every $2 consecutive modelOptions' indices
  #
  # e.g. submit 11 4 22
  #
  #   submits jobs for the following modelOption indices:
  #
  #   - 11:14
  #   - 15:18
  #   - 19:22
  #
  # Important: the HIGH_IDX should be the equal to (k*LOW_IDX) - 1
  #            for some k
  L=$1
  STEP=$2
  H=$3

  SEQ=`seq $L $STEP $H`
  for i in ${SEQ}; do
    thisL=$i
    thisH=$((i+STEP-1))
    OPTS="struct(%modelOptionsIndices%|${thisL}:${thisH})"
    echo modelOptionsIndices: ${thisL}:${thisH}
    subtask "${thisL}_${thisH}"
  done
}

