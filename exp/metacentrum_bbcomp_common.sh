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
METACENTRUM_BINARY=$EXPPATH_SHORT/../metacentrum_bbcomp_task
METACENTRUM_BBCOMP_TASK_SHELL=$EXPPATH_SHORT/../metacentrum_bbcomp_task.sh

# This is (1) make target and (2) when not using MCR compilation:
MATLAB_FCN=metacentrum_bbcomp_task

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
  echo "== Job-submitting for BBComp Metacentrum experiments =="
  echo "The resulting experiment .sh file is expected to be located in exp/experiments/"
  echo ""
  echo "usage: $ exp/experiments/[EXPID].sh"
  echo ""
  echo "options:"
  echo "  -n | --dry-run      no submittion will happen, only print what should be done"
  echo ""
  echo "see also:"
  echo "  exp/bash_settings.sh"
  echo "  exp/metacentrum_bbcomp_common.sh "
  exit 0
fi

if [ "$useMCR" = 1 ]  &&  ! make -q $MATLAB_FCN; then
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
export EXPID
export INST
export OPTS
export MATLAB_FCN

submit() {
  if [ -n "$1" ]; then
    JOBNAME_SUFFIX="__$1"
  else
    JOBNAME_SUFFIX=""
  fi

  if [ "$useMCR" = 1 ]; then
    echo "MCR binary submit: INST=$INST : OPTS=$OPTS"
    if [ "$DRY_RUN" = 0 ]; then
      qsub -N "${EXPID}__${INST}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v INST,OPTS,EXPID,EXPPATH_SHORT "$METACENTRUM_BBCOMP_TASK_SHELL" && echo "submitted ok."
    else
      echo qsub -N "${EXPID}__${INST}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,OPTS,EXPID,EXPPATH_SHORT "$METACENTRUM_BBCOMP_TASK_SHELL"
    fi

  else
    echo "INST=$INST : MATLAB_FCN=$MATLAB_FCN : OPTS=$OPTS"
    if [ "$DRY_RUN" = 0 ]; then
      qsub -N "${EXPID}__${INST}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v INST,OPTS,MATLAB_FCN,EXPID,EXPPATH_SHORT "$METACENTRUM_BBCOMP_TASK_SHELL" && echo "submitted ok."
    else
      echo qsub -N "${EXPID}__${INST}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v INST,OPTS,MATLAB_FCN,EXPID,EXPPATH_SHORT "$METACENTRUM_BBCOMP_TASK_SHELL"
    fi
  fi
  ID=$((ID+1))
}
