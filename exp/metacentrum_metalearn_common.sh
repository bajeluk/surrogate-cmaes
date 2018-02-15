# exp/metacentrum_metalearn_common.sh
#
# support script for metalearning testing experiment submitter
#
# Usage: to be called in experiment submitter via
#
# . ../metacentrum_metalearn_common.sh
#
# see also:
#   exp/experiments/exp_testmodels_common.sh
#   exp/experiments/exp_testmodels_example.sh
#   exp/metacentrum_master_template.sh

# CWD = Directory of this particular file
export EXPPATH_SHORT="$CWD"
METACENTRUM_METALEARN_BINARY=$EXPPATH_SHORT/../metacentrum_metalearn

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
  echo "usage: $ exp/exp_METALEARN_EXPERIMENT.sh"
  echo ""
  echo "options:"
  echo "  -n | --dry-run      no submittion will happen, only print what should be done"
  echo ""
  echo "see also:"
  echo "  exp/bash_settings.sh"
  echo "  exp/metacentrum_testmodels_common.sh "
  exit 0
fi

if [ "$useMCR" = 1 ]  &&  ! make -q metalearn; then
  echo "Warning: the binary $METACENTRUM_METALEARN_BINARY"
  echo "         is out of date."
  echo "         Running MCR compilation (make metalearn) in 5 sec..."
  #sleep 5
  make metalearn
else
  echo "###############################################################"
  echo "Warning: the binary $METACENTRUM_METALEARN_BINARY"
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
export MODEL
export DESIGN
export DATASIZE
export OPTS
export MATLAB_FCN

subtask() {
  if [ -n "$1" ]; then
    JOBNAME_SUFFIX="__$1"
  else
    JOBNAME_SUFFIX=""
  fi

  echo "MCR binary submit: ID=$ID : DIM=$DIM : FUNC=$FUNC : INST=$INST : MODEL=$MODEL : DESIGN=$DESIGN : DATASIZE=$DATASIZE : OPTS=$OPTS : DATASET=`basename $DATASET_PATH`"
  if [ "$DRY_RUN" = 0 ]; then
    qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,MODEL,DESIGN,DATASIZE,OPTS,EXPID,EXPPATH_SHORT,DATASET $EXPPATH_SHORT/../metalearn_binary_metajob.sh
  else
    echo qsub -N "${EXPID}__${ID}${JOBNAME_SUFFIX}" -l "walltime=$QUEUE" -v FUNC,DIM,INST,MODEL,DESIGN,DATASIZE,OPTS,EXPID,EXPPATH_SHORT,DATASET $EXPPATH_SHORT/../metalearn_binary_metajob.sh
  fi
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
    OPTS="struct(%modelOptionsIndices%|{${thisL}:${thisH}})"
    echo modelOptionsIndices: {${thisL}:${thisH}}
    subtask "${thisL}_${thisH}"
  done
}

function submit_model() {
  # submit_model MODEL DIMS FUNCS INSTS DESIGNS DATASIZES
  #
  # read option ids for given model type
  #
  # submits all model options for each data set configuration
  #
  # data set configurations are given by cartesian product
  # of following sequences:
  #   * DIMS      -- dimensionalities
  #   * FUNCS     -- BBOB functions
  #   * INSTS     -- BBOB instances
  #   * DESIGNS   -- sampling designs
  #   * DATASIZES -- sample sizes
  MODEL=`echo $1 | tr '[:upper:]' '[:lower:]'`
  DIMS=$2
  FUNCS=$3
  INSTS=$4
  DESIGNS=$5
  DATASIZES=$6

  FNAME=$( echo $EXPPATH_SHORT/$EXPID/${MODEL}_ids.txt )

  if [ ! -f $FNAME ]; then
    echo "Error: model option indices file $FNAME not found."
    exit 1
  fi

  read OPT_IND < $FNAME
  OPT_IND=( $OPT_IND ) # bash array
  N_OPTS=${#OPT_IND[@]}

  for DIM in $DIMS; do
    for FUNC in $FUNCS; do
      for INST in $INSTS; do
        for DESIGN in $DESIGNS; do
          for DATASIZE in $DATASIZES; do
            LO=${OPT_IND[0]}
            HI=${OPT_IND[(($N_OPTS - 1))]}
            submit_sequence $LO 1 $HI
          done # data sizes
        done # designs
      done # instances
    done # functions
  done # dims
}