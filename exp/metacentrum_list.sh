#!/bin/bash

# metacentrum_list -crqn [EXPID] [ID1] [ID2]
#

# Help function
function HELP {
  # script name
  SCRIPT=`basename ${BASH_SOURCE[0]}`
  #Set fonts for Help.
  NORM=`tput sgr0`
  BOLD=`tput bold`
  REV=`tput smso`
  echo -e "\nUsage:\n$SCRIPT -ecrqn [EXPID] [ID1] [ID2] ...\n"
  echo "where [ID*] are the experiment ID# (accessible from allids.txt)"
  echo "Command line switches are optional. The following switches are recognized."
  echo "  ${REV}-c${NORM}  --List completed tasks."
  echo "  ${REV}-r${NORM}  --List running tasks."
  echo "  ${REV}-q${NORM}  --List queued tasks."
  echo "  ${REV}-n${NORM}  --List not completed, running or queued tasks."
  echo "  ${REV}-e${NORM}  --Tasks that are not running or queued according to qstat, but have no FINISHED message in their stdout."
  echo -e "  ${REV}-h${NORM}  --Displays this help message. No further functions are performed."\\n
  echo -e "Example: ${BOLD}$SCRIPT -cr exp_CMA-ES_01 {100..200}${NORM}"\\n
  exit 1
}

# Filter output of qstat in json format using jq parser
function QSTAT_FILTER {
  # jq (json parser) executable
  local JQ_PATH=$CWD/vendor/jq-linux64
  # path to qstat output
  local QSTAT_FILE=$1
  local EXPID=$2
  local JOB_STATE=$3

  # filter jobs by names matching a regexp and state codes equal to a string; output numerical job suffix
  $JQ_PATH ".Jobs[] | select((.Job_Name | tostring | test(\"^${EXPID}__[0-9]+$\")) and .job_state == \"${JOB_STATE}\") | .Job_Name | sub(\"${EXPID}__\"; \"\")" $QSTAT_FILE;
}

if [ "$#" -eq 0 ]; then
  HELP
fi

# state constants
COMPLETED=0
RUNNING=0
QUEUED=0
NOSTATE=0
ERRORED=0

#TODO: no options - default : crqn

while getopts "ecrqnh:" FLAG; do
  case $FLAG in
    c)  #set option "c"
      COMPLETED=1
      ;;
    r)  #set option "r"
      RUNNING=1
      ;;
    q)  #set option "q"
      QUEUED=1
      ;;
    n)  #set option "n"
      NOSTATE=1
      ;;
    e)  #set option "e"
      ERRORED=1
      ;;
    h)  #show help
      HELP
      ;;
    \?) #unrecognized option - show help
      echo -e \\n"Option -${BOLD}$OPTARG${NORM} not allowed."
      HELP
      ;;
  esac
done
shift $((OPTIND-1))  #This tells getopts to move on to the next argument.

EXPID="$1"
shift;
IDS="$*"
TMPFILE="/tmp/meta_list_$$"

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd $CWD/experiments/$EXPID

ls -1 ${EXPID}_results_*[0-9].mat calculating_* queued_* > $TMPFILE

QSTAT_OUT=""

# completed
if [ $COMPLETED -eq 1 ]; then
  echo "Completed:"
  for i in $IDS; do
    RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
    if grep -q "$RES_FILE" $TMPFILE; then
      echo -n "$i "
    fi
  done
  echo ""
fi

# running
if [ $RUNNING -eq 1 ]; then
  if [ -n $QSTAT_OUT ]; then
    QSTAT_OUT=`tempfile`.json
    qstat -f -F json > $QSTAT_OUT
  fi

  echo "Running:"
  LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID R`
  for i in $IDS; do
    if echo $LIST | grep -q "\"${i}\""; then
      echo -n "$i "
    fi
  done
  echo ""
fi

# queued
if [ $QUEUED -eq 1 ]; then
  if [ -n $QSTAT_OUT ]; then
    QSTAT_OUT=`tempfile`.json
    qstat -f -F json > $QSTAT_OUT
  fi

  echo "Queued (removing queued files is not reliable):"
  LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID Q`
  for i in $IDS; do
    if echo $LIST | grep -q "\"${i}\""; then
      echo -n "$i "
    fi
  done
  echo ""
fi

# errored out
if [ $ERRORED -eq 1 ]; then
  echo $CWD/..
  if [ -n $QSTAT_OUT ]; then
    QSTAT_OUT=`tempfile`.json
    qstat -f -F json > $QSTAT_OUT
  fi

  echo "Errored out:"
  Q_LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID Q`
  R_LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID R`
  for i in $IDS; do
    echo $Q_LIST | grep -q "\"${i}\"" || echo $R_LIST | grep -q "\"${i}\""
    if [ $? -eq 0 ]; then
      # job running or queued, skip
      continue;
    fi

    OUTS=`find ${CWD}/.. -maxdepth 1 -name "${EXPID}__${i}.o*"`
    if [ -z "$OUTS" ]; then
      # outputs not present
      continue;
    fi

    grep -q FINISHED $OUTS
    if [ $? -ne 0 ]; then
      echo -n "${i} ";
    fi
  done
  echo ""
fi



# no state
if [ $NOSTATE -eq 1 ]; then
  echo "No status:"
  for i in $IDS; do
    RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
    if ! grep -q "$RES_FILE"'$\|'"calculating_${i}$"'$\|'"queued_${i}$"'$' $TMPFILE; then
      echo -n "$i "
    fi
  done
  echo ""
fi

#rm $TMPFILE
rm /tmp/file*.json

exit 0
