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
  echo -e "\nUsage:\n$SCRIPT -ecrqng [EXPID] [ID1] [ID2] ...\n"
  echo "where [ID*] are the experiment ID# (accessible from allids.txt)"
  echo "Command line switches are optional. The following switches are recognized."
  echo "  ${REV}-c${NORM}  --List completed tasks."
  echo "  ${REV}-r${NORM}  --List running tasks."
  echo "  ${REV}-q${NORM}  --List queued tasks."
  echo "  ${REV}-e${NORM}  --Tasks that are not completed, running, exiting, or queued according to qstat, but have stdout or temporary result file."
  echo "  ${REV}-n${NORM}  --List not completed, running, queued, exiting, or errored tasks."
  echo "  ${REV}-g${NORM}  --List tasks in groups according to states."
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
  # TODO: filter according to user
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
# json and results usage
USE_JSON=0
USE_RESULTS=0
# show groups
SHOW_GROUPS=0

#TODO: no options - default : g => show all in groups

while getopts "ecrqngh:" FLAG; do
  case $FLAG in
    c)  #set option "c"
      COMPLETED=1
      USE_RESULTS=1
      ;;
    r)  #set option "r"
      RUNNING=1
      USE_JSON=1
      ;;
    q)  #set option "q"
      QUEUED=1
      USE_JSON=1
      ;;
    n)  #set option "n"
      NOSTATE=1
      USE_JSON=1
      USE_RESULTS=1
      ;;
    e)  #set option "e"
      ERRORED=1
      USE_JSON=1
      USE_RESULTS=1
      ;;
    g) #show groups
      SHOW_GROUPS=1
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

# if no groups of states are chosen and show groups option is on, show all groups
if [ $SHOW_GROUPS -eq 1 ] && [ $COMPLETED -eq 0 ] && [ $RUNNING -eq 0 ] && [ $QUEUED -eq 0 ] && [ $NOSTATE -eq 0 ] && [ $ERRORED -eq 0 ]; then
  COMPLETED=1
  RUNNING=1
  QUEUED=1
  NOSTATE=1
  ERRORED=1
  USE_JSON=1
  USE_RESULTS=1
fi

EXPID="$1"
shift;
IDS="$*"
TMPFILE="/tmp/meta_list_$$"

# CWD = Directory of this particular file
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd $CWD/experiments/$EXPID

# no IDS as input -> load from allids.txt
if [ -z "${IDS}" ]; then
  if [ -f allids.txt ]; then
    IDS=$(cat allids.txt)
  else
    echo "Ids not defined (allids.txt does not exist)"
    exit 1
  fi
fi

# list full and temporary results
if [ $USE_RESULTS -eq 1 ]; then
  # stderr of ls is supressed: 2> /dev/null
  ls -1 ${EXPID}_results_*[0-9].mat ${EXPID}_tmp_*[0-9].mat 1> $TMPFILE 2> /dev/null
fi

# list qstat for running, queued, and error states
QSTAT_OUT=""
if [ $USE_JSON -eq 1 ] && [ -n $QSTAT_OUT ]; then
  # check Kerberos credentials
  if ! klist -s; then
    echo "No Kerberos credentials found. Running kinit"
    kinit
    # finish if still no credentials
    if ! klist -s; then exit 1; fi
  fi

  QSTAT_OUT=`tempfile`.json
  # export qstat to json file
  # TODO: export only USERs tasks
  qstat -f -F json > $QSTAT_OUT
  # remove non-printable characters (Ctrl+A - Ctrl+Z)
  sed -i 's/[\x01-\x1A]//g' $QSTAT_OUT

  # prepare lists
  if [ $RUNNING -eq 1 ] || [ $ERRORED -eq 1 ] || [ $NOSTATE -eq 1 ]; then
    R_LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID R`
  fi
  if [ $QUEUED -eq 1 ] || [ $ERRORED -eq 1 ] || [ $NOSTATE -eq 1 ]; then
    Q_LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID Q`
  fi
  if [ $ERRORED -eq 1 ] || [ $NOSTATE -eq 1 ]; then
    E_LIST=`QSTAT_FILTER $QSTAT_OUT $EXPID E`
  fi
fi

# show results as one output
if [ $SHOW_GROUPS -eq 0 ]; then

  for i in $IDS; do

    # completed
    if [ $COMPLETED -eq 1 ]; then
      RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
      if grep -q "$RES_FILE" $TMPFILE; then
        echo -n "$i "
        continue;
      fi
    fi

    # running
    if [ $RUNNING -eq 1 ]; then
      if echo $R_LIST | grep -q "\"${i}\""; then
        echo -n "$i "
        continue;
      fi
    fi

    # queued
    if [ $QUEUED -eq 1 ]; then
      if echo $Q_LIST | grep -q "\"${i}\""; then
        echo -n "$i "
        continue;
      fi
    fi

    # errored and no state
    if [ $ERRORED -eq 1 ] || [ $NOSTATE -eq 1 ]; then
      RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
      TMP_RES_FILE="${EXPID}_tmp_[0-9D_]*_${i}.mat"
      echo $Q_LIST | grep -q "\"${i}\"" || echo $R_LIST | grep -q "\"${i}\"" || echo $E_LIST | grep -q "\"${i}\""
      if [ $? -eq 0 ] || grep -q "$RES_FILE" $TMPFILE; then
        # job running, queued, exiting, or completed -> skip
        continue;
      fi

      OUTS=`find ${CWD}/.. -maxdepth 1 -name "${EXPID}__${i}.o*"`

      if [ -n "$OUTS" ] ||  grep -q "$TMP_RES_FILE" $TMPFILE; then
        # output or tmp file is present => errored 
        if [ $ERRORED -eq 1 ]; then
          echo -n "${i} ";
        else
          continue;
        fi
      else
        # no file is present => no state
        if [ $NOSTATE -eq 1 ]; then
          echo -n "${i} ";
        else
          continue;
        fi
      fi

    fi

  done
  echo ""

# show results in groups
else

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
    echo "Running:"
    for i in $IDS; do
      if echo $R_LIST | grep -q "\"${i}\""; then
        echo -n "$i "
      fi
    done
    echo ""
  fi

  # queued
  if [ $QUEUED -eq 1 ]; then
    echo "Queued:"
    for i in $IDS; do
      if echo $Q_LIST | grep -q "\"${i}\""; then
        echo -n "$i "
      fi
    done
    echo ""
  fi

  # errored out
  if [ $ERRORED -eq 1 ]; then
    echo "Errored out:"
    for i in $IDS; do
      RES_FILE="${EXPID}_results_[0-9D_]*_${i}.mat"
      TMP_RES_FILE="${EXPID}_tmp_[0-9D_]*_${i}.mat"
      echo $Q_LIST | grep -q "\"${i}\"" || echo $R_LIST | grep -q "\"${i}\"" || echo $E_LIST | grep -q "\"${i}\""
      if [ $? -eq 0 ] || grep -q "$RES_FILE" $TMPFILE; then
        # job running, queued, exiting, or completed -> skip
        continue;
      fi

      OUTS=`find ${CWD}/.. -maxdepth 1 -name "${EXPID}__${i}.o*"`
      if [ -n "$OUTS" ] ||  grep -q "$TMP_RES_FILE" $TMPFILE; then
        # output or tmp file is present => errored 
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
      TMP_RES_FILE="${EXPID}_tmp_[0-9D_]*_${i}.mat"
      echo $Q_LIST | grep -q "\"${i}\"" || echo $R_LIST | grep -q "\"${i}\"" || echo $E_LIST | grep -q "\"${i}\""
      if [ $? -eq 0 ] || grep -q "$RES_FILE" $TMPFILE; then
        # job running, queued, or completed -> skip
        continue;
      fi

      OUTS=`find ${CWD}/.. -maxdepth 1 -name "${EXPID}__${i}.o*"`
      if [ -z "$OUTS" ] || ! grep -q "$TMP_RES_FILE" $TMPFILE; then
        # output nor tmp file is present => no state
        echo -n "${i} ";
      fi
    done
    echo ""
  fi

fi

# remove temporary files
if [ $USE_RESULTS -eq 1 ]; then
  rm $TMPFILE
fi
if [ $USE_JSON -eq 1 ]; then
  rm $QSTAT_OUT
fi

exit 0
