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
  echo -e "\nUsage:\n$SCRIPT -crqn [EXPID] [ID1] [ID2] ...\n"
  echo "where [ID*] are the experiment ID# (accessible from allids.txt)"
  echo "Command line switches are optional. The following switches are recognized."
  echo "  ${REV}-c${NORM}  --List completed tasks."
  echo "  ${REV}-r${NORM}  --List running tasks."
  echo "  ${REV}-q${NORM}  --List queued tasks."
  echo "  ${REV}-n${NORM}  --List not completed, running or queued tasks."
  echo -e "  ${REV}-h${NORM}  --Displays this help message. No further functions are performed."\\n
  echo -e "Example: ${BOLD}$SCRIPT -cr exp_CMA-ES_01 {100..200}${NORM}"\\n
  exit 1
}

if [ "$#" -eq 0 ]; then
  HELP
fi

# state constants
COMPLETED=0
RUNNING=0
QUEUED=0
NOSTATE=0

#TODO: no options - default : crqn

while getopts "crqnh:" FLAG; do
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
    if grep -q "calculating_${i}$" $TMPFILE; then
      echo -n "$i "
    fi
  done
  echo ""
fi

# queued
if [ $QUEUED -eq 1 ]; then
  echo "Queued (removing queued files is not reliable):"
  for i in $IDS; do
    if grep -q "queued_${i}$" $TMPFILE; then
      echo -n "$i "
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

exit 0
