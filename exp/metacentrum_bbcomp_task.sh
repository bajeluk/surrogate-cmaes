#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2000mb:scratch_local=5gb

# it suppose the following variables set:
#
#   INST           -- list of instances to process
#   OPTS           -- string with options to be eval()-ed
#   EXPID          -- string with the experiment name
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments

# MATLAB Runtime environment
# setting LD_LIBRARY_PATH moved into bash_settings.sh

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

MATLAB_BINARY_CALL="exp/metacentrum_bbcomp_task"
BBCOMP_PROXY_STARTPORT=`sed -n 's/bbcompParams.proxyPort \?= \?\([0-9]\+\)[^0-9]*/\1/p' $EXPPATH_SHORT/${EXPID}.m`

run_proxy() {
  PROXY_PID=`ps x | sed -n "s/\(^[0-9]\+\).*proxy ${1}$/\1/p"`
  if [ -n "$PROXY_PID" ]; then
    echo "Proxy is already running with PID=$PROXY_PID"
    return $PROXY_PID
  fi

  this_dir=`pwd`
  cd $EXPPATH_SHORT/../vendor/bbcomp/proxy/;
  LD_LIBRARY_PATH=".:$LD_LIBRARY_PATH" ./proxy $1 &
  export PROXY_PID=$!
  cd "$this_dir"
  echo "Proxy PID=$PROXY_PID"
  sleep 2
  if [ -d /proc/$PROXY_PID ]; then
    return $PROXY_PID
  else
    return 0
  fi
}

export SCRATCHDIR
export LOGNAME

if [ -z "$EXPID" ] ; then
  echo "Error: EXPID is empty!"; exit 1
fi
if [ -z "$INST" ]; then
  echo "Warning: INST is empty, default instances will be used."
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

# replace critical characters in $OPTS: '|' with ',' and "%" with "'"
OPTS=`echo $OPTS | tr '%|' "',"`

cd "$EXPPATH_SHORT/../.."

echo "====================="
echo -n "Current dir:    "; pwd
echo -n "Current node:   "; cat "$PBS_NODEFILE"
echo    '$HOME =         ' $HOME
echo    '$MCR_CACHE_ROOT = ' $MCR_CACHE_ROOT
echo    '$DATASET =      ' $DATASET
echo "====================="

######### CALL #########
#
PROXY_PORT=$((BBCOMP_PROXY_STARTPORT+INST))
run_proxy $PROXY_PORT
if [ "$PROXY_PID" -ne 0 ]; then
  echo "Proxy should be running at port ${PROXY_PORT} with PID ${PROXY_PID}..."
else
  echo "Proxy is not successfully started!"
  exit 0;
fi

echo ''
echo '##############'
echo Will be called: $MATLAB_BINARY_CALL \"$EXPID\" \"$EXPPATH_SHORT\" \"$INST\" \"$OPTS\"
echo '##############'

$MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" "$INST" "$OPTS"
#
########################

if [ -d "/proc/$PROXY_PID" ]; then
  echo "Killing BBComp proxy with PID $PROXY_PID ..."
  kill $PROXY_PID
fi

if [ $? -eq 0 ]; then
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== FINISHED ===="
  if [ -n "$SCRATCHDIR" ]; then
    rm -rf $SCRATCHDIR/*
  fi
  exit 0
else
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== ENDED WITH ERROR! ===="
  exit 1
fi
