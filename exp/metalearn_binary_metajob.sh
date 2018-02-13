#!/bin/sh
#PBS -l select=1:ncpus=1:mem=2500mb:scratch_local=1gb:cgroups=cpuacct

# cgroups option limits resources usage on PBS nodes, see:
# https://wiki.metacentrum.cz/wiki/Cgroupy
#
# a bit restrictive: PBS -l select=1:ncpus=1:mem=1500mb:scratch_local=1gb:cl_minos=False:cl_mudrc=False:cl_mandos=False:cl_losgar=False:cl_haldir=False
# very restrictive: PBS -l select=1:ncpus=1:mem=1500mb:scratch_local=1gb:cl_minos=False:cl_mudrc=False:cl_krux=False:cl_perian=False:cl_phi=False:cl_zewura=False:cl_zebra=False:cl_mandos=False:cl_loslab=False:cl_losgar=False:cl_haldir=False

# it supposes the following variables are set:
#
#   FUNC           -- list of integers of functions
#   DIM            -- list of integers of dimensions
#   INST           -- list of instances to process
#   DESIGN         -- design types to test
#   MODEL          -- list of model types to test
#   OPTS           -- string with options to be eval()-ed
#   EXPID          -- string with the experiment name
#   FNAME_TEMPLATE -- template for data set file name
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments

# MATLAB Runtime environment
# setting LD_LIBRARY_PATH moved into bash_settings.sh

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

MATLAB_BINARY_CALL="exp/metacentrum_metalearn"
if [ -z "$DATASET_PATH" ]; then
  DATASET_PATH="$EXPPATH_SHORT/data_metalearn"
  echo "There was no dataset path, setting the dataset path to default: $DATASET_PATH"
elif ! grep -q '/.*$' <<< "$DATASET_PATH"; then
  DATASET_PATH="$EXPPATH_SHORT/$DATASET_PATH"
  echo "There was no path to dataset directory, setting to: $DATASET_PATH"
fi

if [ ! -d DATASET_PATH ]; then
  echo "Error: Data set path \"$DATASET_PATH\" not found."; exit 1
fi

export SCRATCHDIR
export LOGNAME

if [ -z "$FUNC" ] ; then
  echo "Error: Task FUNC numbers are not known"; exit 1
fi
if [ -z "$DIM" ] ; then
  echo "Error: Task DIM numbers are not known"; exit 1
fi
if [ -z "$INST" ]; then
  echo "Warning: INST is empty, default instances will be used."
  INST="[1:5 40:50]";
fi
if [ -z "$DESIGN" ]; then
  echo "Warning: DESIGN empty, default 'lhs' design will be used."
  DESIGN="{'lhs'}"
fi
if [ -z "$MODEL" ]; then
  echo "Error: Model types not known"; exit 1
fi
if [ -z "$FNAME_TEMPLATE" ]; then
  echo "Warning: FNAME_TEMPLATE empty, using a default."
  FNAME_TEMPLATE='data_metalearn_%dD_f%d_inst%d_N%d_design-%s'
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

# replace critical characters in $OPTS: '|' with ',' and ':' with ';'
# and "%" with "'"
OPTS=`echo $OPTS | tr '%|:' "',;"`

####### PREPARE DATA #####
#
cd "$EXPPATH_SHORT/../.."
DST_DIR=$SCRATCHDIR/`basename "$DATASET_PATH"`
mkdir -p $DST_DIR

for dim in $( eval_matlab_array "$DIM" ); do
  mkdir -p $DST_DIR/${DIM}D || :

  for func in $( eval_matlab_array "$FUNC" ); do
    for inst in $( eval_matlab_array "$INST" ); do
      for design in $( eval_matlab_cell "$DESIGN" ); do
        FNAME=$( fprintf $FNAME_TEMPLATE $dim $func $inst $design )
        SRC_FILE="$DATASET_PATH/${DIM}D/$FNAME"

        if [ ! -f $SRC_FILE ]; then
          echo "Error: Data set \"$SRC_FILE\" not found, exiting."; exit 1
        fi
        echo cp $SRC_FILE $DST_DIR/${DIM}D
        cp $SRC_FILE $DST_DIR/${DIM}D
      done
    done
  done
done
DATASET_PATH=$DST_DIR
#
########################

echo "====================="
echo -n "Current dir:    "; pwd
echo -n "Current node:   "; cat "$PBS_NODEFILE"
echo    '$HOME =         ' $HOME
echo    '$MCR_CACHE_ROOT = ' $MCR_CACHE_ROOT
echo    '$DATASET_PATH =   ' $DATASET_PATH
echo "====================="

######### CALL #########
#
echo '##############'
echo Will be called: $MATLAB_BINARY_CALL \"$EXPID\" \"$EXPPATH_SHORT\" \"$FUNC\" \"$DIM\" \"$INST\" \"$OPTS\" \"$DATASET_PATH\"
echo '##############'

$MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" "$FUNC" "$DIM" "$INST" "$OPTS" "$DATASET_PATH"
#
########################

if [ $? -eq 0 ]; then
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== FINISHED ===="
  rm -rf $SCRATCHDIR/*
  exit 0
else
  echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID**  ==== ENDED WITH ERROR! ===="
  exit 1
fi

function eval_matlab_array {
  local S=$1
  local LO=""
  local HI=""
  local EXPANDED=""

  # expand lo:hi expressions
  local COL_EXPR=$( echo "$S" | egrep --color=never -o "[0-9]+:[0-9]+" )

  if [ ! -z $COL_EXPR ]; then
    while read E; do
      LO=$( echo $E | cut -d':' -f1 )
      HI=$( echo $E | cut -d':' -f2 )
      EXPANDED=$( eval echo {$LO..$HI} )
      S=$( echo $S | sed -e s"/$E/$EXPANDED/g" )
    done <<< $( echo $COL_EXPR | tr ' ' '\n' )
  fi

  S=$( echo $S | tr "[]," " " )
  echo $S
}

function eval_matlab_cell {
  echo $1 | tr "{},'" " "
}