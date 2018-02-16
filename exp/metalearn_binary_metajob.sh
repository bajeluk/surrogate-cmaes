#!/bin/bash
#PBS -l select=1:ncpus=10:mem=500mb:scratch_local=1gb:cgroups=cpuacct

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
#   DATASIZE       -- list of data sizes to test
#   OPTS           -- string with options to be eval()-ed
#   EXPID          -- string with the experiment name
#   FNAME_TEMPLATE -- template for data set file name
#   EXPPATH_SHORT  -- usually $APPROOT/exp/experiments


function eval_matlab_array() {
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

function eval_matlab_cell() {
  echo $1 | tr "{},'" " "
}


# MATLAB Runtime environment
# setting LD_LIBRARY_PATH moved into bash_settings.sh

# Load global settings and variables
. $EXPPATH_SHORT/../bash_settings.sh

MATLAB_BINARY_CALL="exp/metacentrum_metalearn"
if [ -z "$DATASET_PATH" ]; then
  DATASET_PATH="$EXPPATH_SHORT/data_metalearn"
  echo "Warning: There was no dataset path, setting the dataset path to default: $DATASET_PATH"
elif ! grep -q '/.*$' <<< "$DATASET_PATH"; then
  DATASET_PATH="$EXPPATH_SHORT/$DATASET_PATH"
  echo "Warning: There was no path to dataset directory, setting to: $DATASET_PATH"
fi

if [ ! -d $DATASET_PATH ]; then
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
if [ -z "$DATASIZE" ]; then
  echo "Error: Data sizes not known"; exit 1
fi
if [ -z "$FNAME_TEMPLATE" ]; then
  echo "Warning: FNAME_TEMPLATE empty, using a default."
  FNAME_TEMPLATE='data_metalearn_%dD_f%d_inst%d_N%d_design-%s.mat'
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
fi

# replace critical characters: '|' with ',' and '@' with ';'
# and "%" with "'"
OPTS=`echo $OPTS | tr '%|@' "',;"`
DATASIZE=`echo $DATASIZE | tr '%|@' "',;"`
DESIGN=`echo $DESIGN | tr '%|@' "',;"`

# debug:
# DATASIZE=`echo $DATASIZE | sed s"/%/''/g"`
# DESIGN=`echo $DESIGN | sed s"/%/''/g"`
# OPTS=`echo $OPTS | sed s"/%/''/g" | sed s"/|/,/g"`

cd "$EXPPATH_SHORT/../.."

####### PREPARE DATA #####
#
DST_DIR="$SCRATCHDIR/`basename $DATASET_PATH`"

for dim in $( eval_matlab_array "$DIM" ); do
  mkdir -p "$DST_DIR/${dim}D" || exit 1

  for func in $( eval_matlab_array "$FUNC" ); do
    for inst in $( eval_matlab_array "$INST" ); do
      for N_expr in $( eval_matlab_cell "$DATASIZE" ); do
        N=$(( N_expr ))
        for design in $( eval_matlab_cell "$DESIGN" ); do
          FNAME=$( printf $FNAME_TEMPLATE $dim $func $inst $N $design )
          SRC_FILE="$DATASET_PATH/${dim}D/$FNAME"

          if [ ! -f $SRC_FILE ]; then
            echo "Error: Data set \"$SRC_FILE\" not found, exiting."; exit 1
          fi
          echo cp $SRC_FILE $DST_DIR/${dim}D
          cp $SRC_FILE $DST_DIR/${dim}D
        done
      done
    done
  done
done
#
########################

DATASET_PATH=$DST_DIR

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
echo "Will be called: $MATLAB_BINARY_CALL \"$EXPID\" \"$EXPPATH_SHORT\" \"$DIM\" \"$FUNC\" \"$INST\" \"$DATASIZE\" \"$DESIGN\" \"$OPTS\" \"$DATASET_PATH\""
echo '##############'

$MATLAB_BINARY_CALL "$EXPID" "$EXPPATH_SHORT" "$DIM" "$FUNC" "$INST" "$DATASIZE" "$DESIGN" "$OPTS" "$DATASET_PATH"
# debug:
# matlab -nodisplay -nodesktop -r "dbstop in metacentrum_metalearn; metacentrum_metalearn('$EXPID', '$EXPPATH_SHORT', '$DIM', '$FUNC', '$INST', '$DATASIZE', '$DESIGN', '$OPTS', '$DATASET_PATH')"
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
