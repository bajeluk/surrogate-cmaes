#!/bin/bash

# split_bbob_results [SOURCE_DIR] [DEST_DIR]
# 

MYALGORITHMS=("GP1-CMAES" "GP5-CMAES" "RF1-CMAES" "RF5-CMAES")
REFALG="CMA-ES"
FUNCTIONS=`seq 1 24`
DIMENSIONS="2 5 10 20"
EXPID="exp_geneEC_08"

if [ -n "$1" ]; then 
  BBOB_RESULTS_DIR=$1
else 
  BBOB_RESULTS_DIR="."
fi
if [ -n "$2" ]; then
  OUTPUT_DIR=$2
else
  OUTPUT_DIR="output"
fi

mkdir -p $BBOB_RESULTS_DIR
mkdir -p $OUTPUT_DIR

echo in :$BBOB_RESULTS_DIR:
echo out:$OUTPUT_DIR:

make_dir_for_each_alg()
{
  NEWDIR=$1
  for dir in ${MYALGORITHMS[*]}; do
    mkdir -p $OUTPUT_DIR/$dir${NEWDIR}
    echo mkdir -p $OUTPUT_DIR/$dir${NEWDIR}
  done
  mkdir -p $OUTPUT_DIR/$REFALG${NEWDIR}
  echo mkdir -p $OUTPUT_DIR/$REFALG${NEWDIR}
}

make_dir_for_each_alg ""

for fun in $FUNCTIONS; do

  make_dir_for_each_alg "/data_f"${fun}

  for dim in $DIMENSIONS; do

    # there has to be at least one such dir:
    if [ -z "`ls -d $BBOB_RESULTS_DIR/${fun}_${dim}D_* 2>/dev/null`" ]; then
      continue
    fi

    # for each directory of the format FUN_DIM_ID/
    for dir in $BBOB_RESULTS_DIR/${fun}_${dim}D_*; do
      ID=`echo ${dir} | sed 's/.*_\([0-9]\+\)$/\1/'`
      ALGNUM=$(( (ID - 1) % 4 ))         # modulo 4
      ALGDIR=${MYALGORITHMS[ALGNUM]}

      echo " fun=${fun}, dim=${dim} (id=${ID}) --> $OUTPUT_DIR/$ALGDIR/"

      # copy the detailed files
      cp ${dir}/data_f${fun}/bbobexp_f* $OUTPUT_DIR/$ALGDIR/data_f${fun}/
      # take the comprehensive info and add it to bbobexp_f#.info
      sed -n "s/'$EXPID[^']*'/'${EXPID}_${ALGDIR}'/;1,3p" ${dir}/bbobexp_f${fun}.info >> $OUTPUT_DIR/$ALGDIR/bbobexp_f${fun}.info
      # newline :)
      echo "" >> $OUTPUT_DIR/$ALGDIR/bbobexp_f${fun}.info

      # test whether cmaes records exist
      if [ -f "${dir}/data_f${fun}/bbobexp-01_f${fun}_DIM${dim}.dat" -a ! -f "$OUTPUT_DIR/$REFALG/data_f${fun}/bbobexp-01_f${fun}_DIM${dim}.dat" ]; then

        echo " cmaes fun=${fun}, dim=${dim} (id = ${ID}) --> $OUTPUT_DIR/$REFALG/data_f${fun}"

        # copy the detailed results: data_f#/bbobexp-01_* files
        cp ${dir}/data_f${fun}/bbobexp-01_f* $OUTPUT_DIR/$REFALG/data_f${fun}/
        # take the comprehensive info and add it to bbobexp_f#.info
        sed -n "s/'$EXPID[^']*'/'${REFALG}'/;4,6p" ${dir}/bbobexp_f${fun}.info >> $OUTPUT_DIR/$REFALG/bbobexp_f${fun}.info
        # newline :)
        echo "" >> $OUTPUT_DIR/$REFALG/bbobexp_f${fun}.info
      fi
    done
  done
done

