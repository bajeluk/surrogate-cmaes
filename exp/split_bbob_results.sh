#!/bin/bash

# split_bbob_results -i [BBOB_RAW_DIR] -o [BBOB_DEST_DIR] -r [REFALG] -e [EXPID] [ALG1] [ALG2]...
# 

#Initialize cmdopts default values
EXPID="exp_SCMAES_01"   # -e
BBOB_RAW_DIR="."        # -i
OUTPUT_DIR='$CWD/../log/bbob_data/$EXPID'  # -o
REFALG="CMA-ES"         # -r
MYALGORITHMS=("ALG1")   # --

FUNCTIONS=`seq 1 24`
DIMENSIONS="2 3 5 10 20"
CWD=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

#Set fonts for Help.
NORM=`tput sgr0`
BOLD=`tput bold`
REV=`tput smso`
SCRIPT=`basename ${BASH_SOURCE[0]}`

#Help function
function HELP {
  echo -e "\nUsage:\nsplit_bbob_results -i [BBOB_RAW_DIR] -o [BBOB_DEST_DIR] -r [REFALG] -e [EXPID] [ALG1] [ALG2]...\n"
  echo "Command line switches are optional. The following switches are recognized."
  echo "  ${REV}-i${NORM}  --Directory with raw BBOB input data. Default is '${BOLD}${BBOB_RAW_DIR}${NORM}'."
  echo "  ${REV}-o${NORM}  --Directory where to put BBOB parsed data. Default is '${BOLD}${OUTPUT_DIR}${NORM}' (will be eval-ed)."
  echo "  ${REV}-r${NORM}  --Referential algorithm name. Default is '${BOLD}${REFALG}${NORM}'."
  echo "  ${REV}-e${NORM}  --EXPID of the current experiment. Default is '${BOLD}${EXPID}${NORM}'."
  echo -e "  ${REV}-h${NORM}  --Displays this help message. No further functions are performed."\\n
  echo -e "Example: ${BOLD}$SCRIPT -i bbob_raw/ -o bbob_data/ -r CMA-ES -e exp_CMA-ES_01${NORM}"\\n
  exit 1
}

if [ "$#" -eq 0 ]; then
  HELP
fi

while getopts "i:o:r:e:h" FLAG; do
  case $FLAG in
    i)  #set option "i"
      BBOB_RAW_DIR=$OPTARG
      echo "-i used: $OPTARG"
      ;;
    o)  #set option "o"
      OUTPUT_DIR=$OPTARG
      echo "-o used: $OPTARG"
      ;;
    r)  #set option "r"
      REFALG=$OPTARG
      echo "-r used: $OPTARG"
      ;;
    e)  #set option "e"
      EXPID=$OPTARG
      echo "-e used: $OPTARG"
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

OUTPUT_DIR=`eval echo $OUTPUT_DIR`
echo "Output dir: '$OUTPUT_DIR'"

i=0
while [ $# -ne 0 ]; do
  MYALGORITHMS[i]=$1
  echo "ALG${i} == '${MYALGORITHMS[i]}'"
  shift
  i=`expr $i + 1`
done

if [ ! -d "$BBOB_RAW_DIR" ]; then
  echo "Input directory '${BBOB_RAW_DIR}' must exist."
  exit 1
fi

mkdir -p $OUTPUT_DIR

echo in :$BBOB_RAW_DIR:
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

merge_resumed_runs()
{
  curdir=`pwd`
  tmpfile=/tmp/merge_info_after_resume_$$.info
  cd "$BBOB_RAW_DIR"
  for d in *D_*; do
    cd $d
    infofile=bbob*.info
    gawk -f "$CWD/merge_info_after_resume.awk" $infofile > $tmpfile
    mv $tmpfile $infofile
    cd ..
  done
  rm $tmpfile
  cd "$curdir"
}

merge_resumed_runs

make_dir_for_each_alg ""

for fun in $FUNCTIONS; do

  make_dir_for_each_alg "/data_f"${fun}

  for dim in $DIMENSIONS; do

    # there has to be at least one such dir:
    if [ -z "`ls -d $BBOB_RAW_DIR/${fun}_${dim}D_* 2>/dev/null`" ]; then
      continue
    fi

    # for each directory of the format FUN_DIM_ID/
    for dir in $BBOB_RAW_DIR/${fun}_${dim}D_*; do
      ID=`echo ${dir} | sed 's/.*_\([0-9]\+\)$/\1/'`
      ALGNUM=$(( (ID - 1) % ${#MYALGORITHMS[*]} ))         # modulo 4
      ALGDIR=${MYALGORITHMS[ALGNUM]}

      echo " fun=${fun}, dim=${dim} (id=${ID}) --> $OUTPUT_DIR/$ALGDIR/"

      # copy the detailed files
      cp ${dir}/data_f${fun}/bbobexp_f* $OUTPUT_DIR/$ALGDIR/data_f${fun}/
      # take the comprehensive info and add it to bbobexp_f#.info
      sed -n "s/'$EXPID[^']*'/'${EXPID}_${ALGDIR}'/;1,3p" ${dir}/bbobexp_f${fun}.info >> $OUTPUT_DIR/$ALGDIR/bbobexp_f${fun}.info
      # newline :)
      # echo "" >> $OUTPUT_DIR/$ALGDIR/bbobexp_f${fun}.info

      # test whether cmaes records exist
      if [ -f "${dir}/data_f${fun}/bbobexp-01_f${fun}_DIM${dim}.dat" -a ! -f "$OUTPUT_DIR/$REFALG/data_f${fun}/bbobexp-01_f${fun}_DIM${dim}.dat" ]; then

        echo " cmaes fun=${fun}, dim=${dim} (id = ${ID}) --> $OUTPUT_DIR/$REFALG/data_f${fun}"

        # copy the detailed results: data_f#/bbobexp-01_* files
        # cp ${dir}/data_f${fun}/bbobexp-01_f* $OUTPUT_DIR/$REFALG/data_f${fun}/

        # filter the detailed results of the additional CMA-ES runs
        # (only the first 15 records are stored): data_f#/bbobexp-01_* files
        $( cd "${dir}/data_f${fun}";
          for f in bbobexp-01_f*; do
            gawk -e 'BEGIN { n = 0 }; /^%/ { n = n + 1 }; { if (n <= 15) { print $0 } }' $f > $OUTPUT_DIR/$REFALG/data_f${fun}/$f
          done
        )

        # take the comprehensive info and add it to bbobexp_f#.info
        # strip the additional CMA-ES entries (above the first 15 entries)
        sed -n "s/'$EXPID[^']*'/'${REFALG}'/;s/\([^,]\+\)\(\(,[^,]\+\)\{15\}\),.*/\1\2/;4,6p" ${dir}/bbobexp_f${fun}.info >> $OUTPUT_DIR/$REFALG/bbobexp_f${fun}.info
        # newline :)
        # echo "" >> $OUTPUT_DIR/$REFALG/bbobexp_f${fun}.info
      fi
    done
  done
done

