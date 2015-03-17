#!/bin/sh
#PBS -l walltime=1d
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l scratch=1gb
#PBS -l matlab=1,matlab_Statistics_Toolbox=1,matlab_Optimization_Toolbox=1

HOSTNAME=optim.wz.cz
USERNAME=optim.wz.cz
PASS=metacentrum
LOGFILENAME=log.txt
LOGFILE="$LOGDIR/$LOGFILENAME"
EXPPATH="$EXPPATH_SHORT/$EXPID"
OUTPUTDIR="$SCRATCHDIR"
RESULTSFILE="$EXPPATH/${EXPID}_results_${FUN}_${DIM}D_${ID}.mat"
FILESTDOUT="$OUTPUTDIR/${EXPID}__log__${ID}.txt"
FILEMANAGER="${EXPPATH}/${EXPID}_manager.sh"
MATLABCALL="matlab"
MATLABPARAMS="-singleCompThread -nodisplay -nodesktop"
MACHINE=`head -1 $PBS_NODEFILE`

module add matlab

if [ -z "$EXPID" ] ; then
  echo "Error: EXPID (experiment ID) is not known"; exit 1
fi
if [ -z "$ID" ] ; then
  echo "Error: Task ID number is not known"; exit 1
fi
if [ -z "$EXPPATH_SHORT" ] ; then
  echo "Error: directory with the experiment is not known"; exit 1
elif [ ! -d "$EXPPATH" ] ; then
  echo "Error: directory with the experiment does not exist"; exit 1
fi
if [ -z "$LOGFILE" ] ; then
  LOGFILE="/storage/plzen1/home/bajeluk/public/log.txt"
  echo "Note: default log file set"
fi

# clean up the lock-file
trap "rm -f $EXPPATH/queued_$ID" TERM EXIT

testMatlabFinished () {
  if [ "$3" -eq 0  -a  -f "$1"  -a  "$1" -nt "$FILEMANAGER" ]; then
    echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] $2 succeeded." >> $LOGFILE
    echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] $2 succeeded." >> $FILESTDOUT
  else
    echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] $2 !!!  FAILED !!!" >> $LOGFILE
    echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] $2 !!!  FAILED !!!" >> $FILESTDOUT
    # tail -n 60 "$FILESTDOUT" | mail -s "Metacentrum: chyba v uloze $PBS_JOBID $PBS_JOBNAME" $USER@arien.metacentrum.cz
  fi
  ftp -n $HOSTNAME <<EOD
user $USERNAME $PASS
lcd $LOGDIR
put $LOGFILENAME
bye
EOD
}

cd "$EXPPATH_SHORT/../.."; ulimit -t unlimited;

echo "###########################################" > $FILESTDOUT
echo "     Matlab call id=${ID}" > $FILESTDOUT
echo "" > $FILESTDOUT
echo "  dim(s): $DIM    f(s): $FUN    N(inst): $NINSTANCES" > $FILESTDOUT
echo "  model: $MODEL" > $FILESTDOUT
echo "###########################################" > $FILESTDOUT

echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] $2 started." >> $LOGFILE

echo "exit(1)" | $MATLABCALL $MATLABPARAMS -r "bbob_test_01($ID, '$EXPID', '$EXPPATH_SHORT'), exit(0)"  > $FILESTDOUT 2>&1
testMatlabFinished "$RESULTSFILE" $ID $?

cp $OUTPUTDIR/* $EXPPATH/ || export CLEAN_SCRATCH=false

# echo `date "+%Y-%m-%d %H:%M:%S"` "  **$EXPID** at [$MACHINE] ==== FINISHED ====" >> $LOGFILE
