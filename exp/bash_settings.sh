# Global settings and variables for surrogate-cmaes project

# MATLAB Runtime environment
module add matlab
MATLAB_PATH=`dirname \`which matlab\``/..
export LD_LIBRARY_PATH=${MATLAB_PATH}/runtime/glnxa64:${MATLAB_PATH}/bin/glnxa64:${MATLAB_PATH}/sys/os/glnxa64:/storage/plzen1/home/bajeluk/bin/mcr_2016a/v901/runtime/glnxa64:/storage/plzen1/home/bajeluk/bin/mcr_2016a/v901/bin/glnxa64:/storage/plzen1/home/bajeluk/bin/mcr_2016a/v901/sys/os/glnxa64:$LD_LIBRARY_PATH

# set SCRACHDIR if it is not set (outside Metacentrum)
if [ "${SCRATCHDIR}XX" = "XX" ]; then
  SCRATCHDIR=/tmp/surrogate-cmaes_$$
  mkdir -p "$SCRATCHDIR"
fi

# Directory where to put the deployed files and where to run
# the experiment from
RUNDIR="$SCRATCHDIR/surrogate-cmaes"

# Binary which is to be called as the individual jobs
MATLAB_BINARY_CALL="exp/metacentrum_task_matlab"
export MCR_CACHE_ROOT=$SCRATCHDIR

# Deployment archive of Matlab source codes (for running and compilation
# on a target machine)
DEPLOY_DIR=deploy
DEPLOY_ARCHIVE=${EXPID}_src.tar

# Files and directories which should be packed into the deployment package
FILES_TO_DEPLOY="exp/*.m exp/*.sh exp/experiments/*.m exp/util exp/log exp/vendor/bbob exp/vendor/saACMES src/ exp/pproc/generateGnuplot* Makefile startup.m"

# allow read-access of newly created files and directories for the group
umask 0027
