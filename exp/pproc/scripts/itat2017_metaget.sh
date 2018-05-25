#!/bin/bash

myrsync_copy='rsync -aiOvh --out-format="%i  %M %9'\'''\'''\''l  %n"'

cd ~/prg/surrogate-cmaes/

MATFILES_ONLY="--include='*results*.mat' --include='scmaes_params.mat' --exclude='*'"
HOST=bajeluk@alfrid.meta.zcu.cz

EXPID=exp_doubleEC_23_adapt05_v3
USER=bajeluk
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/
USER=pitrazby
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/

EXPID=exp_doubleEC_23_adapt05_v2
USER=juranja3
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/

EXPID=exp_doubleEC_23_adapt05
USER=pitrazby
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/
USER=bajeluk
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/
USER=holena
eval $myrsync_copy $MATFILES_ONLY ${HOST}:../$USER/prg/surrogate-cmaes/exp/experiments/$EXPID/ exp/experiments/$EXPID/
