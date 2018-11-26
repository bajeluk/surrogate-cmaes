#!/usr/bin/env bash

EXPPATH=$1

design=ilhs

for modeltype in gp forest xgb; do
  if [ $modeltype = gp ]; then
    optinds=`seq 1 14`
  elif [ $modeltype = forest ]; then
    optinds=`seq 1 15`
  elif [ $modeltype = xgb ]; then
    optinds=`seq 1 5`
  fi

  for dim in 2 5 10; do
    N=$(( 50 * $dim ))

    for fun in `seq 1 24`; do

      for optid in $optinds; do
        fname="$EXPPATH/results/${dim}D/f${fun}/res_${dim}D_f${fun}_inst0_N${N}_design-${design}_model-${modeltype}_opts${optid}.mat"
        if [ ! -f $fname ]; then
          echo missing dim:${dim} fun:${fun} model:${modeltype} optind:${optid}
        fi
      done 
    done

#    if [ $ok -gt 0 ]; then
#      echo $optid
#    fi
  done
done
