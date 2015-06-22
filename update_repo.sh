#!/bin/bash

  echo rm -rf .git doc
rm -rf .git doc
  echo cd ..
cd ..
  echo mv surrogate-cmaes surrogate-cmaes_OLD
mv surrogate-cmaes surrogate-cmaes_OLD
  echo git clone git@github.com:bajeluk/surrogate-cmaes.git
git clone git@github.com:bajeluk/surrogate-cmaes.git
  echo cd surrogate-cmaes_OLD
cd surrogate-cmaes_OLD
  echo cp -r * ../surrogate-cmaes/
cp -r * ../surrogate-cmaes/
  echo cd ..
cd ..
  echo rm -rf surrogate-cmaes_OLD
rm -rf surrogate-cmaes_OLD
  echo cd surrogate-cmaes
cd surrogate-cmaes
  echo git checkout exp/experiments/exp_BIPOP_saACMES_02.m exp/make_bbob_graphs.sh
git checkout exp/experiments/exp_BIPOP_saACMES_02.m exp/make_bbob_graphs.sh
git status
