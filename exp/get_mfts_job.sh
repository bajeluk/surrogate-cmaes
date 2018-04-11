#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb:cgroups=cpuacct

trap 'clean_scratch' TERM EXIT

module add matlab
cd ~/src/surrogate-cmaes

matlab -singleCompThread -nosplash -nodesktop -r 'getDataMetaFeatures; exit';

exit;
