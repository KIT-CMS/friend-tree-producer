#!/bin/bash

executable=$1 
cores=$2
custom_workdir_path=$3
mode=$4
eval $(scram unsetenv -sh)

source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
echo "Running python HiggsAnalysis/friend-tree-producer/scripts/collect_with_RDF.py --executable $executable --cores $cores --custom_workdir_path $custom_workdir_path --mode $mode"
python HiggsAnalysis/friend-tree-producer/scripts/collect_with_RDF.py --executable $executable --cores $cores --custom_workdir_path $custom_workdir_path --mode $mode