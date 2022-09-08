#!/bin/bash

ulimit -s unlimited

executable=$1
cores=$2
custom_workdir_path=$3
mode=$4
path_to_script=$5
eval $(scram unsetenv -sh)

source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
echo "Running python $path_to_script --executable $executable --cores $cores --custom_workdir_path $custom_workdir_path --mode $mode"
python3 $path_to_script --executable $executable --cores $cores --custom_workdir_path $custom_workdir_path --mode $mode