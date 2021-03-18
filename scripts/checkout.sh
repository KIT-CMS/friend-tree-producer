#!/bin/bash


### Setup of CMSSW release
NUM_CORES=10
CMSSW=CMSSW_10_2_14

if [ "$1" == "" ]; then
  echo "$0: Explicit CMSSW version is not provided. Checking out as default $CMSSW"
fi

if [ ! "$1" == "" ]; then
  CMSSW=$1
  echo "$0: Checking out $CMSSW"
fi


export SCRAM_ARCH=slc6_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project $CMSSW; pushd $CMSSW/src
eval `scramv1 runtime -sh`

### Checkout of external software

# SVfit and fastMTT
git clone git@github.com:SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_19_02_2019
git clone git@github.com:SVfit/SVfitTF TauAnalysis/SVfitTF

# FF weights
git clone ssh://git@github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes
cd HTTutilities/Jet2TauFakes
git checkout v0.2.2
cp -r /ceph/htautau/fakefactor_files/nmssm/data_201* . 

git clone git@github.com:KIT-CMS/fake-factor-application.git -b nmssm_analysis HiggsAnalysis/fake-factor-application

# Tau Trigger
git clone git@github.com:KIT-CMS/TauTriggerSFs.git TauAnalysisTools/TauTriggerSFs -b run2_SFs_TriggerFitsForEmbedded_DeepTau_SingleTau

# HH kinematic fitting
git clone git@github.com:janekbechtel/HHKinFit.git HHKinFit/HHKinFit

### Checkout of friend tree producer setup
git clone git@github.com:KIT-CMS/friend-tree-producer.git HiggsAnalysis/friend-tree-producer -b nmssm_analysis
git clone git@github.com:KIT-CMS/grid-control
# Data sources
mkdir HiggsAnalysis/friend-tree-producer/data/input_params
cd HiggsAnalysis/friend-tree-producer/data/input_params
wget https://raw.githubusercontent.com/KIT-CMS/datasets/master/datasets.json -b nmssm
cd -

### Compiling under CMSSW
USER_CXXFLAGS="-Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j $NUM_CORES
popd
