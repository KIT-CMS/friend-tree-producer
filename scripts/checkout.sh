#!/bin/bash


### Setup of CMSSW release
NUM_CORES=10
CMSSW=CMSSW_10_2_22

if [ "$1" == "" ]; then
  echo "$0: Explicit CMSSW version is not provided. Checking out as default $CMSSW"
fi

if [ ! "$1" == "" ]; then
  CMSSW=$1
  echo "$0: Checking out $CMSSW"
fi

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project $CMSSW; pushd $CMSSW/src
eval `scramv1 runtime -sh`

### Checkout of external software

# SVfit and fastMTT
git clone git@github.com:SVfit/ClassicSVfit TauAnalysis/ClassicSVfit -b fastMTT_19_02_2019
git clone git@github.com:SVfit/SVfitTF TauAnalysis/SVfitTF

# MELA
git clone git@github.com:cms-analysis/HiggsAnalysis-ZZMatrixElement ZZMatrixElement -b v2.2.3
# NOTE: The following two lines should be needed following the wiki, but it
# seems to work either way. However, these lines introduce a dependency on AFS.
# Wiki: https://twiki.cern.ch/twiki/bin/view/CMS/MELAProject#Checkout_instructions
#cp ZZMatrixElement/MELA/data/mcfm.xml ../config/toolbox/$SCRAM_ARCH/tools/selected/
#scram setup mcfm
cd ZZMatrixElement/
bash setup.sh -j $NUM_CORES
cd ..

# TODO NN mass

# FF weights
git clone ssh://git@github.com/CMS-HTT/Jet2TauFakes.git HTTutilities/Jet2TauFakes
cd HTTutilities/Jet2TauFakes
git checkout v0.2.2

read -p "lxplus-username: " USERNMLXP

scp -r ${USERNMLXP}@lxplus.cern.ch:/afs/cern.ch/user/j/jbechtel/public/fake-factor-files/* ./

cd ../..
git clone git@github.com:KIT-CMS/fake-factor-application.git HiggsAnalysis/fake-factor-application

cd -
# TODO NN MET

# TODO NN max score

# TODO single-tau HLT

# Tau Trigger
git clone git@github.com:KIT-CMS/TauTriggerSFs.git TauAnalysisTools/TauTriggerSFs -b run2_SFs_TriggerFitsForEmbedded_DeepTau_SingleTau

# HH kinematic fitting
git clone git@github.com:janekbechtel/HHKinFit.git HHKinFit/HHKinFit

### Checkout of friend tree producer setup
git clone git@github.com:KIT-CMS/friend-tree-producer.git HiggsAnalysis/friend-tree-producer
git clone git@github.com:KIT-CMS/grid-control
# Data sources
mkdir HiggsAnalysis/friend-tree-producer/data/input_params
cd HiggsAnalysis/friend-tree-producer/data/input_params
scp ${USERNMLXP}@lxplus.cern.ch:/eos/home-s/swozniew/friend-tree-producer-input-params/* ./
wget https://raw.githubusercontent.com/KIT-CMS/datasets/master/datasets.json
cd -
# Imperial FF weights
mkdir -p HiggsAnalysis/friend-tree-producer/data/imperial_ff
cd HiggsAnalysis/friend-tree-producer/data/imperial_ff
scp -r ${USERNMLXP}@lxplus.cern.ch:/afs/cern.ch/user/g/guttley/public/fake_factors_mssm/fakefactors_ws*_v2.root ./


### Compiling under CMSSW
USER_CXXFLAGS="-Wno-delete-non-virtual-dtor -Wno-error=unused-but-set-variable -Wno-error=unused-variable" scram b -j $NUM_CORES
popd
