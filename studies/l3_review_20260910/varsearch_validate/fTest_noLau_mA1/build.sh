#!/usr/bin/env bash
# Builds the sandbox fTest_ALP_turnOn (Laurent removed at m_a = 1) against the production Background library (read only).
# Flags = production makefile "make vars" (logs/make_vars.txt); rpath -> production lib.
set -eo pipefail
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
CXXF=$(sed -n 's/^CXXFLAGS: *//p' /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/fTest_noLau_mA1/logs/make_vars.txt | sed 's#-Wl,-rpath=[^ ]*#-Wl,-rpath=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/lib#')
LIBS=$(sed -n 's/^LIBS: *//p' /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/fTest_noLau_mA1/logs/make_vars.txt)
g++ -o /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/fTest_noLau_mA1/bin/fTest_ALP_turnOn $CXXF /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/fTest_noLau_mA1/test/fTest_ALP_turnOn.cpp -L/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/lib -lBackgroundProfileFitting $LIBS
