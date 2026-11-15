#!/usr/bin/env bash
# Multi-start fTest on top of the staging library (opened ranges at mA10, 15, 21-29, built by build_staging.sh).
# Same source as the production binary (deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp) plus
# patch_databinned.py (data_obs from the fit-window dataset), patch_multistart.py, patch_extraenv.py (EXTRA_ENV_PDFS) and patch_excludeenv.py (EXCLUDE_ENV_PDFS), no effect when unset; rpath -> the staging library. The library is not rebuilt.
# Input : deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp, Background/lib (staging), logs/make_vars.txt
# Output: Background/fTest_ALP_turnOn.multistart.cpp, Background/bin/${OUTBIN:-fTest_ALP_turnOn.ms}, logs/build_ms.log
#         (OUTBIN lets a new build go next to a binary that is still running)
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
REV=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp
S=$D/Background
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
python3 $D/patch_databinned.py $REV $S/fTest_ALP_turnOn.databinned.cpp
python3 $D/patch_multistart.py $S/fTest_ALP_turnOn.databinned.cpp $S/fTest_ALP_turnOn.ms_only.cpp
python3 $D/patch_extraenv.py $S/fTest_ALP_turnOn.ms_only.cpp $S/fTest_ALP_turnOn.msx_only.cpp
python3 $D/patch_excludeenv.py $S/fTest_ALP_turnOn.msx_only.cpp $S/fTest_ALP_turnOn.multistart.cpp
sed -i "s#\"../../tdrStyle/#\"$F/tdrStyle/#" $S/fTest_ALP_turnOn.multistart.cpp
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/${OUTBIN:-fTest_ALP_turnOn.ms} $CXXF $S/fTest_ALP_turnOn.multistart.cpp -L$S/lib -lBackgroundProfileFitting $LIBS > $D/logs/build_ms.log 2>&1
readelf -d $S/bin/${OUTBIN:-fTest_ALP_turnOn.ms} | grep -i rpath
ls -l $S/bin
