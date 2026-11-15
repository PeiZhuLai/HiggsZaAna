#!/usr/bin/env bash
# Side study (Combine review, mA11 and mA12, option 2): the stored Pow1 is a local minimum. The Lau1 solution lies
# inside the Pow1 ranges with a lower NLL (mA11 478.975 vs 479.81, mA12 525.625 vs 527.117) and the multi-start only
# scans turn-on x width. Here the multi-start also starts the exponent of Pow1, Lau1 and Exp1 from 15/50/85% of its
# range (EXPSCAN_MASSES switch, src/patch_ftest.py); the nominal 15 starts are kept.
# Staging build only; production Background/ is not touched.
#  - library: copy of production Background/{src,interface,makefile}, unchanged (mA11 1 GeV floor, mA12 no floor).
#  - fTest: source of the production binary (../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp) + EXPSCAN_MASSES switch
#    (src/patch_ftest.py); rpath -> staging library.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
# Output: src/fTest_ALP_turnOn.pscan.cpp; Background/{src,interface,lib,obj,bin} under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC0=$D/../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
SRC=$D/src/fTest_ALP_turnOn.pscan.cpp
S=$D/Background
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
python3 $D/src/patch_ftest.py $SRC0 $SRC
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
cmp $S/src/PdfModelBuilder.cc $F/Background/src/PdfModelBuilder.cc
grep -n 'wsFloor = ' $S/src/PdfModelBuilder.cc
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $SRC -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
ls -l $S/lib $S/bin
