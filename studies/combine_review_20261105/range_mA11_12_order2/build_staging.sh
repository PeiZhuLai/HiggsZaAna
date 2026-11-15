#!/usr/bin/env bash
# Side study (Combine review, mA11 and mA12): the Pow1 truth fails the bias requirement at both points (mA11 -0.217,
# deployed 1 GeV floor; mA12 -0.240, production) and the step width / smearing floor does not move it
# (range_mA12_floor{1p5,2,3}, range_mA12_fixsig2). The failing shape is the steep minimum of the erf x m^p form
# (p ~ -9.2), reachable only by the PowerLaw exponent range. Here the envelope additionally gets the next order of
# the Exponential, PowerLaw and Laurent families (Exp3, Pow3, Lau2) when its sideband GOF passes; Pow3 holds two free
# exponents and so both the steep and the shallow shape.
# Staging build only; production Background/ is not touched.
#  - library: copy of production Background/{src,interface,makefile}, unchanged (mA11 1 GeV floor, mA12 no floor).
#  - fTest: source of the production binary (../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp) + ORDER2_MASSES switch
#    (src/patch_ftest.py); rpath -> staging library.
# Input : flashggFinalFit/Background/{src,interface,makefile}, ../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
# Output: src/fTest_ALP_turnOn.order2.cpp; Background/{src,interface,lib,obj,bin} under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
SRC0=$D/../deploy_mA1_2/src/fTest_ALP_turnOn.deploy.cpp
SRC=$D/src/fTest_ALP_turnOn.order2.cpp
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
