#!/usr/bin/env bash
# Staging build for opening the step/smearing ranges at mA10, 15, 21, 22, 23 (wide-truth bias failed there with
# the per-mass windows; same rule as mA24-29). Production Background/ is not touched.
#  - library: HEAD Background/src + interface, with openStepGausRanges() extended to the 5 masses. HEAD differs
#    from the 10-01 production library only by massALP == 0 (merged) branches, checked with git diff d61f2fa e5db9b4.
#  - fTest: built from the same source as the production binary (deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp),
#    rpath -> the staging library (the production binary has DT_RPATH, so LD_LIBRARY_PATH cannot redirect it).
# Input : flashggFinalFit/Background/{src,interface,makefile}, deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp
# Output: Background/{lib,obj,bin} under this directory; logs/build.log
set -eo pipefail
D=$(cd "$(dirname "$0")" && pwd)
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
REV=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_envrevert_20261023/src/fTest_ALP_turnOn.revert.cpp
S=$D/Background
mkdir -p $D/logs $S/lib $S/obj $S/bin
cd $F/.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval $(scram runtime -sh)
rm -rf $S/src $S/interface
cp -rp $F/Background/src $F/Background/interface $F/Background/makefile $S/
python3 - "$S/src/PdfModelBuilder.cc" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
old = 'return massALP == 0 || (massALP >= 24 && massALP <= 29); }'
new = 'return massALP == 0 || massALP == 10 || massALP == 15 || (massALP >= 21 && massALP <= 29); }'
assert s.count(old) == 1, 'openStepGausRanges line not found exactly once'
open(p, 'w').write(s.replace(old, new))
EOF
grep -n 'inline bool openStepGausRanges' $S/src/PdfModelBuilder.cc
( cd $S && make -f makefile lib ) > $D/logs/build.log 2>&1
( cd $S && make -f makefile vars ) > $D/logs/make_vars.txt
CXXF=$(sed -n 's/^CXXFLAGS: *//p' $D/logs/make_vars.txt)
LIBS=$(sed -n 's/^LIBS: *//p' $D/logs/make_vars.txt)
sed "s#\"../../tdrStyle/#\"$F/tdrStyle/#" $REV > $S/fTest_ALP_turnOn.build.cpp
g++ -o $S/bin/fTest_ALP_turnOn $CXXF $S/fTest_ALP_turnOn.build.cpp -L$S/lib -lBackgroundProfileFitting $LIBS >> $D/logs/build.log 2>&1
readelf -d $S/bin/fTest_ALP_turnOn | grep -i rpath
ls -l $S/lib $S/bin
