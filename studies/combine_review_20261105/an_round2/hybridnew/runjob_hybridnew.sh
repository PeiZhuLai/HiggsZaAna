#!/bin/bash
# One HybridNew grid point (LHC-style CLs), xsec card, current production workspace.
# The data are replaced by the b-only Asimov dataset (-D ...:toys/toy_asimov), so the observed data are never used,
# not even for the nuisance profiling of the toy generation. Same command as the earlier review check
# (studies/combine_review_20261105/runjob_hybridnew.sh); only the input/output directories differ.
# Usage: runjob_hybridnew.sh MASS R SEED NTOYS
m=$1; r=$2; seed=$3; ntoys=$4
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
R2=$FF/Combine/Checks/review_round2
ASIMOV=$R2/asimov/higgsCombine_asimov_mA${m}.GenerateOnly.mH125.38.123456.root
OUT=$R2/hybridnew
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` 2>/dev/null
cd "${_CONDOR_SCRATCH_DIR:?no condor scratch}" || exit 2
name=_mA${m}_r${r}_s${seed}
combine -M HybridNew $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 \
  -D ${ASIMOV}:toys/toy_asimov --setParameters MH=125.38 --freezeParameters MH \
  --LHCmode LHC-limits --singlePoint $r --saveHybridResult --clsAcc 0 -T $ntoys -s $seed \
  --cminDefaultMinimizerStrategy 0 --fullBToys -n $name || exit 3
mkdir -p $OUT/mA$m
f=$(ls higgsCombine${name}.HybridNew.mH125.38.${seed}.root) || exit 4
cp $f $OUT/mA$m/ && echo "[OK] $OUT/mA$m/$f"
