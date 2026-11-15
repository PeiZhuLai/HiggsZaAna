#!/bin/bash
# One HybridNew grid point (LHC-style CLs) for the asymptotic-vs-toys check of the Combine review.
# The data are replaced by the Asimov b-only dataset (-D ...:toys/toy_asimov), so the observed
# data are never used, not even for the nuisance-parameter profiling of the toy generation.
# Usage: runjob_hybridnew.sh MASS R SEED NTOYS [SUBDIR] [extra combine options...]
m=$1; r=$2; seed=$3; ntoys=$4; sub=${5:-hybridnew}; shift 5 2>/dev/null; extra="$*"
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
HNA=$FF/Combine/Checks/review_20261105/hybridnew
OUT=$FF/Combine/Checks/review_20261105/$sub
ASIMOV=$HNA/higgsCombine_asimov_mA${m}.GenerateOnly.mH125.38.123456.root
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` 2>/dev/null
cd "${_CONDOR_SCRATCH_DIR:-$(mktemp -d -p $OUT)}" || exit 2
name=_mA${m}_r${r}_s${seed}
combine -M HybridNew $FF/Combine/root_t2w_xsec/${m}_Datacard_leptons.root -m 125.38 \
  -D ${ASIMOV}:toys/toy_asimov --setParameters MH=125.38 --freezeParameters MH \
  --LHCmode LHC-limits --singlePoint $r --saveHybridResult --clsAcc 0 -T $ntoys -s $seed \
  --cminDefaultMinimizerStrategy 0 $extra -n $name || exit 3
mkdir -p $OUT/mA$m
f=$(ls higgsCombine${name}.HybridNew.mH125.38.${seed}.root) || exit 4
cp $f $OUT/mA$m/ && echo "[OK] $OUT/mA$m/$f"
