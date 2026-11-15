#!/bin/bash
# One chunk of the wide-truth bias test: generate toys from the wide-range refit of one function (truth.json), then
# fit the same toys with the nominal envelope twice:
#   narrow : the parameter ranges stored in the workspace, as in the expected limits
#   wide   : the options of Combine/Checks/Bias_nominal/bias_fit_config.sh, as in the AN bias study (Sec. 7.5),
#            which profile with the wider turn-on/width/sigma ranges
# Usage : runjob.sh MASS INJ CHUNK NTOYS     INJ = 0 or q50
# Input : Combine/root_t2w/<m>_Datacard_leptons.root (cross-section card, read only), truth.json
# Output: out/mA<m>_<INJ>/chunk<c>_{narrow,wide}.root (MultiDimFit trees), done marker chunk<c>.done
set -eo pipefail
m="$1"; inj="$2"; c="$3"; n="$4"
H=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/combine_review_20261105/bias_wide_truth
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
card="${FF}/Combine/root_t2w/${m}_Datacard_leptons.root"
out="${H}/out/mA${m}_${inj}"; mkdir -p "$out"
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"
ulimit -s unlimited
wd="${_CONDOR_SCRATCH_DIR:-${H}/work}/mA${m}_${inj}_c${c}_$$"; mkdir -p "$wd"; cd "$wd"
read -r idx setp rng q50 < <(python3 -c "
import json; t=json.load(open('${H}/truth.json'))['${m}']; print(t['index'], t['set'], t['ranges'], t['q50'])")
[[ "$inj" == "q50" ]] && r_inj="$q50" || r_inj=0
[[ "$inj" == "q50" ]] && seed=$(( 700000 + m*1000 + 500 + c )) || seed=$(( 700000 + m*1000 + c ))
echo "mA=${m} inj=${inj} r=${r_inj} chunk=${c} ntoys=${n} seed=${seed} truth index=${idx}"
combine -M GenerateOnly -d "$card" -m 125.38 --expectSignal "$r_inj" --saveToys -t "$n" -s "$seed" -n _gen \
  --setParameters "MH=125.38,pdfindex_cat0_13p6TeV=${idx},${setp}" --freezeParameters MH,pdfindex_cat0_13p6TeV \
  --setParameterRanges "$rng" > gen.log 2>&1
toys=$(ls higgsCombine_gen.GenerateOnly.mH125.38.${seed}.root)
source "${FF}/Combine/Checks/Bias_nominal/bias_fit_config.sh"
wide_opts="$(wide_bias_fit_combine_options)"
# narrow: identical minimizer options, but only the r range is set, so the stored parameter ranges apply
narrow_opts="${wide_opts%%--setParameterRanges*}--setParameterRanges r=-5,5"
for cfg in narrow wide; do
  [[ "$cfg" == narrow ]] && opts="$narrow_opts" || opts="$wide_opts"
  eval combine -M MultiDimFit -d "$card" -m 125.38 -P r --algo singles $opts -t "$n" -n "_${cfg}" --toysFile="$toys" > "fit_${cfg}.log" 2>&1
  cp -p higgsCombine_${cfg}.MultiDimFit.mH125.38.*.root "${out}/chunk${c}_${cfg}.root"
done
date > "${out}/chunk${c}.done"
echo "[OK] mA=${m} inj=${inj} chunk=${c}"
