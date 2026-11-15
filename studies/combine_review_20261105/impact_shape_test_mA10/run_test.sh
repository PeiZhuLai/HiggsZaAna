#!/bin/bash
# Why do the background shape parameters get no interval in the expected impacts? Test at mA10, xsec card,
# Asimov with r = q50, same options as runjob_expectedImpact.sh except for what each test changes:
#   free       discrete index profiled (the production setting)
#   idx0       index frozen at 0 (Exp3, the function the Asimov dataset is generated from)
#   idx0_wide  as idx0, and the Exp3 ranges widened like those of the other families
#   idx1       index frozen at 1 (Pow1)
# Input : Combine/root_t2w_xsec/10_Datacard_leptons.root (read only)
# Output: <test>/impacts.json and the combine logs next to this script. No observed data is fitted.
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src
cd "$SRC" && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scramv1 runtime -sh)" 2>/dev/null
CARD=$SRC/flashggFinalFit/Combine/root_t2w_xsec/10_Datacard_leptons.root
Q50=0.0199219
RANGES='Pow_turnon_p1=95,125:Pow_width_p1=0.1,50:Pow_sigma_p1=0.05,20:Lau_turnon_p1=95,125:Lau_width_p1=0.1,50:Lau_sigma_p1=0.05,20:Bern_gsigma=0.05,20:Bern_step=90,130:Bern_stepWidth=0.1,50:lumi_13p6TeV_Uncorrelated_2022preEE=-5,5'
WIDE=':Exp_turnon_p3=95,125:Exp_width_p3=0.1,50:Exp_sigma_p3=0.05,20'
EXP=Exp_turnon_p3,Exp_width_p3,Exp_sigma_p3,Exp_p1_exp3,Exp_p3_exp3,Exp_cp3_exp3
POW=Pow_turnon_p1,Pow_width_p1,Pow_sigma_p1,Pow_p1_pow1
NORM=shapeBkg_bkg_mass_cat0__norm,CMS_hza_pileup

run() {  # name, frozen parameters, extra --setParameters or "", ranges, parameter list
  local name=$1 frozen=$2 setp=$3 ranges=$4 named=$5
  mkdir -p "$W/$name"; cd "$W/$name" || return 1
  local opts=(-M Impacts -d "$CARD" -m 125.38 --redefineSignalPOIs r --rMin -5 --rMax 5 -t -1 --expectSignal "$Q50"
    --robustFit 1 --cminPreScan --freezeParameters "$frozen"
    --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants
    --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2
    --setParameterRanges "$ranges"
    --cminDefaultMinimizerType Minuit2 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01
    --cminFallbackAlgo Minuit2,0:0.1 --cminFallbackAlgo Minuit2,Migrad,1:0.1
    --cminFallbackAlgo Minuit2,Simplex,0:0.1 --cminFallbackAlgo Minuit2,Migrad,0:1.0 -v 1 --named "$named")
  [[ -n "$setp" ]] && opts+=(--setParameters "$setp")
  echo "[$(date +%H:%M:%S)] $name: initial fit"
  combineTool.py "${opts[@]}" --doInitialFit > initial.log 2>&1; echo "  initial rc=$?"
  echo "[$(date +%H:%M:%S)] $name: fits"
  combineTool.py "${opts[@]}" --doFits --parallel 4 > fits.log 2>&1; echo "  fits rc=$?"
  combineTool.py -M Impacts -d "$CARD" -m 125.38 --redefineSignalPOIs r --named "$named" -o impacts.json > json.log 2>&1; echo "  json rc=$?"
}

run free      MH                           ""                         "$RANGES"       "$EXP,$POW,$NORM"
run idx0      MH,pdfindex_cat0_13p6TeV     pdfindex_cat0_13p6TeV=0    "$RANGES"       "$EXP,$NORM"
run idx0_wide MH,pdfindex_cat0_13p6TeV     pdfindex_cat0_13p6TeV=0    "$RANGES$WIDE"  "$EXP,$NORM"
run idx1      MH,pdfindex_cat0_13p6TeV     pdfindex_cat0_13p6TeV=1    "$RANGES"       "$POW,$NORM"
echo "[$(date +%H:%M:%S)] ALL DONE"
