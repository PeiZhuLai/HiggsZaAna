#!/bin/bash
# mA10 xsec card, r = q50, fixed pdfindex (IMPACT_FIX_INDEX=1), production ranges; strategy 0 and 1 side by side.
# Output: Combine/output_impacts_xsec_q50_fixidx_prodrange{,_s1}/10_impacts.json, logs here
R=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/shellScripts/impact/Condor/runjob_expectedImpact.sh
W="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for s in 0 1; do
  env -i HOME=$HOME USER=$USER PATH=/usr/bin:/bin KRB5CCNAME=$KRB5CCNAME IMPACT_CARD=xsec IMPACT_EXPECT=q50 IMPACT_SKIP_T2W=1 \
    IMPACT_FIX_INDEX=1 IMPACT_RANGES=prod IMPACT_STRATEGY=$s IMPACT_PARALLEL=3 bash "$R" 10 > "$W/mA10_s$s.log" 2>&1 &
done
wait; echo "[$(date +%H:%M:%S)] ALL DONE"
