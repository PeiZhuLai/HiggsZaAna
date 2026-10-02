#!/usr/bin/env bash
# Sequencer for the 2024 DY+jets overlap-veto rerun after stages A-E:
#   1. wait for hza-dyveto-retrain, require VERDICT: PASSED
#   2. dataVmc sideband_rwgt  (chain_dyveto_datavmc.sh; ALP_Optimization reads it)
#   3. flashgg                (chain_dyveto_flashgg.sh: WP -> S1 lnN -> ... -> limits)
#   4. dataVmc nominal        (App F "before" panels)
# Stops at the first failure. Not here: mA3 closure WP scan, impacts/bias (condor), AN.
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
HB=$L/dyveto_driver_heartbeat.txt
REG=$L/PENDING_RESULTS.txt
hb()  { echo "$(date '+%F %T') $*" | tee -a "$HB"; }
stop(){ hb "STOPPING: $*"; echo "$(date '+%F %H:%M') STOPPED dyveto-driver: $*" >> "$REG"; exit 1; }
passed(){ grep -q "^VERDICT: PASSED" "$1" 2>/dev/null; }

hb "driver start (pid $$), waiting for hza-dyveto-retrain"
while [ "$(systemctl --user is-active hza-dyveto-retrain 2>/dev/null)" = "active" ]; do sleep 120; done
passed $L/dyveto_retrain.log || stop "retrain A-E did not pass (dyveto_retrain.log)"
hb "retrain A-E passed"

hb "dataVmc sideband_rwgt start"
FT=sideband_rwgt bash $D/chain_dyveto_datavmc.sh > $L/dyveto_datavmc_sideband_rwgt.log 2>&1
passed $L/dyveto_datavmc_sideband_rwgt.log || stop "dataVmc sideband_rwgt (dyveto_datavmc_sideband_rwgt.log)"
hb "dataVmc sideband_rwgt passed"

hb "flashgg start"
bash $D/chain_dyveto_flashgg.sh > $L/dyveto_flashgg.log 2>&1
passed $L/dyveto_flashgg.log || stop "flashgg (dyveto_flashgg.log)"
hb "flashgg passed"
echo "$(date '+%F %H:%M') DONE dyveto-flashgg (log dyveto_flashgg.log). Next: mA3 closure WP scan, impacts/bias, closure/plots, AN." >> "$REG"

hb "dataVmc nominal start"
FT=nominal bash $D/chain_dyveto_datavmc.sh > $L/dyveto_datavmc_nominal.log 2>&1
passed $L/dyveto_datavmc_nominal.log || stop "dataVmc nominal (dyveto_datavmc_nominal.log)"
hb "dataVmc nominal passed"
hb "DRIVER DONE"
