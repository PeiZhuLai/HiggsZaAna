#!/usr/bin/env bash
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
while true; do
  a=$(systemctl --user is-active hza-flashgg 2>/dev/null); b=$(systemctl --user is-active hza-datavmc2 2>/dev/null)
  {
    echo "[$(date '+%F %T')] host=$(hostname -s) datavmc2=$b flashgg=$a"
    echo "  datavmc2 stage: $(grep -o '\[[0-9]/6\][^=]*' $D/logs_fsrfix/datavmc2.log 2>/dev/null | tail -1)"
    echo "  flashgg stage : $(grep -o '\[[0-9]/9\][^=]*' $D/logs_fsrfix/flashgg.log 2>/dev/null | tail -1)"
    echo "  flashgg last  : $(tail -1 $D/logs_fsrfix/flashgg.log 2>/dev/null | cut -c1-120)"
    echo "  bkg fits done : $(ls $D/../../Background 2>/dev/null >/dev/null; find /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3 -name 'CMS-HGG_mva_13p6TeV_multipdf.root' -newer $D/logs_fsrfix/datavmc2.log 2>/dev/null | wc -l)/30"
  } > $D/logs_fsrfix/flashgg_heartbeat.txt
  [ "$a" != "active" ] && { echo "  FINISHED" >> $D/logs_fsrfix/flashgg_heartbeat.txt; break; }
  sleep 600
done
