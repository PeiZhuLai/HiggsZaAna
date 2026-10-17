#!/usr/bin/env bash
# Step F (plots): regenerate the Plot/ figures and tables of the AN that depend on the m_a 1-3 BDT scores / working
# points, after the production MVA cut of the low_relE2_w123 model. At most 3 at a time (the score patch runs
# beside it with 2; lxplus limit 6). Same environment as an_update_sigrw_20261003/run_plots_sigrw.sh.
#   scan_score_R_significance.py (low-mass R/Z scan, red star = pinned w123 cuts; no --write-json)
#   make_sculpt_R_table_from_wp.py -> Plot/output/latexTable_sculptR.txt
#   plot_MVASigEffVmA.py, plot_mAmigrated{Bar,Hist,Matrix}.py (HZA_SIGNAL_REWEIGHT=1, as for the AN)
#   plot_mA1_doublepeak.py, plot_signal_reweight_{trainvars,mllgg}_bdtwp.py --ma 1,2,3
#   plot_bkgmcSculptingCheck.py, plot_bkgShape_srAdjacentSlices.py, BDT_ma_2D_lib.py
#   table_MVAScore_sigEff_bkgYield.py, table_interpolate_bkgYield_{1,2}.py
# Backup first: backup_prod/anplots_before.tar (Plot/plots dirs + Plot/output tables/json)
# Logs: logs/anplots/<script>.log; status logs/anplots_status.txt; marker logs/anplots.DONE
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
L=$D/logs/anplots; mkdir -p $L
rm -f $D/logs/anplots.DONE $D/logs/anplots_status.txt
[ -s $D/backup_prod/anplots_before.tar ] || tar -C $P -cf $D/backup_prod/anplots_before.tar \
  plots/scan_score_R plots/MVASigEffVmA plots/mAmigratedBar plots/mAmigratedHist plots/mAmigratedMatrix plots/mA1_doublepeak \
  plots/signalReweight_bdtwp plots/bkgmcScupltingCheck plots/BDT_ma_2D $(cd $P && ls output/*.txt output/*.json) || { echo "backup failed"; exit 1; }
md5sum $P/output/MVAcut_points_run3.json > $D/logs/anplots_wp.md5
cd $P
one() {   # name, then the script and its arguments; extra env through ENVX
  local n=$1; shift
  env -u PYTHONHOME PYTHONPATH=$P/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts ${ENVX:-} $PY "$@" > $L/$n.log 2>&1
  echo "$n exit $?" >> $D/logs/anplots_status.txt
}
one scan_score_R scripts/scan_score_R_significance.py &
one sculptR scripts/make_sculpt_R_table_from_wp.py &
ENVX=HZA_SIGNAL_REWEIGHT=1 one MVASigEffVmA scripts/plot_MVASigEffVmA.py &
wait
ENVX=HZA_SIGNAL_REWEIGHT=1 one mAmigratedBar scripts/plot_mAmigratedBar.py &
ENVX=HZA_SIGNAL_REWEIGHT=1 one mAmigratedHist scripts/plot_mAmigratedHist.py &
ENVX=HZA_SIGNAL_REWEIGHT=1 one mAmigratedMatrix scripts/plot_mAmigratedMatrix.py &
wait
one mA1_doublepeak scripts/plot_mA1_doublepeak.py &
one sigrw_trainvars scripts/plot_signal_reweight_trainvars_bdtwp.py --ma 1,2,3 &
one sigrw_mllgg scripts/plot_signal_reweight_mllgg_bdtwp.py --ma 1,2,3 &
wait
one bkgmcSculpting scripts/plot_bkgmcSculptingCheck.py &
one srAdjacent scripts/plot_bkgShape_srAdjacentSlices.py &
one BDT_ma_2D scripts/BDT_ma_2D_lib.py &
wait
one table_MVAScore scripts/table_MVAScore_sigEff_bkgYield.py &
( one table_interp_1 scripts/table_interpolate_bkgYield_1.py; one table_interp_2 scripts/table_interpolate_bkgYield_2.py ) &
wait
md5sum -c $D/logs/anplots_wp.md5 > /dev/null || echo "WP JSON CHANGED" >> $D/logs/anplots_status.txt
date '+%F %T' > $D/logs/anplots.DONE
