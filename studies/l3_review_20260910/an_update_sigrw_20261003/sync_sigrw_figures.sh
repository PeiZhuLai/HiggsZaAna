#!/bin/bash
# Targeted figure sync for the signal-reweight AN update (subset of AN-25-172/sync_figures.sh;
# bias plots are deliberately NOT synced). Same rsync filters as sync_figures.sh.
set -euo pipefail
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
R=$AN/figure_ALP/run3
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
FILT=(--prune-empty-dirs --include='*/' --include='*.pdf' --include='*.png' --include='*.jpg' --include='*.jpeg' --include='*.tex' --include='*.txt' --include='*.json' --exclude='*')
s() { mkdir -p "$2"; echo "[sync] $1 -> $2"; rsync -a --itemize-changes "${FILT[@]}" "$1/" "$2/"; }
s $F/Plots/plot_limits/2_runLimitsPlot   $R/flashggFinalFit/Plots/plot_limits/2_runLimitsPlot
s $F/Combine/output_impacts              $R/flashggFinalFit/Combine/output_impacts
s $F/Signal/outdir_ele                   $R/flashggFinalFit/Signal/outdir_ele
s $F/Signal/outdir_mu                    $R/flashggFinalFit/Signal/outdir_mu
for d in MVASigEffVmA signal_eff_sumw mA1_doublepeak mAmigratedBar mAmigratedMatrix mAmigratedHist; do
  s $P/plots/$d $R/HiggsZaAna/Plot/plots/$d
done
for f in MVASigEffVmA_ele_byYear_5years_quadratic_interp_ma_points.json MVASigEffVmA_muon_byYear_5years_quadratic_interp_ma_points.json \
         sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json; do
  rsync -a --itemize-changes $P/output/$f $R/HiggsZaAna/Plot/output/$f
done
# signal_alias: smodel_<m>_run3_<ch> used by Sec-08-01 (flat copies, as make_alias_dir does)
for ch in ele mu; do
  for f in $R/flashggFinalFit/Signal/outdir_$ch/signalFit/Plots/smodel_*_run3_$ch.*; do
    cp -p "$f" $R/signal_alias/$(basename "$f")
  done
done
echo "[done]"
