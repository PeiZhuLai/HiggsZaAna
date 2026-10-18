#!/usr/bin/env bash
# Copy the w123 limit and effSigma figures used by the PAS from the AN figure tree (synced by sync_w123_figures.sh).
set -euo pipefail
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172/figure_ALP/run3/flashggFinalFit
PAS=/afs/cern.ch/work/p/pelai/HZa/PAS/NPS-25-014/figure_ALP/run3/flashggFinalFit
for f in Plots/plot_limits/2_runLimitsPlot/compare_wilson_run2.pdf Plots/plot_limits/2_runLimitsPlot/compare_xs_run2.pdf \
         Plots/plot_limits/2_runLimitsPlot/Limits_Wilson.pdf Plots/plot_limits/2_runLimitsPlot/Limits_XS.pdf \
         Signal/outdir_ele/signalFit/Plots/effSigmaVmA_ele.pdf Signal/outdir_mu/signalFit/Plots/effSigmaVmA_mu.pdf; do
  mkdir -p $(dirname $PAS/$f); cmp -s $AN/$f $PAS/$f && echo "same    $f" || { cp -p $AN/$f $PAS/$f; echo "updated $f"; }
done
