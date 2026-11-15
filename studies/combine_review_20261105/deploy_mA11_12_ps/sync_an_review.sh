#!/usr/bin/env bash
# Targeted copy of the products changed by the range deployments of the Combine review into the AN working tree:
# the subset of AN-25-172/sync_figures.sh that covers the background fits (blind plots), the envelope diagnostics
# tables, the bias plots and the q50 fixed-index impacts. Same rsync filters and alias layout as sync_figures.sh;
# nothing else of the AN is touched (the full sync_figures.sh would also pull unrelated SF/TnP products).
# Default is a dry run (rsync -n); DRY_RUN=0 copies.
# Input : flashggFinalFit/Background/ALP_BkgModel_ReReco (fit_results_run3), /eos/home-p/pelai/HZa/output_tables/bkg_envelope,
#         /eos/home-p/pelai/HZa/output_plots/plots_bias, flashggFinalFit/Combine/output_impacts{,_xsec}_q50_fixidx_prodrange
# Output: AN-25-172/figure_ALP/run3/{flashggFinalFit/...,background_alias,bias_alias}; log logs/sync_an_review.log
set -euo pipefail
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172; R=$AN/figure_ALP/run3
FF=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
DRY=${DRY_RUN:-1}; N=(); [ "$DRY" = 0 ] || N=(-n)
F=(--prune-empty-dirs --include='*/' --include='*.pdf' --include='*.png' --include='*.jpg' --include='*.jpeg'
   --include='*.tex' --include='*.txt' --include='*.json' --exclude='*')
s() { echo "[sync] $1 -> $2"; mkdir -p "$2"; rsync -a --itemize-changes "${N[@]}" "${F[@]}" "$1/" "$2/" | grep -v '/$' || true; }
s $FF/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind $R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind
s /eos/home-p/pelai/HZa/output_tables/bkg_envelope $R/flashggFinalFit/Background/output_tables
s /eos/home-p/pelai/HZa/output_plots/plots_bias $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
s $FF/Combine/output_impacts_q50_fixidx_prodrange $R/flashggFinalFit/Combine/output_impacts_q50_fixidx_prodrange
s $FF/Combine/output_impacts_xsec_q50_fixidx_prodrange $R/flashggFinalFit/Combine/output_impacts_xsec_q50_fixidx_prodrange
[ "$DRY" = 0 ] || { echo "dry run: aliases not rebuilt"; exit 0; }
# aliases as in sync_figures.sh: background_alias <- AllFitResults_blind (allPdfs_*, bkgplot_*), bias_alias <- plots_bias
B=$R/background_alias; mkdir -p $B; rm -f $B/*.pdf $B/*.png
for f in $R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind/{allPdfs,bkgplot}_*.{pdf,png}; do
  [ -f "$f" ] && cp -p "$f" $B/; done
mkdir -p $R/bias_alias
cp -p $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias/*.{pdf,png} $R/bias_alias/
echo "aliases rebuilt"
