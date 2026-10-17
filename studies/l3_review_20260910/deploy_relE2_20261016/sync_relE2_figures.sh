#!/usr/bin/env bash
# Targeted AN figure sync for the low_swapPt_relE2 deployment (subset of AN-25-172/sync_figures.sh with the same rsync
# filters and alias rules). Only the products regenerated for this deployment are synced; the SF web summaries,
# HZgamma trigger plots and the rest of Plot/plots are not touched.
# Log: logs/sync_relE2_figures.log (rsync --itemize-changes); list of changed AN files: logs/sync_relE2_changed.txt
set -euo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_relE2_20261016
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
R=$AN/figure_ALP/run3
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
MVA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/plots_MVA
FILT=(--prune-empty-dirs --include='*/' --include='*.pdf' --include='*.png' --include='*.jpg' --include='*.jpeg' --include='*.tex' --include='*.txt' --include='*.json' --exclude='*')
s() { mkdir -p "$2"; echo "[sync] $1 -> $2"; rsync -a --itemize-changes "${FILT[@]}" "$1/" "$2/"; }
f() { mkdir -p "$(dirname "$2")"; rsync -a --itemize-changes "$1" "$2"; }
for d in variables_dataVmc/plot_UL_CR variables_dataVmc/plot_UL_mva train_vars_sigVbkg optimize_run3UL signalReweight_bdtwp \
         scan_score_R MVASigEffVmA mAmigratedBar mAmigratedHist mAmigratedMatrix mA1_doublepeak signal_eff_sumw \
         bkgmcScupltingCheck BDT_ma_2D; do
  s $P/plots/$d $R/HiggsZaAna/Plot/plots/$d
done
for x in latexTable_MVAcut_eff_bkg.txt latexTable_significanceVcat.txt latexTable_sculptR.txt latexTable_bkg_interpolated_yields.txt \
         MVAcut_points_run3.json sigEfficiencyVmA_ele_byYear_5years_quadratic_interp_ma_points.json \
         sigEfficiencyVmA_muon_byYear_5years_quadratic_interp_ma_points.json; do
  [ -f $P/output/$x ] && f $P/output/$x $R/HiggsZaAna/Plot/output/$x
done
s $MVA/run3_lowmass $R/HiggsZaAna/HZaMVA/plots_MVA/run3_lowmass
for m in 1 2 3; do s $F/Background/ALP_BkgModel_ReReco/fit_results_run3/$m $R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/$m; done
s $F/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind $R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind
s $F/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_unblind $R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_unblind
s $F/Background/output_tables $R/flashggFinalFit/Background/output_tables
s $F/Signal/outdir_ele $R/flashggFinalFit/Signal/outdir_ele
s $F/Signal/outdir_mu $R/flashggFinalFit/Signal/outdir_mu
s $F/Combine/Checks/Bias_nominal/plots_bias $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
s $F/Combine/output_impacts $R/flashggFinalFit/Combine/output_impacts
s $F/Combine/output_impacts_xsec $R/flashggFinalFit/Combine/output_impacts_xsec
s $F/Plots/plot_limits $R/flashggFinalFit/Plots/plot_limits
# pseudo-data closure: summary figure and the unblind background-fit panels of m_a 1-3
PC=$F/pseudodata_closure
f $PC/plots/sculpting_closure_pseudodata.pdf $R/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots/sculpting_closure_pseudodata.pdf
f $PC/plots/sculpting_closure_pseudodata.png $R/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots/sculpting_closure_pseudodata.png
for m in 1 2 3; do for k in bkgplot allPdfs; do
  x=$PC/fit_results_run3_pseudodata/AllFitResults_unblind/${k}_$m.pdf; [ -f $x ] && f $x $R/pseudodata_bkg_alias/${k}_$m.pdf
done; done
# aliases (flat copies, same rules as sync_figures.sh)
alias_dir() {
  python3 - "$@" <<'PY'
from pathlib import Path
import shutil, sys
alias_dir = Path(sys.argv[1]).resolve(); sources = [Path(p).resolve() for p in sys.argv[2:]]
exts = {".pdf", ".png", ".jpg", ".jpeg", ".tex", ".txt", ".json"}
alias_dir.mkdir(parents=True, exist_ok=True); seen = {}
for src in sources:
    if src.exists():
        for path in sorted(src.rglob("*")):
            if path.is_file() and path.suffix.lower() in exts:
                seen.setdefault(path.name, path)
n = 0
for name, src in sorted(seen.items()):
    dst = alias_dir / name
    if not dst.exists() or dst.read_bytes() != src.read_bytes():
        shutil.copy2(src, dst); n += 1; print(f"[alias] {src} -> {dst}")
print(f"[alias] {alias_dir.name}: {n} files updated of {len(seen)}")
PY
}
alias_dir $R/selection_alias $R/HiggsZaAna/Plot/plots $R/HiggsZaAna/Plot/output
alias_dir $R/optimize_alias $R/HiggsZaAna/Plot/plots/optimize_run3UL $R/HiggsZaAna/Plot/plots $R/HiggsZaAna/Plot/output $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
alias_dir $R/bias_alias $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
alias_dir $R/impact_alias $R/flashggFinalFit/Combine/output_impacts
alias_dir $R/limit_alias $R/flashggFinalFit/Plots/plot_limits
alias_dir $R/signal_alias $R/flashggFinalFit/Signal/outdir_ele $R/flashggFinalFit/Signal/outdir_mu
for m in 1 2 3; do for k in allPdfs bkgplot; do for e in pdf png; do
  x=$R/flashggFinalFit/Background/ALP_BkgModel_ReReco/fit_results_run3/AllFitResults_blind/${k}_$m.$e; [ -f $x ] && f $x $R/background_alias/${k}_$m.$e
done; done; done
echo "[done]"
