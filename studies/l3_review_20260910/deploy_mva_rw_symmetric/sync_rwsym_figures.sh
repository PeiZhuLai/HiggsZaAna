#!/usr/bin/env bash
# Targeted AN figure sync after the m_a = 20 signal-fit fix (deploy_mA20_sigfit_powell) and the symmetric
# CMS_hza_mva_reweight deployment (deploy_mva_rw_symmetric); copy of deploy_mA1_low100bern4/sync_mA1_figures.sh
# restricted to the regenerated products: electron signal model (effSigmaVmA, smodel_20_*), bias plots (new axis
# label), impacts, limit plots and the pseudo-data closure summary. Background fits are not touched.
set -euo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_mva_rw_symmetric
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
R=$AN/figure_ALP/run3
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
MVA=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/plots_MVA
FILT=(--prune-empty-dirs --include='*/' --include='*.pdf' --include='*.png' --include='*.jpg' --include='*.jpeg' --include='*.tex' --include='*.txt' --include='*.json' --exclude='*')
s() { mkdir -p "$2"; echo "[sync] $1 -> $2"; rsync -a --itemize-changes "${FILT[@]}" "$1/" "$2/"; }
f() { mkdir -p "$(dirname "$2")"; rsync -a --itemize-changes "$1" "$2"; }
f $F/Signal/outdir_ele/signalFit/Plots/effSigmaVmA_ele.pdf $R/flashggFinalFit/Signal/outdir_ele/signalFit/Plots/effSigmaVmA_ele.pdf
for x in $F/Signal/outdir_ele/signalFit/Plots/smodel_20_*; do f $x $R/flashggFinalFit/Signal/outdir_ele/signalFit/Plots/$(basename $x); done
s $F/Combine/Checks/Bias_nominal/plots_bias $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
s $F/Combine/output_impacts $R/flashggFinalFit/Combine/output_impacts
s $F/Combine/output_impacts_xsec $R/flashggFinalFit/Combine/output_impacts_xsec
s $F/Plots/plot_limits $R/flashggFinalFit/Plots/plot_limits
# pseudo-data closure summary figure
PC=$F/pseudodata_closure
f $PC/plots/sculpting_closure_pseudodata.pdf $R/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots/sculpting_closure_pseudodata.pdf
f $PC/plots/sculpting_closure_pseudodata.png $R/HiggsZaAna/flashggFinalFit/pseudodata_closure/plots/sculpting_closure_pseudodata.png
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
alias_dir $R/optimize_alias $R/HiggsZaAna/Plot/plots/optimize_run3UL $R/HiggsZaAna/Plot/plots $R/HiggsZaAna/Plot/output $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
alias_dir $R/bias_alias $R/flashggFinalFit/Combine/Checks/Bias_nominal/plots_bias
alias_dir $R/impact_alias $R/flashggFinalFit/Combine/output_impacts
alias_dir $R/limit_alias $R/flashggFinalFit/Plots/plot_limits
alias_dir $R/signal_alias $R/flashggFinalFit/Signal/outdir_ele $R/flashggFinalFit/Signal/outdir_mu
echo "[done]"
