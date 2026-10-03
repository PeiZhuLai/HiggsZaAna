#!/bin/bash
# Local validation of the signal normalization N in 1_prepare_dataVmc.py:
# mA5, era 2022preEE, SR region, new script vs the pre-change backup (.bak_sigrwnorm_20261003).
# Outputs are moved from Plot/plots/variables_dataVmc into this directory afterwards.
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003/datavmc_norm
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
J=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json
cd $P
export PYTHONPATH="${PYTHONPATH:-}:$P/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts"
export PYTHONUNBUFFERED=1
common=(-y run3 -m --ln --histOnly --backend loop --optimizeBranches --samples "M5@${ERA:-2022preEE}" --region 1 --useSidebandReweight --sidebandReweightJson $J --noSidebandReweightUnc)
$PY -u scripts/1_prepare_dataVmc.py "${common[@]}" --outputTag sigrwnorm_val_new > $W/val_new.log 2>&1 &
p1=$!
$PY -u scripts/1_prepare_dataVmc.py.bak_sigrwnorm_20261003 "${common[@]}" --outputTag sigrwnorm_val_old > $W/val_old.log 2>&1 &
p2=$!
wait $p1; echo "new rc=$?"
wait $p2; echo "old rc=$?"
mv $P/plots/variables_dataVmc/ALP_plot_run3_UL_SR_sigrwnorm_val_new.root $P/plots/variables_dataVmc/ALP_plot_run3_UL_SR_sigrwnorm_val_old.root $W/
rmdir $P/plots/variables_dataVmc/plot_UL_SR/sigrwnorm_val_new $P/plots/variables_dataVmc/plot_UL_SR/sigrwnorm_val_old 2>/dev/null
ls -la $W/*.root
