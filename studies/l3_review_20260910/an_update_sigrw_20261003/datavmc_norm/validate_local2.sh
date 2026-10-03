#!/bin/bash
# Second local validation (after the --optimizeBranches fix): mA5, era ${ERA:-2022preEE}, SR.
#   new       : current 1_prepare_dataVmc.py (R with all reweight branches enabled, times N), with --optimizeBranches
#   oldnoopt  : pre-change backup WITHOUT --optimizeBranches (correct R, no N)
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003/datavmc_norm
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
J=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json
E=${ERA:-2022preEE}
cd $P
export PYTHONPATH="${PYTHONPATH:-}:$P/lib:/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts"
export PYTHONUNBUFFERED=1
common=(-y run3 -m --ln --histOnly --samples "M5@$E" --region 1 --useSidebandReweight --sidebandReweightJson $J --noSidebandReweightUnc)
$PY -u scripts/1_prepare_dataVmc.py "${common[@]}" --backend auto --optimizeBranches --outputTag sigrwnorm_val_new_$E > $W/val_new_$E.log 2>&1 &
p1=$!
$PY -u scripts/1_prepare_dataVmc.py.bak_sigrwnorm_20261003 "${common[@]}" --backend loop --outputTag sigrwnorm_val_oldnoopt_$E > $W/val_oldnoopt_$E.log 2>&1 &
p2=$!
wait $p1; echo "new rc=$?"
wait $p2; echo "oldnoopt rc=$?"
for t in new oldnoopt; do
  mv $P/plots/variables_dataVmc/ALP_plot_run3_UL_SR_sigrwnorm_val_${t}_$E.root $W/
  rmdir $P/plots/variables_dataVmc/plot_UL_SR/sigrwnorm_val_${t}_$E 2>/dev/null
done
ls -la $W/*_$E.root
