#!/bin/bash
# Regenerate the signal-efficiency (Fig. mva_sig_eff_run3) and migration plots with the normalized signal reweight.
# Input: /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix ; Output: Plot/plots/{MVASigEffVmA,mAmigrated*}
# Logs: an_update_sigrw_20261003/run_*.log ; done marker run_plots_sigrw.DONE
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
run() { HZA_SIGNAL_REWEIGHT=1 env -u PYTHONHOME PYTHONPATH=$PWD/lib $PY scripts/$1 > $W/run_${1%.py}.log 2>&1; echo "$1 exit $?" >> $W/run_plots_sigrw.status; }
rm -f $W/run_plots_sigrw.status $W/run_plots_sigrw.DONE
run plot_MVASigEffVmA.py & run plot_mAmigratedBar.py & run plot_mAmigratedMatrix.py & run plot_mAmigratedHist.py &
wait
date > $W/run_plots_sigrw.DONE
