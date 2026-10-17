#!/usr/bin/env bash
# mA1 closure bootstrap (MC-statistics spread of the closure Z / r_hat), pseudo-data without the BDT training DY:
# new low-mass model at the cuts of 2x and 1x the production data-sideband count, production model at its working
# point. 30 replicas each, 3 parallel, production mA1 datacard. Starts after logs/scans.DONE. Marker: logs/boot.DONE
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
DC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt
cd $V; rm -f logs/boot.DONE
until [ -e logs/scans.DONE ]; do sleep 30; done
bash boot_closure.sh prod 1 0.963 30 $DC 3 > logs/boot_prod_drop_mA1.log 2>&1
bash boot_closure.sh new 1 0.97804 30 $DC 3 > logs/boot_new_drop_mA1_097804.log 2>&1
bash boot_closure.sh new 1 0.96443 30 $DC 3 > logs/boot_new_drop_mA1_096443.log 2>&1
cat /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_boot/*/summary.txt > logs/boot.DONE
