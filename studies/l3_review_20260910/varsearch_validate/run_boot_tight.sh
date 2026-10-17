#!/usr/bin/env bash
# mA1 closure bootstrap of the new low-mass model at the tight cuts (0.7x and 0.5x the production data-sideband
# count), 30 replicas each, both cuts at once with 3 parallel each (6 processes). Marker: logs/boot_tight.DONE
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
cd $V; rm -f logs/boot_tight.DONE
bash boot_closure.sh new 1 0.98207 30 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt 3 > logs/boot_new_drop_mA1_098207.log 2>&1 &
bash boot_closure.sh new 1 0.98566 30 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt 3 > logs/boot_new_drop_mA1_098566.log 2>&1 &
wait
cat /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_boot/new_drop_mA1_cut0982*/summary.txt /eos/home-p/pelai/HZa/sandbox_varsearch_v1/closure_boot/new_drop_mA1_cut0985*/summary.txt > logs/boot_tight.DONE
