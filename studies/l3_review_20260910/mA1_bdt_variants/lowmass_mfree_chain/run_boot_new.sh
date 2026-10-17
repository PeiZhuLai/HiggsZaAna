#!/usr/bin/env bash
# win_mfree closure bootstrap with the sandbox datacards, 30 replicas per mass, 3 parallel. Marker: logs/boot_new.DONE
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain; rm -f logs/boot_new.DONE
bash boot_closure.sh new 1 0.9645 30 /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt 3 > logs/boot_new_mA1.log 2>&1
bash boot_closure.sh new 2 0.975 30 /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain/flashggFinalFit/Datacard/output_Datacard_leptons/2_pruned_datacard_leptons.txt 3 > logs/boot_new_mA2.log 2>&1
bash boot_closure.sh new 3 0.985 30 /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain/flashggFinalFit/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt 3 > logs/boot_new_mA3.log 2>&1
cat /eos/home-p/pelai/HZa/sandbox_lowmass_mfree/closure_boot/new_mA*/summary.txt > logs/boot_new.DONE
