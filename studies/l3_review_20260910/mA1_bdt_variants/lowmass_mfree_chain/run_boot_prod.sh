#!/usr/bin/env bash
# production-model closure bootstrap, 30 replicas per mass, 3 parallel. Marker: logs/boot_prod.DONE
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain; rm -f logs/boot_prod.DONE
bash boot_closure.sh prod 1 0.963 30 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt 3 > logs/boot_prod_mA1.log 2>&1
bash boot_closure.sh prod 2 0.975 30 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/2_pruned_datacard_leptons.txt 3 > logs/boot_prod_mA2.log 2>&1
bash boot_closure.sh prod 3 0.988 30 /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/3_pruned_datacard_leptons.txt 3 > logs/boot_prod_mA3.log 2>&1
cat /eos/home-p/pelai/HZa/sandbox_lowmass_mfree/closure_boot/prod_mA*/summary.txt > logs/boot_prod.DONE
