#!/usr/bin/env bash
# Low-mass model low_nodR (variable-search low model without min_lg_deltaR): closure + data-GOF scans for mA1-3 at
# the cuts of eval/closure_cuts_nodr.txt (production data-sideband count x 2, 1.5, 1, 0.7, 0.5), then the mA1
# bootstrap at x1 and x0.5 (30 replicas each). Pseudo-data without the BDT training DY (DROP=1). At most 6
# processes. Marker: logs/nodr.DONE
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate
DC=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Datacard/output_Datacard_leptons/1_pruned_datacard_leptons.txt
B=/eos/home-p/pelai/HZa/sandbox_varsearch_v1
cd $V; rm -f logs/nodr.DONE
cuts() { awk -v m=$1 '$1==m {$1=""; print}' eval/closure_cuts_nodr.txt; }
cut_n() { awk -v m=$1 -v i=$2 '$1==m {print $(i+1)}' eval/closure_cuts_nodr.txt; }
scan() { DROP=1 NPAR=3 bash scan_closure.sh nodr $1 "$(cuts $1)" > logs/scan_nodr_drop_mA$1.log 2>&1; }
( scan 1; scan 3 ) & ( scan 2 ) & wait
bash boot_closure.sh nodr 1 $(cut_n 1 3) 30 $DC 3 > logs/boot_nodr_drop_mA1_x1.log 2>&1 &
bash boot_closure.sh nodr 1 $(cut_n 1 5) 30 $DC 3 > logs/boot_nodr_drop_mA1_x05.log 2>&1 &
wait
cat $B/closure_scan/nodr_drop_mA*/summary.txt $B/closure_boot/nodr_drop_mA1_cut*/summary.txt > logs/nodr.DONE
