#!/usr/bin/env bash
# 2026-09-30 05:3x: mva_reweight (S1 scheme 2) and photon_csev were typed a_h -> symmetrized kappas.
# calcSystematics.py now forces a_w (.bak_awforce_20260930). Wait for the running flashgg
# (hza-dyveto-flashgg-resume) to pass, keep its limits as the symmetrized reference, then rerun
# stages 7-9 (datacards -> limits -> plots).
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
F=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit
HB=$L/dyveto_driver_heartbeat.txt
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }
until grep -q "^VERDICT: PASSED\|^STOPPING" $L/dyveto_flashgg_s4.log 2>/dev/null; do sleep 120; done
grep -q "^VERDICT: PASSED" $L/dyveto_flashgg_s4.log || { hb "redo-aw: flashgg s4 did not pass, not redoing"; exit 1; }
S=$D/snapshot_dyveto_symS1lnN_20260930; mkdir -p $S
cp -r $F/Combine/output_combine_results $S/ ; cp -r $F/Datacard/output_Datacard_leptons $S/
hb "redo-aw: symmetrized-lnN limits kept in $S; rerunning stages 7-9 with a_w"
START_STAGE=7 bash $D/chain_dyveto_flashgg.sh > $L/dyveto_flashgg_s7_aw.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_flashgg_s7_aw.log && hb "redo-aw: stages 7-9 PASSED" || { hb "STOPPING: redo-aw stages 7-9 (dyveto_flashgg_s7_aw.log)"; exit 1; }
grep -h "CMS_hza_mva_reweight\|CMS_hza_photon_csev_2023postBPix" $F/Datacard/output_Datacard_leptons/{1,5,20}_pruned_datacard_leptons.txt | awk '{print $1,$2,$3,$4}' >> $L/dyveto_flashgg_s7_aw.log
echo "$(date '+%F %H:%M') DONE redo datacards/limits with a_w for mva_reweight + csev (log dyveto_flashgg_s7_aw.log; symmetrized reference in snapshot_dyveto_symS1lnN_20260930)" >> $L/PENDING_RESULTS.txt
