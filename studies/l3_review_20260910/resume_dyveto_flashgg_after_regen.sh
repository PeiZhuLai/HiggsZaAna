#!/usr/bin/env bash
# 2026-09-30 01:5x: MVA-cut data for mA2 (output_2023preBPix, segfault on read) and mA28
# (output_2023postBPix, read fails at entry 1/49) are corrupt on EOS (persistent on retry).
# 1. regenerate mA2 + mA28 data (regen_mvacut_data.sh: scratch -> full read -> EOS -> full read)
# 2. per-file full-read audit of all 250 MVA-cut files (one process per file: a segfault names the file)
# 3. flashgg from Tree2WS (START_STAGE=4), then dataVmc nominal
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
HB=$L/dyveto_driver_heartbeat.txt
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }
hb "regen MVA-cut data mA2 mA28"
bash $D/regen_mvacut_data.sh 2 28 > $L/regen_mvacut_2_28.log 2>&1 || { hb "STOPPING: regen failed (regen_mvacut_2_28.log)"; exit 1; }
bash $D/audit_mvacut_perfile.sh > $L/mva_audit_perfile_after_regen.log 2>&1
grep -q "bad 0$" $L/mva_audit_perfile_after_regen.log || { hb "STOPPING: audit after regen: $(head -1 $L/mva_audit_perfile_after_regen.log)"; exit 1; }
hb "audit after regen: $(head -1 $L/mva_audit_perfile_after_regen.log)"
hb "flashgg resume at stage 4"
START_STAGE=4 bash $D/chain_dyveto_flashgg.sh > $L/dyveto_flashgg_s4.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_flashgg_s4.log || { hb "STOPPING: flashgg stage>=4 (dyveto_flashgg_s4.log)"; exit 1; }
hb "flashgg passed"
echo "$(date '+%F %H:%M') DONE dyveto-flashgg (logs dyveto_flashgg.log + dyveto_flashgg_s4.log). Next: mA3 closure WP scan, impacts/bias, closure/plots, AN." >> $L/PENDING_RESULTS.txt
hb "dataVmc nominal start"
FT=nominal bash $D/chain_dyveto_datavmc.sh > $L/dyveto_datavmc_nominal.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_datavmc_nominal.log || { hb "STOPPING: dataVmc nominal"; exit 1; }
hb "dataVmc nominal passed"; hb "DRIVER DONE"
