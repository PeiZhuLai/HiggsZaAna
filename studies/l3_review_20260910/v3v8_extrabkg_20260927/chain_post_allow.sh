#!/usr/bin/env bash
#
# Post-production chain for the V3/V8 extra-background samples (2024):
#   [0] wait for the HiggsDNA driver (logs/driver_<MODE>.done) -- the driver itself
#       monitors the condor jobs, resubmits failures and merges (--merge_outputs)
#   [1] completeness: job dirs vs *_summary_job<N>.json (the completion marker at fpo=1;
#       zero selected events writes a summary but no parquet)
#   [2] three-stage check, stage 1+2: chunk rows (find, excluding merged_nominal.parquet)
#       == merged rows, and merged mtime >= newest chunk mtime (merge-before-rechunk race)
#   [3] p2root (--split, like run3_bdt_inputs_fsrfix) and BDT scoring (no split, like the
#       Bkg lines of the run3_bdt_scored_fsrfix joblist), same converter + same models
#   [4] three-stage check, stage 3: ROOT inclusive entries == merged rows
#   [5] yields at the 30 working points (compute_extrabkg_yields.py)
# Heartbeat: $D/logs/heartbeat_<MODE>.txt (AFS close-to-open: look at it, not at the log).
#
# Usage: MODE=test|full [ALLOW_NOTRUN=sample/job_N,...] bash chain_post_allow.sh
set -uo pipefail
MODE="${MODE:-full}"
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v3v8_extrabkg_20260927
PBASE=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_extraBkg
RBASE=/eos/home-p/pelai/HZa/root_P2Root
case "$MODE" in
  test) IN=$PBASE/Bkg_MC_extraBkgTEST
        OUT_IN=$RBASE/run3_bdt_inputs_fsrfix_extraBkgTEST
        OUT_SC=$RBASE/run3_bdt_scored_fsrfix_extraBkgTEST ;;
  full) IN=$PBASE/Bkg_MC_extraBkg2024
        OUT_IN=$RBASE/run3_bdt_inputs_fsrfix_extraBkg
        OUT_SC=$RBASE/run3_bdt_scored_fsrfix_extraBkg ;;
esac
SAMPLES="TTto2L2Nu TTG TTGG_Run3 ZGG"
HB=$D/logs/heartbeat_${MODE}.txt
PY=/eos/home-p/pelai/App/Conda/.conda/envs/hza_ana/bin/python
hb() { echo "$(date '+%F %T') $*" >> "$HB"; }

hb "chain_post start MODE=$MODE"
# [0] wait for driver (checked every 5 min; the heartbeat records condor state each time)
while [ ! -f "$D/logs/driver_${MODE}.done" ]; do
  q=$(condor_q -constraint "regexp(\"eos_logs/$(basename $IN)/\", Cmd)" -af JobStatus 2>/dev/null | sort | uniq -c | tr '\n' ' ')
  hb "[0] waiting for driver; condor JobStatus counts: ${q:-none}"
  sleep 300
done
hb "[0] driver done: $(cat $D/logs/driver_${MODE}.done)"

# [1] + [2]
$PY - "$IN" $SAMPLES > "$D/logs/reconcile_${MODE}.txt" 2>&1 <<'PYEOF'
import os, sys, glob
import pyarrow.parquet as pq
base = sys.argv[1]; rc = 0
for s in sys.argv[2:]:
    d = os.path.join(base, s + "_2024")
    jobs = sorted(x for x in os.listdir(d) if x.startswith("job_"))
    notrun = [j for j in jobs if not glob.glob(os.path.join(d, j, "*_summary_%s.json" % j.replace("job_", "job")))]
    # 2026-09-29: explicitly accepted missing jobs (ALLOW_NOTRUN="sample/job_N,..."). TTGG_Run3 job_23
    # retired after repeated xrdcp failures through the INFN redirector; the merge normalizes with the
    # sum of weights of the completed jobs only, so the remaining 61/62 jobs stay correctly scaled.
    allow = set(x for x in os.environ.get("ALLOW_NOTRUN", "").split(",") if x)
    accepted = [j for j in notrun if "%s/%s" % (s, j) in allow]
    notrun = [j for j in notrun if j not in accepted]
    if accepted: print("%-16s accepted missing jobs: %s" % (s, ",".join(accepted)))
    # find, not a shell glob; exclude merged_nominal.parquet itself (else ~2x)
    chunks = []
    for root, _, files in os.walk(d):
        for f in files:
            if f.endswith("_nominal.parquet") and f != "merged_nominal.parquet":
                chunks.append(os.path.join(root, f))
    crow = sum(pq.ParquetFile(f).metadata.num_rows for f in chunks)
    mp = os.path.join(d, "merged_nominal.parquet")
    # no chunk at all (every job selected zero events) legitimately gives no merged file
    mrow = pq.ParquetFile(mp).metadata.num_rows if os.path.exists(mp) else (0 if not chunks else -1)
    late = [f for f in chunks if os.path.exists(mp) and os.path.getmtime(f) > os.path.getmtime(mp)]
    ok = (not notrun) and crow == mrow and not late
    rc |= (not ok)
    print("%-16s jobs=%4d ran=%4d chunks=%4d chunk_rows=%8d merged_rows=%8d late_chunks=%d %s"
          % (s, len(jobs), len(jobs) - len(notrun), len(chunks), crow, mrow, len(late), "OK" if ok else "FAIL"))
sys.exit(rc)
PYEOF
rc12=$?
hb "[1-2] completeness + chunk==merged rc=$rc12 -> logs/reconcile_${MODE}.txt"
cat "$D/logs/reconcile_${MODE}.txt" >> "$HB"
[ $rc12 -ne 0 ] && { hb "STOP: reconcile failed"; echo "CHAIN_DONE rc=12" > $D/logs/chain_${MODE}.done; exit 12; }

# [3] p2root + scoring, 1 local process at a time (lxplus foreground budget)
# --sideband-reweight-mode always: the auto rule keys on a path part named exactly
# "Bkg_MC", which these dirs deliberately do not have (eos_logs collision). DY got the
# branches through that rule; 'factor' (used for the yields) is unaffected either way.
for s in $SAMPLES; do
  src=$IN/${s}_2024/merged_nominal.parquet
  [ -f "$src" ] || { hb "[3] $s: no merged parquet (zero selected events) -- skipped"; continue; }
  bash $D/score_local.sh "$src" "$OUT_IN/$s/2024.root" --split --sideband-reweight-mode always > $D/logs/p2root_${MODE}_$s.log 2>&1
  hb "[3] p2root(--split) $s rc=$?"
  bash $D/score_local.sh "$src" "$OUT_SC/$s/2024.root" --sideband-reweight-mode always > $D/logs/score_${MODE}_$s.log 2>&1
  hb "[3] scoring $s rc=$?"
done

# [4] merged rows == ROOT inclusive entries, scores not all-NaN, ROOT newer than merged
$PY - "$IN" "$OUT_IN" "$OUT_SC" $SAMPLES > "$D/logs/stage3_${MODE}.txt" 2>&1 <<'PYEOF'
import os, sys
import numpy as np, uproot, pyarrow.parquet as pq
IN, OI, OS = sys.argv[1:4]; rc = 0
for s in sys.argv[4:]:
    mp = os.path.join(IN, s + "_2024", "merged_nominal.parquet")
    if not os.path.exists(mp):
        print("%-10s no merged parquet (zero selected) -- nothing to convert" % s); continue
    mrow = pq.ParquetFile(mp).metadata.num_rows
    for lab, base in (("inputs", OI), ("scored", OS)):
        rp = os.path.join(base, s, "2024.root")
        if not os.path.exists(rp):
            print("%-10s %-7s MISSING %s" % (s, lab, rp)); rc = 1; continue
        t = uproot.open(rp)["inclusive"]; n = t.num_entries
        fresh = os.path.getmtime(rp) >= os.path.getmtime(mp)
        sc = t["MVA_Score_mA_M5"].array(library="np") if n else np.array([])
        nan_frac = float(np.mean(~np.isfinite(sc))) if n else 0.0
        ok = (n == mrow) and fresh and (n == 0 or nan_frac < 1.0)
        rc |= (not ok)
        print("%-10s %-7s merged=%8d root_inclusive=%8d fresh=%s nanfrac(M5)=%.3f %s"
              % (s, lab, mrow, n, fresh, nan_frac, "OK" if ok else "FAIL"))
sys.exit(rc)
PYEOF
rc4=$?
hb "[4] stage-3 root==merged rc=$rc4 -> logs/stage3_${MODE}.txt"
cat "$D/logs/stage3_${MODE}.txt" >> "$HB"

# [5] yields
if [ "$MODE" = "full" ]; then
  $PY $D/compute_extrabkg_yields.py --new-base "$OUT_SC" --out-prefix $D/yields_extrabkg > $D/logs/yields_full.log 2>&1
else
  $PY $D/compute_extrabkg_yields.py --new-base "$OUT_SC" --out-prefix $D/logs/yields_TEST > $D/logs/yields_test.log 2>&1
fi
rc5=$?
hb "[5] yields rc=$rc5"
echo "CHAIN_DONE rc12=$rc12 rc4=$rc4 rc5=$rc5" > $D/logs/chain_${MODE}.done
hb "chain_post done"
