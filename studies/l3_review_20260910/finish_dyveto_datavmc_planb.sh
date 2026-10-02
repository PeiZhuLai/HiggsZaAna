#!/usr/bin/env bash
# Plan B for the dyveto dataVmc sideband_rwgt run (2026-09-30 00:4x).
# The 3 dyll_2024 jobs (12778398.5/.30/.55) ran >3.6 h in the Python loop because the scored
# DY 2024 files had no weight_sideband_rwgt (p2root matched only an exact "Bkg_MC" dir; fixed).
# DY 2024 was rescored with the fixed p2root into run3_bdt_scored_fsrfix_dyrescore_tmp and verified:
# all 280 common branches identical, entries equal, only the sideband branches are new.
#   1. stop the driver, remove the 3 slow jobs (being replaced, not held)
#   2. swap the rescored files into run3_bdt_scored_fsrfix (old kept as *.nobranch_20260930), re-hadd DYJetsToLL/2024
#   3. rerun the 3 dataVmc jobs (fast path), check all 75 jobs, merge + draw (same checks as chain_dyveto_datavmc.sh)
#   4. continue the driver's remaining steps: flashgg, then dataVmc nominal
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
L=$D/logs_fsrfix
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/Condor
V=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/variables_dataVmc
LS=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/logs_split
SC=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix
TMP=/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix_dyrescore_tmp
LOG=$L/dyveto_datavmc_sideband_rwgt.log
HB=$L/dyveto_driver_heartbeat.txt
FT=sideband_rwgt
hb(){ echo "$(date '+%F %T') $*" | tee -a "$HB"; }
die(){ hb "STOPPING (planB): $*"; echo "VERDICT: FAILED (planB: $*)" >> "$LOG"; exit 1; }

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana
set -u

hb "planB: stop driver, replace the 3 slow dyll_2024 dataVmc jobs"
systemctl --user stop hza-dyveto-driver 2>/dev/null
condor_rm 12778398.5 12778398.30 12778398.55 >/dev/null 2>&1
sleep 30

for s in DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau; do
  [ -s $SC/$s/2024.root.nobranch_20260930 ] || mv $SC/$s/2024.root $SC/$s/2024.root.nobranch_20260930 || die "mv $s"
  cp -f $TMP/$s/2024.root $SC/$s/2024.root || die "cp $s"
done
[ -s $SC/DYJetsToLL/2024.root.nobranch_20260930 ] || mv $SC/DYJetsToLL/2024.root $SC/DYJetsToLL/2024.root.nobranch_20260930
( cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Parquet2Rootfile/Condor && bash 4_prepaare_2024DYJetsToLL.sh ) > $L/planb_hadd_dyll.log 2>&1
python - <<PY || die "rebuilt DYJetsToLL/2024 check"
import uproot, sys
t = uproot.open("$SC/DYJetsToLL/2024.root")["inclusive"]
ok = t.num_entries == 188313 and "weight_sideband_rwgt" in t.keys()
print("DYJetsToLL/2024 entries", t.num_entries, "has weight_sideband_rwgt", "weight_sideband_rwgt" in t.keys())
sys.exit(0 if ok else 1)
PY
hb "planB: DY 2024 swapped (entries 188313, sideband branches present)"

T1=$(date -d "2026-09-29 20:47:00" +%s)
cd $P
sed -n '6p;31p;56p' dataVmc_jobs.txt > dataVmc_planb_jobs.txt
grep -c "dyll_2024" dataVmc_planb_jobs.txt | grep -q '^3$' || die "planB jobs file"
sed -e "s|$P/dataVmc_jobs.txt|$P/dataVmc_planb_jobs.txt|" -e 's|"HZa_dataVmc"|"HZa_dataVmc_planb"|' dataVmc.submit > dataVmc_planb.submit
out=$(condor_submit dataVmc_planb.submit 2>&1); CL=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
[ -n "$CL" ] || die "submit: $out"
hb "planB: 3 dyll_2024 jobs resubmitted, cluster $CL"
while [ "$(condor_q $CL -format '%d\n' ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  [ "$(condor_q $CL -constraint 'JobStatus==5' -format '%d\n' ProcId 2>/dev/null | wc -l)" -gt 0 ] && hb "planB: held: $(condor_q $CL -constraint 'JobStatus==5' -af HoldReason | head -1 | cut -c1-120)"
  sleep 120
done

check_jobs() {
  local jobs="$1" since="$2" bad=0 n=0 region tag stag samples suffix part log
  while read -r region tag stag samples; do
    [ -z "${region:-}" ] && continue
    n=$((n+1))
    case "$region" in SR) suffix=UL_SR;; CR) suffix=UL_CR;; mva) suffix=UL_mva;; esac
    part="$V/ALP_plot_run3_${suffix}_${tag}_part_${stag}.root"
    log="$LS/${tag}_${region}_${stag}.log"
    if [ ! -s "$part" ] || [ "$(stat -c %Y "$part")" -lt "$since" ]; then echo "  STALE/MISSING partial: $part"; bad=$((bad+1)); continue; fi
    if [ ! -s "$log" ] || [ "$(stat -c %Y "$log")" -lt "$since" ]; then echo "  STALE/MISSING log: $log"; bad=$((bad+1)); continue; fi
    grep -q "\[FAILED\]" "$log" && { echo "  FAILED in log: $log"; bad=$((bad+1)); continue; }
    grep -q "Fetch conda tarball" "$log" || { echo "  NOT LOCALIZED: $log"; bad=$((bad+1)); continue; }
    grep -q "Activate conda env" "$log" && { echo "  ACTIVATED FROM EOS: $log"; bad=$((bad+1)); continue; }
    grep -q "\[DONE\]" "$log" || { echo "  NO [DONE]: $log"; bad=$((bad+1)); continue; }
  done < "$jobs"
  echo "  checked $n jobs, problems $bad"
  [ "$n" -gt 0 ] && [ "$bad" -eq 0 ]
}
{
echo "=============== [planB 5/6] every job: fresh partial, clean log, localized env ==============="
check_jobs dataVmc_jobs.txt "$T1"
} >> $LOG 2>&1 || die "job check (see $LOG)"
grep -q "falling back" $LS/${FT}_SR_dyll_2024.log && hb "planB: WARNING dyll_2024 still used the loop"

echo "=============== [planB 6/6] merge + draw ===============" >> $LOG
T2=$(date +%s)
export DATA_VMC_MAKE_JOBS_ARGS="--final-tags $FT"
bash 3_merge_dataVmc_condor.sh >> $LOG 2>&1
echo "raw exit: $? (not a verdict)" >> $LOG
ok=1
for r in SR CR mva; do
  f="$V/ALP_plot_run3_UL_${r}_${FT}.root"
  if [ -s "$f" ] && [ "$(stat -c %Y "$f")" -ge "$T2" ]; then echo "  merged OK  $f" >> $LOG; else echo "  merged STALE/MISSING  $f" >> $LOG; ok=0; fi
done
for r in SR CR mva; do
  n=$(find "$V/plot_UL_${r}/${FT}" -name '*.pdf' -newermt "@$T2" 2>/dev/null | wc -l)
  echo "  fresh pdfs in plot_UL_${r}/${FT}: $n" >> $LOG
  [ "$n" -eq 0 ] && ok=0
done
echo "  merge+draw took $(( $(date +%s) - T2 )) s" >> $LOG
[ $ok -eq 1 ] || die "merge/draw (see $LOG)"
echo "VERDICT: PASSED" >> $LOG
hb "dataVmc sideband_rwgt passed (planB; merge+draw $(( $(date +%s) - T2 )) s)"

hb "flashgg start"
bash $D/chain_dyveto_flashgg.sh > $L/dyveto_flashgg.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_flashgg.log || { hb "STOPPING: flashgg (dyveto_flashgg.log)"; exit 1; }
hb "flashgg passed"
echo "$(date '+%F %H:%M') DONE dyveto-flashgg (log dyveto_flashgg.log). Next: mA3 closure WP scan, impacts/bias, closure/plots, AN." >> $L/PENDING_RESULTS.txt

hb "dataVmc nominal start"
FT=nominal bash $D/chain_dyveto_datavmc.sh > $L/dyveto_datavmc_nominal.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_datavmc_nominal.log || { hb "STOPPING: dataVmc nominal"; exit 1; }
hb "dataVmc nominal passed"
hb "DRIVER DONE"
