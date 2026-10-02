#!/usr/bin/env bash
# 2024 DY+jets overlap-veto rerun, stages A-E (plan: doc/HZa/plan_dy2024_overlap_rerun_20260927.md).
#   A  install the vetoed DY 2024 into run3_bdt_inputs_fsrfix, rebuild DYJetsToLL + All_Bkg
#   B  sideband reweight: nominal (h_m, live JSON) + S1 alternatives z_m_hsb, z_m_narrow_hsb
#   C  retrain the three BDTs (single 16-feature for dataVmc, low-mass, high-mass)
#   D  rescore all 1207 samples on condor into a FRESH run3_bdt_scored_fsrfix
#   E  hadd the scored DYJetsToLL/2024 and verify it
# Stops after E: dataVmc -> working points -> flashgg is the next chain.
#
# Old products are kept, never overwritten:
#   inputs  : the 6 replaced files -> run3_bdt_inputs_fsrfix_preDYveto_20260929/
#   scored  : whole dir renamed    -> run3_bdt_scored_fsrfix_preDYveto_20260929
#   JSON    : sideband_run3_iterative.json.bak_preDYveto_20260929
#   models  : *.bak_preDYveto_20260929 (scripts/ and using/)
#   joblist : joblist.tsv.bak_preDYveto_20260929
# Every stage is judged on its product (entries, mtimes, n_features), not exit codes.
#
# USAGE  START_STAGE=A|B|C|D|E  bash chain_dyveto_retrain.sh   (run under systemd-run)
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
REPO=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
L=$D/logs_fsrfix
R=/eos/home-p/pelai/HZa/root_P2Root
IN=$R/run3_bdt_inputs_fsrfix
NEWDY=$R/run3_bdt_inputs_fsrfix_dyveto
OLDIN=$R/run3_bdt_inputs_fsrfix_preDYveto_20260929
SC=$R/run3_bdt_scored_fsrfix
OLDSC=$R/run3_bdt_scored_fsrfix_preDYveto_20260929
RW=$REPO/HZaMVA/reweights
S=$REPO/HZaMVA/scripts
U=$REPO/HZaMVA/using
C=$REPO/Parquet2Rootfile/Condor
PQNEW=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024
TAG=bak_preDYveto_20260929
HB=$L/dyveto_retrain_heartbeat.txt
REG=$L/PENDING_RESULTS.txt
DY3="DYJetsTo2E DYJetsTo2Mu DYJetsTo2Tau"
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3

hb()  { echo "$(date '+%F %T') $*" | tee -a "$HB"; }
die() { hb "STOPPING: $*"; echo "VERDICT: FAILED ($*)"; echo "$(date '+%F %H:%M') STOPPED dyveto-retrain: $*" >> "$REG"; exit 1; }
order() { case "$1" in A) echo 1;; B) echo 2;; C) echo 3;; D) echo 4;; E) echo 5;; esac; }
START=$(order "${START_STAGE:-A}")
run() { [ "$START" -le "$(order "$1")" ]; }

set +u
source /eos/home-p/pelai/App/Anaconda/Miniconda/Install/miniconda3/etc/profile.d/conda.sh
conda activate higgs-alp-ana || { echo "conda activate failed"; exit 1; }
set -u
hb "chain start (pid $$, host $(hostname), START_STAGE=${START_STAGE:-A})"

# entries of the 'inclusive' tree
nent() { python -c "import uproot,sys; print(uproot.open(sys.argv[1])['inclusive'].num_entries)" "$1"; }

# ------------------------------------------------------------------------------------ A
if run A; then
hb "[A] gate: DY reprod orch.done must say DONE"
grep -q "^DONE" $D/dy2024_reprod_20260927/logs/orch.done || die "orch.done is not DONE"
for s in $DY3; do [ -s $NEWDY/$s/2024.root ] || die "missing $NEWDY/$s/2024.root"; done

hb "[A] move the 6 replaced files to $OLDIN"
for p in DYJetsTo2E/2024.root DYJetsTo2Mu/2024.root DYJetsTo2Tau/2024.root DYJetsToLL/2024.root DYJetsToLL/run3.root All_Bkg/run3.root; do
  mkdir -p $OLDIN/$(dirname $p)
  if [ -s $IN/$p ] && [ ! -s $OLDIN/$p ]; then mv $IN/$p $OLDIN/$p || die "mv $p"; fi
  [ -s $OLDIN/$p ] || die "old $p not preserved"
done
for s in $DY3; do cp -f $NEWDY/$s/2024.root $IN/$s/2024.root || die "cp $s"; done
mkdir -p $IN/All_Bkg
( cd $REPO/Parquet2Rootfile && HZA_P2ROOT_BASE=$IN bash 2_prepare_rootfile.sh --only dyjets ) > $L/dyveto_prepare_dyjets.log 2>&1
# all-bkg module rm -rf's All_Bkg first; the old run3.root is already in $OLDIN
( cd $REPO/Parquet2Rootfile && HZA_P2ROOT_BASE=$IN bash 2_prepare_rootfile.sh --only all-bkg ) > $L/dyveto_prepare_allbkg.log 2>&1

hb "[A] verify entries"
e2=$(nent $IN/DYJetsTo2E/2024.root); m2=$(nent $IN/DYJetsTo2Mu/2024.root); t2=$(nent $IN/DYJetsTo2Tau/2024.root)
ll=$(nent $IN/DYJetsToLL/2024.root)
[ "$ll" -eq $((e2+m2+t2)) ] || die "DYJetsToLL/2024 $ll != $e2+$m2+$t2"
llr=0; for y in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do llr=$((llr+$(nent $IN/DYJetsToLL/$y.root))); done
[ "$(nent $IN/DYJetsToLL/run3.root)" -eq "$llr" ] || die "DYJetsToLL/run3 != sum of eras"
ab=$(nent $IN/All_Bkg/run3.root); dg=$(nent $IN/DYGto2LG/run3.root)
[ "$ab" -eq $((llr+dg)) ] || die "All_Bkg $ab != DYJetsToLL $llr + DYGto2LG $dg"
oab=$(nent $OLDIN/All_Bkg/run3.root)
hb "[A] PASSED: DYJetsToLL/2024 $ll (2E $e2, 2Mu $m2, 2Tau $t2); All_Bkg $ab (old $oab, ratio $(python -c "print('%.3f'%($ab/$oab))"))"
fi

# ------------------------------------------------------------------------------------ B
if run B; then
T=$(date +%s)
hb "[B] sideband reweight: nominal + z_m_hsb + z_m_narrow_hsb"
cp -n $RW/sideband_run3_iterative.json $RW/sideband_run3_iterative.json.$TAG
cd $S
derive() {  # $1 sideband  $2 out json  $3 plot dir
  env -i HOME=$HOME PATH=/usr/bin:/bin $PY $S/1_make_sideband_reweight.py --sideband "$1" -o "$2" --plot-dir "$3" > $L/dyveto_rw_$1.log 2>&1
}
derive h_m            $RW/sideband_run3_iterative.json                 $REPO/HZaMVA/plots_MVA/sideband_reweight_run3 &
derive z_m_hsb        $RW/sideband_run3_iterative_fsrfix_dyveto_zmhsb.json       $D/s1_altsideband_20260927/validation_pdfs_dyveto/z_m_hsb &
derive z_m_narrow_hsb $RW/sideband_run3_iterative_fsrfix_dyveto_zmnarrowhsb.json $D/s1_altsideband_20260927/validation_pdfs_dyveto/z_m_narrow_hsb &
wait
for j in sideband_run3_iterative.json sideband_run3_iterative_fsrfix_dyveto_zmhsb.json sideband_run3_iterative_fsrfix_dyveto_zmnarrowhsb.json; do
  [ -s $RW/$j ] && [ "$(stat -c %Y $RW/$j)" -ge "$T" ] || die "reweight $j not written (logs dyveto_rw_*.log)"
  grep -q "All_Bkg/run3.root" $L/dyveto_rw_*.log || die "reweight did not read All_Bkg"
done
cp -f $RW/sideband_run3_iterative.json $RW/sideband_run3_iterative_fsrfix_dyveto.json
cmp -s $RW/sideband_run3_iterative.json $RW/sideband_run3_iterative.json.$TAG && die "new nominal JSON identical to old one"
hb "[B] PASSED (live JSON replaced; old kept as .$TAG)"
fi

# ------------------------------------------------------------------------------------ C
if run C; then
hb "[C] retrain three models"
for f in $S/model_Za_BDT_run3.pkl $S/model_Za_BDT_{lowmass,highmass}_run3.{pkl,json,meta.json} \
         $U/model_Za_BDT_run3.pkl $U/model_Za_BDT_{lowmass,highmass}_run3.{pkl,json,meta.json}; do
  [ -f $f ] && cp -n $f $f.$TAG
done
bash $D/chain_train_bdt.sh > $L/dyveto_train.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_train.log || die "single-model training (dyveto_train.log)"
bash $D/chain_train_lowhigh.sh > $L/dyveto_train_lowhigh.log 2>&1
grep -q "^VERDICT: PASSED" $L/dyveto_train_lowhigh.log || die "low/high training (dyveto_train_lowhigh.log)"
hb "[C] PASSED: $(grep '^\[result\]' $L/dyveto_train_lowhigh.log | tr '\n' ' ' | cut -c1-400)"
fi

# ------------------------------------------------------------------------------------ D
if run D; then
hb "[D] scoring: repoint DY 2024 parquet, fresh scored dir"
cp -n $C/joblist.tsv $C/joblist.tsv.$TAG
for s in $DY3; do
  sed -i "s|^[^\t]*/${s}_2024/merged_nominal.parquet\t|$PQNEW/${s}_2024/merged_nominal.parquet\t|" $C/joblist.tsv
  [ "$(grep -c "^$PQNEW/${s}_2024/merged_nominal.parquet" $C/joblist.tsv)" -eq 1 ] || die "joblist repoint $s"
done
[ "$(wc -l < $C/joblist.tsv)" -eq 1207 ] || die "joblist is not 1207 lines"
if [ -d $SC ] && [ ! -d $OLDSC ]; then mv $SC $OLDSC || die "mv scored dir"; fi
[ -d $OLDSC ] || die "old scored dir not preserved"
cut -f2 $C/joblist.tsv | xargs -n1 dirname | sort -u | xargs mkdir -p
[ -f $L/scoring.log ] && mv $L/scoring.log $L/scoring.log.$TAG
# SKIP_FULL_SCORING=1: resume after the full 1207 already ran (only retry what failed)
if [ "${SKIP_FULL_SCORING:-0}" != 1 ]; then
  bash $D/chain_scoring.sh > $L/scoring.log 2>&1
fi
# fresh directory: anything present was written by this stage; retry what failed (EOS I/O etc.)
for round in 1 2 3; do
  nbad=$(python $D/scoring_failed_rows.py $C/joblist.tsv 0 $C/joblist_resub.tsv)
  hb "[D] after round $((round-1)): $nbad of 1207 not OK"
  [ "$nbad" -eq 0 ] && break
  [ "$nbad" -gt 400 ] && die "$nbad scoring jobs failed -- systematic problem, not retrying (scoring.log)"
  # 2026-09-29: the retries are the big samples (Data/DYG 2024 need ~6.1 GB and >2 h), so 8 GB +
  # workday; condor_submit must run in $C because "queue ... from joblist_resub.tsv" is relative.
  sed -e 's/^request_memory .*/request_memory  = 8000MB/' -e 's/^+JobFlavour .*/+JobFlavour     = "workday"/' \
      $C/2_submit_resub.sub > $C/2_submit_resub_retry.sub
  out=$(cd $C && condor_submit 2_submit_resub_retry.sub 2>&1); cl=$(grep -oP 'submitted to cluster \K[0-9]+' <<< "$out")
  [ -n "$cl" ] || die "resub submit failed: $out"
  hb "[D] resub round $round: $nbad jobs cluster $cl (8 GB, workday)"
  while [ "$(condor_q $cl -format '%d\n' ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
    for id in $(condor_q $cl -constraint 'JobStatus==5 && RequestMemory < 16000' -af:j ClusterId 2>/dev/null | sed 's/ .*//'); do
      hb "[D] $id held: $(condor_q $id -af HoldReason | cut -c1-120) -> RequestMemory 16000 + release"
      condor_qedit $id RequestMemory 16000 >/dev/null; condor_release $id >/dev/null
    done
    sleep 300
  done
done
python $D/gate5_scored.py --joblist $C/joblist.tsv --since 0 > $L/dyveto_gate5.txt 2>&1 || die "GATE 5 on full set (dyveto_gate5.txt)"
hb "[D] PASSED: $(grep 'checked:' $L/dyveto_gate5.txt)"
fi

# ------------------------------------------------------------------------------------ E
if run E; then
hb "[E] hadd scored DYJetsToLL/2024"
( cd $C && bash 4_prepaare_2024DYJetsToLL.sh ) > $L/dyveto_scored_dyll.log 2>&1
e2=$(nent $SC/DYJetsTo2E/2024.root); m2=$(nent $SC/DYJetsTo2Mu/2024.root); t2=$(nent $SC/DYJetsTo2Tau/2024.root)
ll=$(nent $SC/DYJetsToLL/2024.root)
[ "$ll" -eq $((e2+m2+t2)) ] || die "scored DYJetsToLL/2024 $ll != $e2+$m2+$t2"
for s in $DY3; do
  pr=$(python -c "import pyarrow.parquet as pq,sys; print(pq.ParquetFile(sys.argv[1]).metadata.num_rows)" $PQNEW/${s}_2024/merged_nominal.parquet)
  [ "$(nent $SC/$s/2024.root)" -eq "$pr" ] || die "scored $s/2024 != new parquet rows $pr"
done
hb "[E] PASSED: scored DYJetsToLL/2024 $ll"
fi

echo "$(date '+%F %H:%M') DONE dyveto-retrain A-E (heartbeat $HB). Next: dataVmc -> WP -> flashgg" >> "$REG"
hb "VERDICT: PASSED"
echo "VERDICT: PASSED"
