#!/usr/bin/env bash
# User-approved cleanup (2026-09-29, categories 1-4):
#   1  pre-FSR-fix ROOT         (EOS home root_P2Root *_nominal*, root_MVAcut *_preFSRfix_* and mA3 WP variants)
#   2  pre-FSR-fix parquet      (EOS project parquet_DNA, parquet_DNA_tmp, parquet_DNA_tmp_fsrfix/{Bkg_MC,Data};
#                                EOS home parquet_DNA)
#   3  chunk parquets of the CURRENT production, merged_*.parquet kept -- only for sample dirs that pass
#      chunk_rows == merged_rows and merged mtime >= max(chunk mtime) for EVERY merged file (cleanup_chunks.py)
#   4  test dirs + AFS condor job logs of finished productions
# Kept on purpose: HZa_merged (Sig_MC_MLNANO_all), Zee/Zmmg control, parquet_eveto, TnP, *_preDYveto_20260929,
#   the vetoed DY production, everything the current chains read.
# `eos rm -r` refuses big trees but returns rc=0 (memory ref_eos_rm_recursive_query_limit), so every target
# is checked with `eos stat` afterwards and split into children when it survived.
#
# USAGE  bash cleanup_old_20260929.sh            # dry run: list targets with sizes
#        APPLY=1 bash cleanup_old_20260929.sh    # delete
set -uo pipefail
D=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
LOG=$D/logs_fsrfix/cleanup_20260929.log
APPLY=${APPLY:-0}
PRJ=root://eosproject-h.cern.ch
USR=root://eosuser.cern.ch
P=/eos/project/h/htozg-dy-privatemc/pelai/HZa
H=/eos/user/p/pelai/HZa
say(){ echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

T1="$USR $H/root_P2Root/run3_bdt_scored_nominal
$USR $H/root_P2Root/run3_bdt_scored_nominal.old_preSignalReweight_20260706
$USR $H/root_P2Root/run3_bdt_inputs_nominal
$USR $H/root_MVAcut/sig_preFSRfix_20260924
$USR $H/root_MVAcut/data_preFSRfix_20260924
$USR $H/root_MVAcut/data_variants_mA_M3_before_wp0940_20260926_2133
$USR $H/root_MVAcut/data_variants_mA_M3_before_wp0992_20260927_0308
$USR $H/root_MVAcut/data_variants_mA_M3_wp0940
$USR $H/root_MVAcut/sig_variants_mA_M3_before_wp0940_20260926_2133
$USR $H/root_MVAcut/sig_variants_mA_M3_before_wp0992_20260927_0308
$USR $H/root_MVAcut/sig_variants_mA_M3_wp0940"
T2="$PRJ $P/parquet_DNA
$PRJ $P/parquet_DNA_tmp
$PRJ $P/parquet_DNA_tmp_fsrfix/Bkg_MC
$PRJ $P/parquet_DNA_tmp_fsrfix/Data
$USR $H/parquet_DNA"
T4="$PRJ $P/parquet_DNA_tmp_fsrfix_extraBkg/Bkg_MC_extraBkgTEST
$USR $H/root_P2Root/run3_bdt_inputs_fsrfix_extraBkgTEST
$USR $H/root_P2Root/run3_bdt_scored_fsrfix_extraBkgTEST
$USR $H/smoketest_cutflow
$USR $H/smoketest_cutflow2"
EL=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs
# finished HZa productions only; Sig_MC_MLNANO_all (HZa_merged) and the zee_zmmg logs stay
AFSLOGS=$(ls -d $EL/*/ $EL/resub3_*.log 2>/dev/null | grep -v "Sig_MC_MLNANO_all")

size(){ EOS_MGM_URL=$1 timeout 120 eos ls -lh "$(dirname "$2")" 2>/dev/null | awk -v n="$(basename "$2")" '$NF==n{print $5,$6}'; }
exists(){ EOS_MGM_URL=$1 eos stat "$2" >/dev/null 2>&1; }

rmtree(){   # $1 mgm  $2 path
  local m=$1 p=$2 c
  EOS_MGM_URL=$m eos rm -r "$p" >> "$LOG" 2>&1
  exists $m "$p" || return 0
  say "  split: $p survived eos rm -r"
  local files; files=$(mktemp -p $D/logs_fsrfix rmfiles.XXXX)
  EOS_MGM_URL=$m eos ls -l "$p" 2>/dev/null | while read -r l; do
    c=$(awk '{print $NF}' <<< "$l")
    case "$l" in d*) rmtree $m "$p/$c" ;; *) echo "rm \"$p/$c\"" >> "$files" ;; esac
  done
  [ -s "$files" ] && EOS_MGM_URL=$m eos < "$files" >> "$LOG" 2>&1
  rm -f "$files"
  EOS_MGM_URL=$m eos rmdir "$p" >> "$LOG" 2>&1 || EOS_MGM_URL=$m eos rm -r "$p" >> "$LOG" 2>&1
  exists $m "$p" && { say "  FAILED to remove $p"; return 1; }
  return 0
}

run_list(){  # $1 label  $2 list
  local n=0 bad=0 m p
  while read -r m p; do
    [ -z "${p:-}" ] && continue
    if ! exists $m "$p"; then say "[$1] already gone: $p"; continue; fi
    say "[$1] $(size $m "$p")  $p"
    n=$((n+1))
    if [ "$APPLY" = 1 ]; then rmtree $m "$p" && say "[$1]   removed" || bad=$((bad+1)); fi
  done <<< "$2"
  say "[$1] targets=$n failed=$bad"
}

say "===== cleanup start APPLY=$APPLY ====="
run_list 1 "$T1"
run_list 2 "$T2"
run_list 4 "$T4"
for d in $AFSLOGS; do
  say "[4-afs] $(du -sh "$d" 2>/dev/null | cut -f1)  $d"
  [ "$APPLY" = 1 ] && { rm -rf "$d"; [ -e "$d" ] && say "  FAILED $d"; }
done
say "[3] chunk parquets: reconcile then delete (cleanup_chunks.py APPLY=$APPLY)"
/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3 $D/cleanup_chunks.py $APPLY >> "$LOG" 2>&1
say "[3] cleanup_chunks.py rc=$? (report logs_fsrfix/cleanup_chunks_report.txt)"
say "===== cleanup end ====="
