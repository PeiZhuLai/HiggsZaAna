#!/usr/bin/env bash
# Repoint the remaining downstream consumers from the pre-FSR-fix ROOT directories to
# the FSR-fix ones. Literal substitution, not an env var: these are plot/table scripts
# where the path sits inside f-strings and dict literals, and a literal swap is what
# `git diff` can show and `git revert` can undo.
#   run3_bdt_inputs_nominal -> run3_bdt_inputs_fsrfix
#   run3_bdt_scored_nominal -> run3_bdt_scored_fsrfix
# Every file is backed up and every .py is re-parsed; a file that fails to parse is
# restored rather than left broken.
set -uo pipefail
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna
FILES="
Parquet2Rootfile/A_prepare_rootfile_NFlow.sh
Parquet2Rootfile/Condor/4_prepaare_2024DYJetsToLL.sh
Parquet2Rootfile/Local/1_run_P2Root_BDT.sh
Parquet2Rootfile/Local/2_prepaare_2024DYJetsToLL.sh
Parquet2Rootfile/verify_two_model_scoring.py
Plot/lib/Analyzer_Configs.py
Plot/scripts/add_merged_flag.py
Plot/scripts/BDT_ma_2D_lib.py
Plot/scripts/make_sculpt_R_table_from_wp.py
Plot/scripts/plot_dREffBar.py
Plot/scripts/plot_dREff.py
Plot/scripts/plot_hza_angular_distributions.py
Plot/scripts/plot_mAmigratedBar.py
Plot/scripts/plot_mAmigratedHist.py
Plot/scripts/plot_mAmigratedMatrix.py
Plot/scripts/plot_MVASigEffVmA.py
Plot/scripts/plot_preselectSigEffSumwVmA.py
Plot/scripts/plot_sigGenInfo.py
Plot/scripts/plot_signal_reweight_mllgg.py
Plot/scripts/plot_training_variables.py
Plot/scripts/plot_trigEffVlepPt.py
Plot/scripts/scan_score_R_significance.py
Plot/scripts/signal_eff_sumw.py
Plot/scripts/table_interpolate_bkgYield_1.py
"
tot_before=0; tot_after=0; nfiles=0
for f in $FILES; do
  [ -f "$f" ] || { echo "MISSING  $f"; continue; }
  b=$(grep -c "run3_bdt_inputs_nominal\|run3_bdt_scored_nominal" "$f")
  [ "$b" -eq 0 ] && { echo "no-op    $f"; continue; }
  cp -n "$f" "${f}.bak_nominal_20260921"
  sed -i 's|run3_bdt_inputs_nominal|run3_bdt_inputs_fsrfix|g; s|run3_bdt_scored_nominal|run3_bdt_scored_fsrfix|g' "$f"
  a=$(grep -c "run3_bdt_inputs_nominal\|run3_bdt_scored_nominal" "$f" || true)
  ok=1
  case "$f" in
    *.py) python3 -c "import ast;ast.parse(open('$f').read())" 2>/dev/null || ok=0 ;;
    *.sh) bash -n "$f" 2>/dev/null || ok=0 ;;
  esac
  if [ "$ok" -eq 0 ]; then
    echo "BROKEN   $f -- restored"; cp "${f}.bak_nominal_20260921" "$f"; continue
  fi
  printf "patched  %-52s %d -> %d\n" "$f" "$b" "$a"
  tot_before=$((tot_before+b)); tot_after=$((tot_after+a)); nfiles=$((nfiles+1))
done
echo
echo "files patched: $nfiles   stale refs: $tot_before -> $tot_after"
# NB: `[ cond ] && { ...; }` as the last line returns false when cond is false, so the
# script would exit 1 on success. Use an explicit if.
if [ "$nfiles" -eq 0 ]; then
  echo "ERROR: nothing was patched -- the input list is empty or already clean"
  exit 1
fi
exit 0
