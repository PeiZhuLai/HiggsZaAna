#!/bin/bash
# Copy the AN-included figures produced by this rerun (sideband_rwgt dataVmc plots with the
# normalized signal, and the optimize_run3UL plots) into the AN tree, then cmp each one.
# The previous AN copies are kept in $W/an_fig_backup/. scan_score_R is NOT synced: its inputs
# (stored scores, unreweighted signal) do not depend on the dataVmc histograms.
set -uo pipefail
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/an_update_sigrw_20261003/datavmc_norm
AN=/afs/cern.ch/work/p/pelai/HZa/AN/AN-25-172
P=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot
T2=$(cat $W/T2_merge)
n=0; bad=0; same=0
while read -r inc; do
  case "$inc" in *scan_score_R*) continue;; esac
  rel=${inc#figure_ALP/run3/HiggsZaAna/Plot/}
  for ext in pdf png; do
    src=$P/$rel.$ext; dst=$AN/$inc.$ext
    [ -f "$src" ] || { [ "$ext" = pdf ] && { echo "MISSING source $src"; bad=$((bad+1)); }; continue; }
    [ "$(stat -c %Y "$src")" -ge "$T2" ] || { echo "STALE source (older than the merge) $src"; bad=$((bad+1)); continue; }
    if [ -f "$dst" ]; then
      mkdir -p "$W/an_fig_backup/$(dirname "$inc")"; cp -p "$dst" "$W/an_fig_backup/$inc.$ext"
      cmp -s "$src" "$dst" && same=$((same+1))
    fi
    mkdir -p "$(dirname "$dst")"; cp -p "$src" "$dst"
    cmp -s "$src" "$dst" || { echo "CMP FAILED $dst"; bad=$((bad+1)); continue; }
    n=$((n+1))
  done
done < $W/an_figure_includes.txt
echo "synced $n files (byte-identical to source after copy), $same were already identical before, problems $bad"
[ $bad -eq 0 ]
