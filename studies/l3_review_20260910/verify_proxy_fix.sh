#!/usr/bin/env bash
# Watch for the first chunks produced after the proxy copy was refreshed.
# Baseline is taken at start; success = more than +30 chunks within 30 minutes.
set -uo pipefail
S=/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix
L=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/logs_fsrfix/verify_proxy_fix.log
c(){ find $S/Bkg_MC $S/Data -name '*_nominal.parquet' ! -name 'merged_nominal.parquet' 2>/dev/null | wc -l; }
b=$(c); echo "$(date '+%F %T') baseline=$b" >> "$L"
end=$(( $(date +%s) + 1800 ))
while [ "$(date +%s)" -lt "$end" ]; do
  sleep 180
  n=$(c); echo "$(date '+%F %T') chunks=$n (+$((n-b)))" >> "$L"
  if [ "$n" -gt $((b+30)) ]; then echo "$(date '+%F %T') PROXY FIX CONFIRMED: $b -> $n" >> "$L"; exit 0; fi
done
echo "$(date '+%F %T') NOT CONFIRMED after 30 min: $b -> $(c)" >> "$L"
