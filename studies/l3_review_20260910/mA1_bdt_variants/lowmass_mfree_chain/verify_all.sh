#!/usr/bin/env bash
# Full-branch verification of every re-scored file (rescore_lowmass.py --verify-only), 4 in parallel.
# Bad outputs are moved aside to <file>.bad_<time> and re-produced with run_rescore.sh (which now runs
# the same full comparison before renaming its output). Marker: logs/verify_all.DONE
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain
O=/eos/home-p/pelai/HZa/sandbox_lowmass_mfree/run3_bdt_scored
PY=/eos/home-p/pelai/App/Conda/.conda/envs/higgs-alp-ana/bin/python3
cd $C; rm -f logs/verify_all.DONE
v() { $PY rescore_lowmass.py --sample $1 --era $2 --out-base $O --verify-only 2>&1 | tail -n 1; }
export -f v; export C O PY
xargs -P4 -L1 bash -c 'v $0 $1' < joblist_rescore.txt > logs/verify_round1.txt
grep -c "VERIFY OK" logs/verify_round1.txt > logs/verify_counts.txt
grep -v "VERIFY OK" logs/verify_round1.txt > logs/verify_bad1.txt
ts=$(date +%Y%m%d_%H%M)
while read -r line; do f=$(echo "$line" | grep -oE "$O/[^:]+\.root"); [ -n "$f" ] && mv "$f" "$f.bad_$ts"; done < logs/verify_bad1.txt
if [ -s logs/verify_bad1.txt ]; then bash run_rescore.sh; fi
xargs -P4 -L1 bash -c 'v $0 $1' < joblist_rescore.txt > logs/verify_round2.txt
echo "round1_ok=$(cat logs/verify_counts.txt) round1_bad=$(wc -l < logs/verify_bad1.txt) round2_ok=$(grep -c 'VERIFY OK' logs/verify_round2.txt) round2_bad=$(grep -vc 'VERIFY OK' logs/verify_round2.txt)" > logs/verify_all.DONE
