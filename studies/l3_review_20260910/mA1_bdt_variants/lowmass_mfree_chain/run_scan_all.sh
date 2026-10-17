#!/usr/bin/env bash
# Two lanes of scan_closure.sh, 3 cuts in parallel each (6 processes on lxplus). Marker: logs/scan_all.DONE
C=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/mA1_bdt_variants/lowmass_mfree_chain
cd $C; rm -f logs/scan_all.DONE; export NPAR=3
( bash scan_closure.sh new 1 "0.955 0.9645 0.970 0.975 0.980" > logs/scan_new_mA1.log 2>&1
  bash scan_closure.sh prod 2 "0.975" > logs/scan_prod_mA2.log 2>&1 ) &
( bash scan_closure.sh new 2 "0.975 0.9805 0.985 0.988 0.990" > logs/scan_new_mA2.log 2>&1
  bash scan_closure.sh new 3 "0.980 0.985 0.988 0.992" > logs/scan_new_mA3.log 2>&1
  bash scan_closure.sh prod 3 "0.992" > logs/scan_prod_mA3.log 2>&1 ) &
wait
date '+%F %T' > logs/scan_all.DONE
