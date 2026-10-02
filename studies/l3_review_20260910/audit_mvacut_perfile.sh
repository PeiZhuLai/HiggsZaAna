#!/usr/bin/env bash
# per-file full read (each file in its own process, so a segfault names the file)
W=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910
B=/eos/home-p/pelai/HZa/root_MVAcut
set +u; source /cvmfs/cms.cern.ch/cmsset_default.sh; cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src && eval "$(scramv1 runtime -sh)"; set -u 2>/dev/null
{ for m in $(seq 1 30); do echo $B/data/mA_M$m/run3.root; for y in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do echo $B/data/mA_M$m/output_$y.root; done; done
  for m in 1 2 3 4 5 6 7 8 9 10 15 20 25 30; do for y in 2022preEE 2022postEE 2023preBPix 2023postBPix 2024; do echo $B/sig/mA_M$m/output_$y.root; done; done; } > $W/logs_fsrfix/mva_audit_files.txt
cat $W/logs_fsrfix/mva_audit_files.txt | xargs -P4 -I{} sh -c 'timeout 600 python3 '$W'/audit_one_mvacut.py {} >/dev/null 2>&1 && echo "OK {}" || echo "BAD rc=$? {}"' > $W/logs_fsrfix/mva_audit_perfile.txt
echo "checked $(wc -l < $W/logs_fsrfix/mva_audit_perfile.txt) bad $(grep -c ^BAD $W/logs_fsrfix/mva_audit_perfile.txt)"
grep ^BAD $W/logs_fsrfix/mva_audit_perfile.txt
