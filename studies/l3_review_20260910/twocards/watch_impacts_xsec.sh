#!/usr/bin/env bash
# Waits for the cross-section-card impact cluster ($1), checks the 30 plots and json in Combine/output_impacts_xsec,
# writes a DONE/STOP line in the varsearch registry. Held jobs are left to periodic_release; all-held -> STOP.
CL=$1; C=/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine
REG=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/logs/PENDING_RESULTS.txt
ST=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/twocards/logs/impacts_xsec.status
while [ "$(condor_q $CL -af ProcId 2>/dev/null | wc -l)" -gt 0 ]; do
  s=$(condor_q $CL -af JobStatus 2>/dev/null); h=$(echo "$s" | grep -c '^5$'); r=$(echo "$s" | grep -c '^2$'); t=$(echo "$s" | grep -c .)
  echo "$(date '+%F %T') impacts xsec cluster $CL: $t in queue, $r running, $h held" > $ST
  [ "$t" -gt 0 ] && [ "$h" -eq "$t" ] && { echo "$(date '+%F %H:%M') STOP  xsec impacts: all $t remaining jobs held ($C/impact_logs)" >> $REG; exit 1; }
  sleep 600; done
miss=""; for m in $(seq 1 30); do [ -s $C/output_impacts_xsec/${m}_Expected_Impacts.pdf ] && [ -s $C/output_impacts_xsec/${m}_impacts.json ] || miss="$miss $m"; done
if [ -z "$miss" ]; then echo "$(date '+%F %H:%M') DONE  xsec impacts (30/30 plots)  -> $C/output_impacts_xsec" >> $REG
else echo "$(date '+%F %H:%M') STOP  xsec impacts missing mA:$miss ($C/impact_logs/impact_xsec_mA*.err)" >> $REG; fi
echo "$(date '+%F %T') finished" > $ST
