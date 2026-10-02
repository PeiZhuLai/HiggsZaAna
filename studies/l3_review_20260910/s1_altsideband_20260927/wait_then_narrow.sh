#!/bin/bash
S=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/s1_altsideband_20260927
until grep -q DONE $S/logs/derive_z_m.log 2>/dev/null; do sleep 15; done
$S/run_derive.sh z_m_narrow /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative_fsrfix_zmnarrow.json > $S/logs/derive_z_m_narrow.log 2>&1
