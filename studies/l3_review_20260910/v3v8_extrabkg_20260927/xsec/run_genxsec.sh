#!/bin/bash
# GenXSecAnalyzer on 3 MiniAOD files per sample, to cross-check the xs used for
# ZGG (McM 0.2147 pb), TTGG (HZgamma catalog 0.02391 vs McM 0.56) and TTG (catalog 4.629).
export X509_USER_PROXY=/tmp/x509up_u175325
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/CMSSW_15_0_14/src && eval $(scram runtime -sh) && cd - >/dev/null
for f in files_*.txt; do
  s=${f#files_}; s=${s%.txt}
  echo "=== $s $(date)"
  cmsRun genxsec_cfg.py $f > genxsec_$s.log 2>&1
  grep -A3 "Final cross section\|After filter: final cross section\|Filter efficiency (event-level)" genxsec_$s.log | head -20
done
echo DONE
