#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
export X509_USER_PROXY=/afs/cern.ch/user/p/pelai/.x509up_pelai
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor
for era in 2022preEE 2022postEE 2023preBPix 2023postBPix; do
  echo "########## RESUME $era ##########"
  FILES_PER_JOB=8 bash resume_missing_data.sh "$era" 2>&1 | grep -E "RESUB|SKIP|submitted to cluster|jobs submitted"
done
echo "########## RESUME 完成 ##########"
