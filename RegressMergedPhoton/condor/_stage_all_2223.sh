#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
export X509_USER_PROXY_AFS=/afs/cern.ch/user/p/pelai/.x509up_pelai
cd /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor
CAT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json
export FILES_PER_JOB=8
TOT_FILES=0; TOT_JOBS=0
printf "%-42s %8s %6s\n" "TAG" "FILES" "JOBS"
for era in 2022preEE 2022postEE 2023preBPix 2023postBPix; do
  mapfile -t entries < <(python3 - "$CAT" "$era" <<'PY'
import json,sys,re
d=json.load(open(sys.argv[1])); era=sys.argv[2]
for path in d["Data_MINI"]["files"].get(era,[]):
    m=re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+",path)
    if not m: continue
    primary,_,run,rest=m.group(1),m.group(2),m.group(3),m.group(4)
    vers="".join(re.findall(r"[-_]v\d+",rest)).replace("-","").replace("_","")
    print(f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}|{path}")
PY
)
  for e in "${entries[@]}"; do
    tag="${e%%|*}"; das="${e##*|}"
    out=$(FILES_PER_JOB=8 bash submit_mass_point.sh "$tag" "$das" data "$era" 2>&1)
    nf=$(echo "$out" | grep -oE "Got [0-9]+ files" | grep -oE "[0-9]+")
    nj=$(echo "$out" | grep -oE "submit [0-9]+ condor jobs" | grep -oE "[0-9]+")
    nf=${nf:-0}; nj=${nj:-0}
    printf "%-42s %8s %6s\n" "$tag" "$nf" "$nj"
    TOT_FILES=$((TOT_FILES+nf)); TOT_JOBS=$((TOT_JOBS+nj))
  done
done
echo "-----------------------------------------------------------"
printf "%-42s %8s %6s\n" "TOTAL (36 datasets)" "$TOT_FILES" "$TOT_JOBS"
