#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
export X509_USER_PROXY_AFS=/afs/cern.ch/user/p/pelai/.x509up_pelai
CD=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/RegressMergedPhoton/condor
cd "$CD"
CAT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json
export FILES_PER_JOB=8

# 1) clear the 3 fixed Muon1 stale empty stages so they re-query DAS
for t in Data_2023preBPix_Muon1_RunC_v1v1 Data_2023preBPix_Muon1_RunC_v2v1 Data_2023preBPix_Muon1_RunC_v3v1; do
  rm -f "$CD/stage/$t/files.txt" "$CD/stage/$t/args.txt"
  echo "cleared stale stage: $t"
done

# 2) re-stage the 3 fixed datasets (fresh DAS query w/ corrected path)
for v in v1 v2 v3; do
  das=$(python3 -c "import json;d=json.load(open('$CAT'));print([p for p in d['Data_MINI']['files']['2023preBPix'] if '/Muon1/Run2023C-22Sep2023_${v}-v1' in p][0])")
  bash submit_mass_point.sh "Data_2023preBPix_Muon1_RunC_${v}v1" "$das" data 2023preBPix 2>&1 | grep -E "Got|Will submit"
done

# 3) submit ALL staged datasets (skip any with empty/zero args)
echo "=========== SUBMITTING ALL ==========="
TOTJOBS=0; NCL=0
for era in 2022preEE 2022postEE 2023preBPix 2023postBPix; do
  mapfile -t tags < <(python3 - "$CAT" "$era" <<'PY'
import json,sys,re
d=json.load(open(sys.argv[1])); era=sys.argv[2]
for path in d["Data_MINI"]["files"].get(era,[]):
    m=re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+",path)
    if not m: continue
    primary,_,run,rest=m.group(1),m.group(2),m.group(3),m.group(4)
    vers="".join(re.findall(r"[-_]v\d+",rest)).replace("-","").replace("_","")
    print(f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}")
PY
)
  for tag in "${tags[@]}"; do
    jdl="$CD/stage/$tag/job.jdl"; args="$CD/stage/$tag/args.txt"
    if [ ! -s "$args" ]; then echo "  [SKIP empty] $tag"; continue; fi
    nj=$(wc -l < "$args")
    res=$(condor_submit "$jdl" 2>&1 | tail -1)
    echo "  [$tag] $nj jobs -> $res"
    TOTJOBS=$((TOTJOBS+nj)); NCL=$((NCL+1))
  done
done
echo "======================================"
echo "TOTAL submitted: $NCL clusters, $TOTJOBS jobs"
