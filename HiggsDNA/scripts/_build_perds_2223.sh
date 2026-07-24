#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh >/dev/null 2>&1
CAT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/metadata/samples/zgamma_tutorial_sample_manager_full.json
OUT=/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/metadata/samples/za_merged_data_2223_perds.json
python3 - "$CAT" "$OUT" <<'PY'
import json,sys,re,subprocess
cat,out=sys.argv[1],sys.argv[2]
d=json.load(open(cat))
res={}; missing=[]
for era in ["2022preEE","2022postEE","2023preBPix","2023postBPix"]:
    for mini in d["Data_MINI"]["files"].get(era,[]):
        m=re.match(r"^/([^/]+)/Run(\d{4})([A-Z])([^/]*)/.+",mini)
        if not m: continue
        primary,_,run,rest=m.group(1),m.group(2),m.group(3),m.group(4)
        vers="".join(re.findall(r"[-_]v\d+",rest)).replace("-","").replace("_","")
        tag=f"Data_{era}_{primary}_Run{run}_{vers or 'novers'}"
        # child NanoAOD of this MiniAOD, pick plain -NanoAODv15-
        ch=subprocess.run(["dasgoclient","-query",f"child dataset={mini}"],
                          capture_output=True,text=True,timeout=120).stdout.split()
        nano=[c for c in ch if re.search(r"-NanoAODv15-v\d+/NANOAOD$",c)]
        if not nano:
            # fallback: any *NanoAODv15* NANOAOD excluding BTV/JME
            nano=[c for c in ch if "NanoAODv15" in c and c.endswith("/NANOAOD")
                  and "BTV" not in c and "JME" not in c]
        if nano:
            res[tag]={"fpo":1,"files":{era:[sorted(nano)[0]]}}
        else:
            missing.append((tag,mini,ch[:4]))
json.dump(res,open(out,"w"),indent=1)
print(f"wrote {len(res)} tags -> {out}")
if missing:
    print(f"\n⚠️ {len(missing)} tags WITHOUT plain NanoAODv15:")
    for t,mi,ch in missing: print(f"  {t}\n    mini={mi}\n    children={ch}")
PY
