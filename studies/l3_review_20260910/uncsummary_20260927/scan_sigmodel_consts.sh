#!/bin/bash
# Generalization of HiggsDNA/signal_rerun_logs/fnuf2022/scan_fnuf_consts.sh (the source of the AN's
# FNUF peak/resolution numbers) to all shape nuisances, on the CURRENT (FSR-fix) signal models.
# GeV conversion as in the original: |const_mean| * max(mean_g*), |const_sigma| * max(sigma_g*).
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src; eval `scramv1 runtime -sh` >/dev/null 2>&1
python3 - <<'PY'
import ROOT, os, json, statistics as st
ROOT.gErrorIgnoreLevel = ROOT.kFatal
base="/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Signal"
mAs=[1,2,3,4,5,6,7,8,9,10,15,20,25,30]
years=["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
rows=[]
for ch in ["ele","mu"]:
  for mA in mAs:
    for y in years:
      fp=f"{base}/outdir_{ch}/signalFit/output/{mA}_CMS-HGG_sigfit_{y}_{ch}_Hm125.root"
      f=ROOT.TFile(fp); w=f.Get("wsig_13p6TeV")
      pref=f"GG2H_{y}_cat0_13p6TeV"
      means=[w.obj(f"mean_g{i}_{pref}").getVal() for i in range(6) if w.obj(f"mean_g{i}_{pref}")]
      sigmas=[w.obj(f"sigma_g{i}_{pref}").getVal() for i in range(6) if w.obj(f"sigma_g{i}_{pref}")]
      rec=dict(ch=ch,mA=mA,year=y,mu0=max(means),sg0=max(sigmas),mtime=os.path.getmtime(fp))
      for syst in ["FNUF_13p6TeVscaleCorr","Material_13p6TeVscaleCorr",f"ElectronScale_13p6TeVscale_{y}",f"ElectronSmear_13p6TeVsmear_{y}",
                   f"MuonScale_13p6TeVscale_{y}",f"MuonSmear_13p6TeVsmear_{y}",f"PhotonScale_13p6TeVscale_{y}",f"PhotonSmear_13p6TeVsmear_{y}"]:
        key=syst.split("_")[0]
        for stub in ["mean","sigma","rate"]:
          o=w.obj(f"const_{pref}_{stub}_{syst}")
          rec[f"{key}_{stub}"]=o.getVal() if o else 0.0
      rows.append(rec); f.Close()
json.dump(rows,open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/uncsummary_20260927/sigmodel_consts.json","w"),indent=1)
import time
print("models:",len(rows)," oldest model mtime:",time.ctime(min(r['mtime'] for r in rows)))
for key in ["FNUF","Material","ElectronScale","ElectronSmear","MuonScale","MuonSmear","PhotonScale","PhotonSmear"]:
  for sel,lab in [(lambda r:r['mA']>=2,"mA>=2"),(lambda r:r['mA']==1,"mA=1")]:
    rr=[r for r in rows if sel(r)]
    pk=[abs(r[f"{key}_mean"])*r['mu0'] for r in rr]; rs=[abs(r[f"{key}_sigma"])*r['sg0'] for r in rr]; rt=[r[f"{key}_rate"] for r in rr]
    ncap=sum(1 for v in rt if v>=0.0499)
    print(f"{key:14s} {lab:6s} peak GeV median {st.median(pk):.3f} max {max(pk):.3f} | res GeV median {st.median(rs):.3f} max {max(rs):.3f} | model rate const median {100*st.median(rt):.2f}% max {100*max(rt):.2f}% (at 5% cap: {ncap}/{len(rt)})")
PY
