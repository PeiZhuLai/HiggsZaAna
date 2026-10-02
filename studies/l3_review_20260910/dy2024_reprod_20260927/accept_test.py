"""Small-batch acceptance for the 2024 DY+jets re-production with the MC overlap veto.

Input : test_files.json (sample -> [old job, old rows, input file])
        NEW  : parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyvetoTEST/<s>_2024/job_*/ (config, summary, parquet)
        OLD  : parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/<s>_2024/job_<N>/output_job_<N>_nominal.parquet (no veto)
        NanoAOD of the same input file (xrootd, GenPart branches only)
        condor event logs HiggsDNA/eos_logs/Bkg_MC_dyvetoTEST/<s>_2024/job_*/*.log (memory / runtime)
Output: stdout (-> logs/accept_test.txt); last line ACCEPT_OK or ACCEPT_FAIL

Criteria
 (a) every test job has its summary json (completion marker at fpo=1; parquet may be absent)
 (b) tagger definition recomputed from NanoAOD (GenPart pdgId 22, pT>10, any eta,
     isPrompt|fromHardProcess, dR>0.05 to every other hard-process particle with pT>5):
     NO surviving new event has such a photon. Cross-check: MCOverlapTagger.overlap_selection
     called with the job's own file URL gives the same per-event cut (name matching works).
 (c) survival new/old for the same files; expected ~0.65-0.80 per the diagnosis. Also
     new set == old set restricted to n_tag==0 (the veto only REMOVES events).
Extra: weight_central new/old on matched events (SF payloads changed on 2026-09-27; a ratio != 1
       means the new sample carries newer SFs than the rest of the production).
"""
import os, sys, json, glob, re
import numpy as np, awkward as ak, uproot, pandas as pd, pyarrow.parquet as pq

S = os.path.dirname(os.path.abspath(__file__))
NEW = sys.argv[1] if len(sys.argv) > 1 else "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyvetoTEST"
OLD = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
LOGS = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/eos_logs/" + os.path.basename(NEW)
BR = ["run", "luminosityBlock", "event", "GenPart_pdgId", "GenPart_pt", "GenPart_eta",
      "GenPart_phi", "GenPart_mass", "GenPart_status", "GenPart_statusFlags", "GenPart_genPartIdxMother"]
KEY = ["run", "luminosityBlock", "event"]

def n_tag(gp):
    sel = (gp.pdgId == 22) & (gp.pt > 10) & (((gp.statusFlags & 0x1) != 0) | ((gp.statusFlags & 0x100) != 0))
    ph = gp[sel]
    tr = gp[(gp.pdgId != 22) & (gp.pt > 5) & ((gp.statusFlags & 0x100) != 0)]
    deta = ph.eta[:, :, None] - tr.eta[:, None, :]
    dphi = (ph.phi[:, :, None] - tr.phi[:, None, :] + np.pi) % (2 * np.pi) - np.pi
    keep = ak.fill_none(ak.all(np.sqrt(deta**2 + dphi**2) > 0.05, axis=-1), True)
    return ak.to_numpy(ak.num(ph[keep], axis=1))

def cols(p, want):
    if not os.path.exists(p):
        return pd.DataFrame(columns=want)
    have = pq.ParquetFile(p).schema_arrow.names
    return pq.read_table(p, columns=[c for c in want if c in have]).to_pandas()

def condor_mem(jobdir):
    out = []
    for lg in glob.glob(jobdir + "/*.log"):
        txt = open(lg, errors="ignore").read()
        m = re.findall(r"Memory \(MB\)\s*:\s*(\d+)\s+(\d+)\s+(\d+)", txt)
        t = re.findall(r"Total Remote Usage\s*\n?\s*Usr 0 (\d+):(\d+):(\d+)", txt)
        wall = re.findall(r"Run Remote Usage", txt)
        if m:
            out.append("%s usage=%sMB req=%s" % (os.path.basename(lg), m[-1][0], m[-1][1]))
    return "; ".join(out) if out else "no-condor-log-memory"

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA")
from higgs_dna.taggers.mc_overlap_tagger import MCOverlapTagger

tf = json.load(open(os.path.join(S, "test_files.json")))
ok = True
tot_old = tot_new = 0
for s, entries in tf.items():
    d = f"{NEW}/{s}_2024"
    jobs = {}
    for cfgp in glob.glob(f"{d}/job_*/{s}_2024_config_job*.json"):
        c = json.load(open(cfgp))
        for f in c["files"]:
            jobs[os.path.basename(f)] = (c["job_id"], f)
    for oldjob, oldrows, f in entries:
        b = os.path.basename(f)
        if b not in jobs:
            print(f"[{s}] {b}: NO new job uses this file -> FAIL"); ok = False; continue
        j, furl = jobs[b]
        jd = f"{d}/job_{j}"
        summ = glob.glob(f"{jd}/*_summary_job{j}.json")
        a_ok = bool(summ)
        newp = f"{jd}/output_job_{j}_nominal.parquet"
        oldp = f"{OLD}/{s}_2024/job_{oldjob}/output_job_{oldjob}_nominal.parquet"
        new = cols(newp, KEY + ["weight_central", "n_iso_photons"])
        old = cols(oldp, KEY + ["weight_central", "n_iso_photons"])
        # NanoAOD, tagger definition
        with uproot.open(furl + ":Events", timeout=900) as t:
            arr = t.arrays(BR, library="ak", how="zip")
            tagger = MCOverlapTagger(is_data=False, year="2024")
            cut_code = ak.to_numpy(tagger.overlap_selection(furl, t))
        nt = n_tag(arr.GenPart)
        nano = pd.DataFrame({"run": ak.to_numpy(arr.run), "luminosityBlock": ak.to_numpy(arr.luminosityBlock),
                             "event": ak.to_numpy(arr.event), "n_tag": nt, "cut_code": cut_code})
        code_agrees = bool(np.all(nano.cut_code == (nano.n_tag == 0)))
        mn = new.merge(nano, on=KEY, how="left")
        mo = old.merge(nano, on=KEY, how="left")
        b_bad = int((mn.n_tag > 0).sum()); b_unmatched = int(mn.n_tag.isna().sum())
        old_keep = set(map(tuple, mo[mo.n_tag == 0][KEY].values.tolist()))
        new_set = set(map(tuple, new[KEY].values.tolist()))
        set_eq = old_keep == new_set
        w = new.merge(old, on=KEY, suffixes=("_new", "_old"))
        wr = (w.weight_central_new / w.weight_central_old).to_numpy() if len(w) else np.array([np.nan])
        surv = len(new) / len(old) if len(old) else float("nan")
        tot_old += len(old); tot_new += len(new)
        this_ok = a_ok and b_bad == 0 and b_unmatched == 0 and code_agrees and set_eq
        ok &= this_ok
        print(f"[{s}] newjob={j} oldjob={oldjob} file={b}\n"
              f"   (a) summary={'yes' if a_ok else 'NO'}  new parquet={'yes' if os.path.exists(newp) else 'no'}\n"
              f"   nano events={len(nano)}  frac n_tag>0={np.mean(nt>0):.3f}  tagger-code cut == (n_tag==0): {code_agrees} "
              f"(code keeps {int(cut_code.sum())}/{len(cut_code)})\n"
              f"   (b) new rows={len(new)}  with iso gen photon={b_bad}  unmatched to nano={b_unmatched}\n"
              f"   (c) old rows={len(old)}  old with n_tag==0={len(old_keep)}  survival new/old={surv:.3f}  new==old&(n_tag==0): {set_eq}\n"
              f"   weight_central new/old on {len(w)} matched: min={np.nanmin(wr):.4f} median={np.nanmedian(wr):.4f} max={np.nanmax(wr):.4f}\n"
              f"   condor: {condor_mem(LOGS + f'/{s}_2024/job_{j}')}\n"
              f"   -> {'OK' if this_ok else 'FAIL'}")
print(f"TOTAL old={tot_old} new={tot_new} survival={tot_new/max(tot_old,1):.3f}")
print("ACCEPT_OK" if ok else "ACCEPT_FAIL")
sys.exit(0 if ok else 1)
