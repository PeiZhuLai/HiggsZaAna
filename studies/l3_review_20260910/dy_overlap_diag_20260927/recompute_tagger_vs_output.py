"""Recompute MCOverlapTagger count (tagger def) and the za_tagger_resolved
n_iso_photons (output def) from the NanoAOD input of one production job, and
match to the job's output parquet by (run, lumi, event).

Input : job config + output parquet in parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/<sample>/job_<N>/
        NanoAOD via xrootd (files listed in the job config)
Output: printed table (stdout). READ-ONLY.
"""
import sys, json, glob
import numpy as np, awkward as ak, uproot, pyarrow.parquet as pq

B = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
BR = ["run", "luminosityBlock", "event", "GenPart_pdgId", "GenPart_pt", "GenPart_eta",
      "GenPart_phi", "GenPart_statusFlags"]

def dr_clean(ph, tr, mdr=0.05):
    deta = ph.eta[:, :, None] - tr.eta[:, None, :]
    dphi = (ph.phi[:, :, None] - tr.phi[:, None, :] + np.pi) % (2 * np.pi) - np.pi
    keep = ak.all(np.sqrt(deta**2 + dphi**2) > mdr, axis=-1)
    return ak.fill_none(keep, True)

def count(gp, ptcut, etacut=None):
    sel = (gp.pdgId == 22) & (gp.pt > ptcut) & (((gp.statusFlags & 0x1) != 0) | ((gp.statusFlags & 0x100) != 0))
    if etacut is not None:
        sel = sel & (abs(gp.eta) < etacut)
    ph = gp[sel]
    tr = gp[(gp.pdgId != 22) & (gp.pt > 5) & ((gp.statusFlags & 0x100) != 0)]
    return ak.num(ph[dr_clean(ph, tr)], axis=1), ph[dr_clean(ph, tr)]

def main(sample, job):
    d = f"{B}/{sample}/job_{job}"
    cfg = json.load(open(glob.glob(f"{d}/*_config_job{job}.json")[0]))
    pqt = pq.read_table(f"{d}/output_job_{job}_nominal.parquet",
                        columns=["run", "luminosityBlock", "event", "n_iso_photons"]).to_pandas()
    rows = []
    for fn in cfg["files"]:
        with uproot.open(fn + ":Events", timeout=600) as t:
            a = t.arrays(BR, library="ak", how="zip")
        gp = a.GenPart
        n_tag, _ = count(gp, 10.0)            # MCOverlapTagger Run3 def
        n_out, _ = count(gp, 15.0, 2.6)       # za_tagger_resolved.overlap_removal def
        n_tag15, _ = count(gp, 15.0)          # pt>15, no eta
        _, ph10 = count(gp, 10.0)
        lead = ak.fill_none(ak.max(ph10.pt, axis=1), -1)
        rows.append(dict(run=ak.to_numpy(a.run), lumi=ak.to_numpy(a.luminosityBlock),
                         event=ak.to_numpy(a.event), n_tag=ak.to_numpy(n_tag),
                         n_out=ak.to_numpy(n_out), n_tag15=ak.to_numpy(n_tag15),
                         lead=ak.to_numpy(lead)))
    import pandas as pd
    nano = pd.concat([pd.DataFrame(r) for r in rows])
    print(f"== {sample} job {job}: nano events={len(nano)}  file={cfg['files'][0].split('/')[-5]}")
    print(f"   nano: frac n_tag>0 = {np.mean(nano.n_tag>0):.4f}   frac n_out>0 = {np.mean(nano.n_out>0):.4f}")
    m = pqt.merge(nano, left_on=["run", "luminosityBlock", "event"], right_on=["run", "lumi", "event"], how="left")
    print(f"   parquet rows={len(pqt)} matched={m.n_tag.notna().sum()}")
    m = m[m.n_tag.notna()]
    print(f"   parquet: n_tag>0 {np.sum(m.n_tag>0)}  n_tag==0 {np.sum(m.n_tag==0)}")
    print(f"   parquet n_iso_photons == recomputed output-def : {np.mean(m.n_iso_photons==m.n_out):.4f}")
    z = m[m.n_iso_photons == 0]
    print(f"   parquet n_iso_photons==0 : {len(z)} ({len(z)/len(m):.3f}); of these n_tag>0: {np.sum(z.n_tag>0)}; "
          f"lead tagger photon 10-15 GeV: {np.sum((z.lead>10)&(z.lead<=15))}; lead>15 (so |eta|>2.6 is reason): {np.sum(z.lead>15)}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
