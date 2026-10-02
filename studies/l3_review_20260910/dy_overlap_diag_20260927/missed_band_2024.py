"""For 2024 DY+jets (veto NOT applied in production): among preselected parquet events,
compare the tagger definition (pT>10, no eta) with the stored n_iso_photons (pT>15,|eta|<2.6),
to size what an offline n_iso_photons==0 veto would miss.
Input : parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC/<sample>/job_<N>/ + NanoAOD via xrootd. READ-ONLY.
Output: stdout.
"""
import sys, glob, json
import numpy as np, awkward as ak, uproot, pyarrow.parquet as pq, pandas as pd
sys.path.insert(0, ".")
from recompute_tagger_vs_output import B, BR, count
sample, jobs = sys.argv[1], sys.argv[2].split(",")
agg = dict(n=0, tag=0, out=0, tag_not_out=0, w=0., wtag=0., wout=0.)
for job in jobs:
    d = f"{B}/{sample}/job_{job}"
    cfg = json.load(open(glob.glob(f"{d}/*_config_job{job}.json")[0]))
    cols = [c for c in ["run", "luminosityBlock", "event", "n_iso_photons", "weight_central"]
            if c in pq.ParquetFile(f"{d}/output_job_{job}_nominal.parquet").schema_arrow.names]
    p = pq.read_table(f"{d}/output_job_{job}_nominal.parquet", columns=cols).to_pandas()
    if "weight_central" not in p: p["weight_central"] = 1.0
    for fn in cfg["files"]:
        with uproot.open(fn + ":Events", timeout=600) as t:
            a = t.arrays(BR, library="ak", how="zip")
        n_tag, _ = count(a.GenPart, 10.0)
        nano = pd.DataFrame(dict(run=ak.to_numpy(a.run), luminosityBlock=ak.to_numpy(a.luminosityBlock),
                                 event=ak.to_numpy(a.event), n_tag=ak.to_numpy(n_tag)))
        m = p.merge(nano, on=["run", "luminosityBlock", "event"], how="inner")
        agg["n"] += len(m); agg["tag"] += int((m.n_tag > 0).sum()); agg["out"] += int((m.n_iso_photons > 0).sum())
        agg["tag_not_out"] += int(((m.n_tag > 0) & (m.n_iso_photons == 0)).sum())
        agg["w"] += m.weight_central.sum(); agg["wtag"] += m.weight_central[m.n_tag > 0].sum()
        agg["wout"] += m.weight_central[m.n_iso_photons > 0].sum()
    print(sample, "after job", job, agg, flush=True)
print(f"FINAL {sample}: preselected={agg['n']}  tagger-def>0: {agg['tag']} ({agg['tag']/agg['n']:.3f}, w {agg['wtag']/agg['w']:.3f})  "
      f"stored n_iso>0: {agg['out']} ({agg['out']/agg['n']:.3f}, w {agg['wout']/agg['w']:.3f})  missed by offline veto: {agg['tag_not_out']}")
