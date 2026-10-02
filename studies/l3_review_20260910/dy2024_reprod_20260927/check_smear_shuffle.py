"""Matched events new vs old: fraction with identical lead-photon pT / lead-lepton pT / weight_central.
Tests the hypothesis that the seeded (seed=123) per-file smearing sequence shifts once the veto removes
events before the smearing, so surviving events get different smear values (threshold migration)."""
import os, json, glob, numpy as np, pandas as pd, pyarrow.parquet as pq
NEW="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyvetoTEST"
OLD="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
K=["run","luminosityBlock","event"]; V=["ALP_lead_photon_pt","ALP_lead_photon_mvaID","Z_lead_lepton_pt","weight_central","Z_lead_lepton_id"]
tf=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_files.json")))
fr=[]
for s,ents in tf.items():
    jm={os.path.basename(json.load(open(c))["files"][0]):json.load(open(c))["job_id"] for c in glob.glob(f"{NEW}/{s}_2024/job_*/*config*.json")}
    for oj,_,f in ents:
        j=jm[os.path.basename(f)]
        n=pq.read_table(f"{NEW}/{s}_2024/job_{j}/output_job_{j}_nominal.parquet",columns=K+V).to_pandas()
        o=pq.read_table(f"{OLD}/{s}_2024/job_{oj}/output_job_{oj}_nominal.parquet",columns=K+V).to_pandas()
        m=n.merge(o,on=K,suffixes=("_n","_o"))
        r=lambda c: np.mean(np.isclose(m[c+"_n"],m[c+"_o"],rtol=1e-6))
        print(f"{s:13s} {os.path.basename(f)[:8]} matched={len(m):3d} same pho_pt={r('ALP_lead_photon_pt'):.2f} same mvaID={r('ALP_lead_photon_mvaID'):.2f} same lep_pt={r('Z_lead_lepton_pt'):.2f} same weight={r('weight_central'):.2f} "
              f"median |dpt/pt| photon={np.median(abs(m.ALP_lead_photon_pt_n/m.ALP_lead_photon_pt_o-1)):.4f}")
