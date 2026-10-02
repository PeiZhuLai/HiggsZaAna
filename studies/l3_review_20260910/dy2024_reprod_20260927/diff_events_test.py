"""Events that differ between NEW (veto) and OLD (no veto) test outputs beyond the veto itself.
Input : test_files.json, new/old job parquets.  Output: stdout. Uses the old n_iso_photons? no --
we compare sets on (run,lumi,event) and print photon/lepton kinematics of the symmetric difference."""
import os, json, glob, pandas as pd, pyarrow.parquet as pq
NEW="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyvetoTEST"
OLD="/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
K=["run","luminosityBlock","event"]
C=K+["ALP_lead_photon_pt","ALP_sublead_photon_pt","ALP_lead_photon_mvaID","ALP_sublead_photon_mvaID","Z_mass","H_mass","ALP_mass","Z_lead_lepton_pt","Z_sublead_lepton_pt","n_iso_photons","weight_central"]
tf=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),"test_files.json")))
pd.set_option("display.width",250); pd.set_option("display.max_columns",30)
for s,ents in tf.items():
    jm={}
    for c in glob.glob(f"{NEW}/{s}_2024/job_*/*config*.json"):
        cc=json.load(open(c)); jm[os.path.basename(cc["files"][0])]=cc["job_id"]
    for oj,_,f in ents:
        j=jm[os.path.basename(f)]
        n=pq.read_table(f"{NEW}/{s}_2024/job_{j}/output_job_{j}_nominal.parquet",columns=C).to_pandas()
        o=pq.read_table(f"{OLD}/{s}_2024/job_{oj}/output_job_{oj}_nominal.parquet",columns=C).to_pandas()
        m=n.merge(o,on=K,how="outer",suffixes=("_new","_old"),indicator=True)
        only_new=m[m._merge=="left_only"]
        print(f"== {s} file {os.path.basename(f)}: new={len(n)} old={len(o)} both={int((m._merge=='both').sum())} new_only={len(only_new)}")
        if len(only_new):
            print(only_new[[c for c in m.columns if c.endswith('_new')]+K].to_string())
            # is it in old at all? no (left_only). show nothing else
