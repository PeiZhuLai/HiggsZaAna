"""Sum per-job HiggsDNA summary n_events (events passing the MC overlap veto, before any
analysis selection) and n_events_selected for the 2024 DY samples, new (vetoed) vs old
(un-vetoed) production. Input: summary json files on EOS. Output: stdout."""
import os, json, glob, sys
from concurrent.futures import ThreadPoolExecutor
BASE = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
sets = {
  "new": BASE + "/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024",
  "old": BASE + "/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC",
}
samples = ["DYJetsTo2E_2024", "DYJetsTo2Mu_2024", "DYJetsTo2Tau_2024", "DYGto2LG_10to100_2024"]
def rd(f):
    try:
        d = json.load(open(f)); return d.get("n_events", 0), d.get("n_events_selected", {}).get("nominal", 0), 1
    except Exception as e:
        return 0, 0, 0
for tag, b in sets.items():
    for s in samples:
        d = os.path.join(b, s)
        if not os.path.isdir(d): continue
        fs = []
        for j in os.listdir(d):
            if j.startswith("job_"):
                fs += glob.glob(os.path.join(d, j, "*_summary_job*.json"))
        with ThreadPoolExecutor(32) as ex:
            r = list(ex.map(rd, fs))
        print(tag, s, "jobs_with_summary=%d n_events=%d n_selected=%d" % (sum(x[2] for x in r), sum(x[0] for x in r), sum(x[1] for x in r)), flush=True)
