import sys, numpy as np, uproot
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import SidebandReweighter
RW = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/"
new = SidebandReweighter.from_json(RW + "sideband_run3_iterative.json")
old = SidebandReweighter.from_json(RW + "sideband_run3_iterative.json.bak_preDYveto_20260929")
for m in (1, 5, 15, 30):
    tot = {}
    for era in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"):
        a = uproot.open(f"/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/mA_M{m}/{era}.root")["inclusive"].arrays(library="pd")
        a = a.rename(columns={"pho1PIso_noCorr": "pho1ECALIso", "pho2PIso_noCorr": "pho2ECALIso"})
        w = a["weight"].to_numpy(float)
        h = a.drop(columns=[c for c in ("param",) if c in a]).copy()           # event-hash
        t = a.copy(); t["param"] = (t["ALP_m"] - m) / t["H_m"]                  # true mass
        for k, rw, fr in (("old_hash", old, h), ("new_hash", new, h), ("new_true", new, t)):
            tot[k] = tot.get(k, 0) + np.sum(w * rw.weights_for_dataframe(fr.copy()))
        tot["none"] = tot.get("none", 0) + np.sum(w)
    print("mA%-2d  sum(w*R)/sum(w):  old JSON+hash %.3f   new JSON+hash %.3f   new JSON+true mass %.3f"
          % (m, tot["old_hash"] / tot["none"], tot["new_hash"] / tot["none"], tot["new_true"] / tot["none"]))
