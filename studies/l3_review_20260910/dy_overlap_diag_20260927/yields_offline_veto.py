"""DY+jets / DY+gamma yields in 115<H_mass<135 with/without offline overlap veto via
n_iso_photons, at preselection and after the adopted BDT WP (MVAcut_points_run3.json).

Input : /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix/<sample>/<era>.root (inclusive)
        /afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json
Output: stdout table + yields_offline_veto.json in this directory. READ-ONLY on inputs.
"""
import json, os
import numpy as np, uproot

R = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
WP = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
MAS = [1, 3, 5, 10, 20, 30]
cuts = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(WP))["results"]}
DYJ = {"2022preEE": ["DYJetsToLL"], "2022postEE": ["DYJetsToLL"], "2023preBPix": ["DYJetsToLL"],
       "2023postBPix": ["DYJetsToLL"], "2024": ["DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]}
out = {}
def load(samples, era):
    brs = ["H_mass", "factor", "n_iso_photons"] + [f"MVA_Score_mA_M{m}" for m in MAS]
    parts = [uproot.open(f"{R}/{s}/{era}.root")["inclusive"].arrays(brs, library="np") for s in samples]
    return {b: np.concatenate([p[b] for p in parts]) for b in brs}

def row(a, keep):
    win = (a["H_mass"] > 115) & (a["H_mass"] < 135)
    w = a["factor"]
    res = {"presel": [w[win].sum(), w[win & keep].sum()]}
    for m in MAS:
        s = a[f"MVA_Score_mA_M{m}"] > cuts[m]
        res[f"mA{m}"] = [w[win & s].sum(), w[win & s & keep].sum(), int((win & s).sum()), int((win & s & keep).sum())]
    return res

tot = {}
for era in DYJ:
    j = load(DYJ[era], era); g = load(["DYGto2LG_10to100"], era)
    rj = row(j, j["n_iso_photons"] == 0); rg = row(g, g["n_iso_photons"] > 0)
    out[era] = {"DYjets": rj, "DYgamma": rg}
    print(f"\n=== {era}   (cols: nominal / offline-veto / removed%)   DY+jets veto: keep n_iso==0 ; DY+gamma: keep n_iso>0")
    for k in ["presel"] + [f"mA{m}" for m in MAS]:
        a, b = rj[k][:2]; c, d = rg[k][:2]
        extra = f"  [raw DYj {rj[k][2]}->{rj[k][3]}, DYg {rg[k][2]}->{rg[k][3]}]" if k != "presel" else ""
        print(f"  {k:7s} DY+jets {a:10.2f} -> {b:10.2f} ({100*(a-b)/a if a else 0:5.1f}%) | DY+gamma {c:10.2f} -> {d:10.2f} ({100*(c-d)/c if c else 0:5.1f}%) | total {a+c:10.2f} -> DYjveto-only {b+c:10.2f} ({100*(a-b)/(a+c) if a+c else 0:5.1f}%){extra}")
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "yields_offline_veto.json"), "w"), indent=1)
