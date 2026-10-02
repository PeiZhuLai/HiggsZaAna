#!/usr/bin/env python3
"""READ-ONLY. Emulate fixed split variants on the scored INCLUSIVE trees (after own MVA cut, per lepton):
  INC : whole inclusive tree (no split) syst/nominal yield ratio
  HASH: key-stable split bucket = hash(run,lumi,event)%10 >= 7 applied to BOTH nominal and syst, own MVA cut,
        no join to the nominal pass set (what a proper fix would give)
Also checks whether Electron_scale_X and Electron_smear_X are identical (lepton pT arrays).
"""
import sys, json, os
import numpy as np
import uproot

SCORED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
CUTS = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
cuts = {int(r["mA"]): float(r["MVAcut"]) for r in json.load(open(CUTS))["results"]}
SYSTS = ["FNUF", "Material", "Electron_scale", "Electron_smear", "Muon_scale", "Muon_smear", "Photon_scale", "Photon_smear"]


def hbucket(df):
    # splitmix64-style mix of (run, lumi, event) -> deterministic, order independent
    x = (df.run.to_numpy().astype(np.uint64) * np.uint64(0x9E3779B97F4A7C15)
         ^ df.luminosityBlock.to_numpy().astype(np.uint64) * np.uint64(0xBF58476D1CE4E5B9)
         ^ df.event.to_numpy().astype(np.uint64) * np.uint64(0x94D049BB133111EB))
    x ^= x >> np.uint64(31)
    return (x % np.uint64(10)).astype(int)


def yields(df, mva, cut, test_only):
    d = df[df[mva] > cut]
    if test_only:
        d = d[hbucket(d) >= 7]
    return {lep: float(d[d[col] == 2].weight.sum()) for lep, col in (("ele", "n_electrons"), ("mu", "n_muons"))}


for cell in sys.argv[1:]:
    ma, era = cell.split(":")
    m = int(ma.split("_M")[1]); mva = "MVA_Score_%s" % ma; cut = cuts[m]
    cols = ["run", "luminosityBlock", "event", mva, "weight", "n_electrons", "n_muons",
            "Z_lead_lepton_pt", "Z_sublead_lepton_pt"]
    nom = uproot.open(os.path.join(SCORED, ma, "%s.root" % era))["inclusive"].arrays(cols, library="pd")
    yi, yh = yields(nom, mva, cut, False), yields(nom, mva, cut, True)
    print("=" * 90); print("CELL %s %s cut=%.3f  nominal INC ele/mu %.4f/%.4f  HASH ele/mu %.4f/%.4f"
                           % (ma, era, cut, yi["ele"], yi["mu"], yh["ele"], yh["mu"]))
    print("%-22s %9s %9s | %9s %9s   (ratio syst/nominal - 1, in %%)" % ("corr", "INC ele", "INC mu", "HASH ele", "HASH mu"))
    store = {}
    for s in SYSTS:
        for ud in ("up", "down"):
            corr = "%s_%s" % (s, ud)
            d = uproot.open(os.path.join(SCORED, "%s_%s" % (ma, corr), "%s.root" % era))["inclusive"].arrays(cols, library="pd")
            store[corr] = d
            a, b = yields(d, mva, cut, False), yields(d, mva, cut, True)
            print("%-22s %+9.3f %+9.3f | %+9.3f %+9.3f" % (corr, 100 * (a["ele"] / yi["ele"] - 1), 100 * (a["mu"] / yi["mu"] - 1),
                                                         100 * (b["ele"] / yh["ele"] - 1), 100 * (b["mu"] / yh["mu"] - 1)))
    for ud in ("up", "down"):
        a, b = store["Electron_scale_%s" % ud], store["Electron_smear_%s" % ud]
        same = len(a) == len(b) and np.array_equal(a.Z_lead_lepton_pt.to_numpy(), b.Z_lead_lepton_pt.to_numpy())
        ae = a[a.n_electrons == 2].set_index(["run", "luminosityBlock", "event"])
        ne = nom[nom.n_electrons == 2].set_index(["run", "luminosityBlock", "event"])
        j = ae.join(ne, rsuffix="_n", how="inner")
        rel = (j.Z_lead_lepton_pt / j.Z_lead_lepton_pt_n - 1).abs()
        print("Electron_scale_%s vs Electron_smear_%s identical lead-lep pT arrays: %s ; |scale/nom-1| ele median %.2e max %.2e"
              % (ud, ud, same, rel.median(), rel.max()))
    sys.stdout.flush()
