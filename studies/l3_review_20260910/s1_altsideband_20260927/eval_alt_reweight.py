#!/usr/bin/env python3
"""L3 review S1: effect of alternative sideband reweightings on the signal after the BDT WP.

Mirrors how the reweighting enters the signal in the current chain:
  * the BDT (trained with the NOMINAL h_m reweighting, signal weighted with the TRUE-mass
    param) is NOT retrained; the stored scores of run3_bdt_scored_fsrfix are used;
  * the signal sample is the 'test' tree (30% split) as in apply_bdt_sig.py, the WP is
    MVA_Score_mA_M<m> > MVAcut (Plot/output/MVAcut_points_run3.json), channels n_electrons==2 /
    n_muons==2, yield weight 'weight', fit variable H_mass (renamed CMS_hza_mass downstream);
  * the reweighting R enters as a per-event weight on the signal. Its overall normalization on
    the signal is arbitrary (shape-only correction), so the yield is quoted as the BDT-WP
    efficiency eff_R = sum_{pass & fit window} w R / sum_{all preselected} w R, i.e. R is
    normalized to preserve the preselected signal yield per channel.
Variants: none (R=1, = what the datacard yield uses), nominal (h_m), and the alternatives.
The param step uses the signal TRUE mass, (ALP_m - mA)/H_m, as in the training
(hza_features.signal_weight); the event-hash definition used by apply_bdt_sig.py for the
current nuisance is evaluated as a cross-check (suffix _hash).
Also reads the currently stored nuisance weight_mva_reweight_Up/Down from root_MVAcut/sig.
"""
import json, os, sys, math, argparse
import numpy as np, pandas as pd, uproot

sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from sideband_reweight import SidebandReweighter  # noqa: E402

RW = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights"
SCORED = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
MVACUT = "/eos/home-p/pelai/HZa/root_MVAcut/sig"
WPJSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
FIT_WINDOW = {1: (115.0, 133.0)}     # signalFit.py _resolve_fit_range; default 100-180
DEFAULT_WINDOW = (100.0, 180.0)
NBOOT = 100

# S1_RW_TAG selects the reweight set: "fsrfix" (2026-09-27 study) or "fsrfix_dyveto" (after the
# 2024 DY+jets overlap-veto retrain). Variants whose JSON does not exist are skipped.
_T = os.environ.get("S1_RW_TAG", "fsrfix")
VARIANTS = {
    "nominal": f"{RW}/sideband_run3_iterative_{_T}.json",
    "zm": f"{RW}/sideband_run3_iterative_{_T}_zm.json",
    "zmnarrow": f"{RW}/sideband_run3_iterative_{_T}_zmnarrow.json",
    "zmhsb": f"{RW}/sideband_run3_iterative_{_T}_zmhsb.json",
    "zmnarrowhsb": f"{RW}/sideband_run3_iterative_{_T}_zmnarrowhsb.json",
}

READ = ["var_dR_Za", "var_dR_g1g2", "var_dR_g1Z", "pho1IetaIeta55", "pho2IetaIeta55",
        "pho1PIso_noCorr", "pho2PIso_noCorr", "ALP_calculatedPhotonIso", "pho1R9", "pho2R9",
        "var_PtaOverMh", "pho1Pt_oHm", "pho2Pt_oHm", "H_pt_oHm", "Z_m", "H_m", "ALP_m",
        "H_mass", "event", "weight", "n_electrons", "n_muons"]


def wp_cuts():
    d = json.load(open(WPJSON))
    return {int(r["mA"]): float(r["MVAcut"]) for r in d["results"]}


def sigma_eff(x, w, frac=0.683):
    if len(x) < 5 or np.sum(w) <= 0:
        return float("nan")
    o = np.argsort(x); x = x[o]; c = np.cumsum(w[o]); c /= c[-1]
    best = np.inf
    j = 0
    for i in range(len(x)):
        target = (c[i - 1] if i > 0 else 0.0) + frac
        j = max(j, i)
        while j < len(x) and c[j] < target:
            j += 1
        if j >= len(x):
            break
        best = min(best, x[j] - x[i])
    return 0.5 * best


def wstats(x, w):
    sw = np.sum(w)
    if sw <= 0:
        return float("nan"), float("nan")
    m = np.sum(w * x) / sw
    return float(m), float(math.sqrt(max(np.sum(w * (x - m) ** 2) / sw, 0.0)))


def stored_nuisance(m, era):
    out = {}
    p = f"{MVACUT}/mA_M{m}/output_{era}.root"
    if not os.path.exists(p):
        return out
    with uproot.open(p) as f:
        for lep in ("ele", "mu"):
            k = f"DiphotonTree/ggh_125_Za_{lep}_13p6TeV_cat0"
            if k not in f:
                continue
            b = f[k].arrays(["weight", "weight_mva_reweight_Up", "weight_mva_reweight_Down"], library="np")
            if len(b["weight"]) == 0:
                continue
            out[lep] = (float(b["weight_mva_reweight_Up"][0]), float(b["weight_mva_reweight_Down"][0]),
                        float(np.sum(b["weight"])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masses", default=",".join(map(str, MASSES)))
    ap.add_argument("--eras", default=",".join(ERAS))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    masses = [int(x) for x in args.masses.split(",")]
    eras = args.eras.split(",")
    cuts = wp_cuts()
    rws = {k: SidebandReweighter.from_json(v) for k, v in VARIANTS.items() if os.path.exists(v)}
    print("variants:", list(rws), flush=True)
    rows = []
    rng = np.random.default_rng(20260927)
    for m in masses:
        lo, hi = FIT_WINDOW.get(m, DEFAULT_WINDOW)
        cut = cuts[m]
        for era in eras:
            p = f"{SCORED}/mA_M{m}/{era}.root"
            with uproot.open(p) as f:
                a = f["test"].arrays(READ + [f"MVA_Score_mA_M{m}"], library="pd")
            a = a.rename(columns={"pho1PIso_noCorr": "pho1ECALIso", "pho2PIso_noCorr": "pho2ECALIso"})
            w = a["weight"].to_numpy(float)
            score = a[f"MVA_Score_mA_M{m}"].to_numpy(float)
            mass = a["H_mass"].to_numpy(float)
            passed = score > cut
            inwin = (mass >= lo) & (mass <= hi)
            R = {"none": np.ones(len(a))}
            for k, rw in rws.items():
                fr = a.copy()
                fr["param"] = (fr["ALP_m"].to_numpy(float) - m) / fr["H_m"].to_numpy(float)
                R[k] = rw.weights_for_dataframe(fr)
                fr2 = a.copy()          # event-hash param (as apply_bdt_sig.py)
                R[k + "_hash"] = rw.weights_for_dataframe(fr2)
            stored = stored_nuisance(m, era)
            nev = len(a)
            boot = rng.poisson(1.0, size=(NBOOT, nev)).astype(float)
            for ch, chm in (("ele", a["n_electrons"].to_numpy() == 2),
                            ("mu", a["n_muons"].to_numpy() == 2),
                            ("all", (a["n_electrons"].to_numpy() == 2) | (a["n_muons"].to_numpy() == 2))):
                sel = chm & passed & inwin
                row = {"mA": m, "era": era, "ch": ch, "wp": cut, "win_lo": lo, "win_hi": hi,
                       "n_presel": int(chm.sum()), "n_pass": int(sel.sum()),
                       "yield_none": float(np.sum(w[sel]))}
                effs = {}
                beffs = {}
                for k, r in R.items():
                    wr = w * r
                    den = np.sum(wr[chm]); num = np.sum(wr[sel])
                    effs[k] = num / den if den > 0 else float("nan")
                    bnum = boot[:, sel] @ wr[sel]; bden = boot[:, chm] @ wr[chm]
                    beffs[k] = bnum / bden
                    mu_, sd_ = wstats(mass[sel], wr[sel])
                    row[f"eff_{k}"] = effs[k]
                    row[f"mean_{k}"] = mu_
                    row[f"rms_{k}"] = sd_
                    row[f"seff_{k}"] = sigma_eff(mass[sel], wr[sel])
                    row[f"rmean_R_{k}"] = float(np.sum(wr[chm]) / np.sum(w[chm])) if np.sum(w[chm]) > 0 else float("nan")
                for k in R:
                    if k in ("none",) or k.startswith("nominal"):
                        continue
                    ref = "nominal_hash" if k.endswith("_hash") else "nominal"
                    if ref not in effs:
                        continue
                    row[f"dy_{k}"] = effs[k] / effs[ref] - 1.0
                    row[f"dyerr_{k}"] = float(np.std(beffs[k] / beffs[ref] - 1.0))
                for k in ("nominal", "nominal_hash"):
                    if k in effs:
                        row[f"dy_{k}_vs_none"] = effs[k] / effs["none"] - 1.0
                        row[f"dyerr_{k}_vs_none"] = float(np.std(beffs[k] / beffs["none"] - 1.0))
                if ch in stored:
                    up, dn, _ = stored[ch]
                    row["cur_up"] = up - 1.0; row["cur_dn"] = dn - 1.0
                elif ch == "all" and stored:
                    tot = sum(v[2] for v in stored.values())
                    row["cur_up"] = sum((v[0] - 1.0) * v[2] for v in stored.values()) / tot
                    row["cur_dn"] = sum((v[1] - 1.0) * v[2] for v in stored.values()) / tot
                rows.append(row)
            r_all = [x for x in rows if x["mA"] == m and x["era"] == era and x["ch"] == "all"][0]
            print(f"mA{m} {era}: npass={r_all['n_pass']} " + " ".join(
                f"{k}={100*r_all.get('dy_'+k, float('nan')):+.2f}%" for k in ("zm", "zmnarrow", "zmhsb", "zmnarrowhsb"))
                + f" nom_vs_none={100*r_all.get('dy_nominal_vs_none', float('nan')):+.2f}% cur=+{100*r_all.get('cur_up', float('nan')):.2f}/{100*r_all.get('cur_dn', float('nan')):.2f}%",
                flush=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
