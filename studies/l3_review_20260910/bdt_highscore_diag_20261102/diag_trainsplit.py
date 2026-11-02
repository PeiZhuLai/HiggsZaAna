#!/usr/bin/env python3
"""Is the high-score data excess of AN Fig. 104 (m_a = 3, 4, 5 GeV) caused by BDT training events
in the background MC? The trainers use the 'train' tree of run3_bdt_inputs_fsrfix/All_Bkg (50% of
the events); 'validation' + 'test' (50%) are never trained on. Every scored background event is
matched to its split by (event, H_m, Z_m), and the excess at score > 0.5 (MC scaled to data at
score < 0.5, as in the AN plots) is compared for MC = all events, train split only, non-train only.
Read-only; writes diag_trainsplit.txt / .json next to this script.
"""
import json, math, os
import numpy as np
import uproot
from scipy import stats

B = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
ALLBKG = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix/All_Bkg/run3.root"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
BKG = {e: ["DYGto2LG_10to100", "DYJetsToLL"] for e in ERAS[:4]}
BKG["2024"] = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
MASSES = [3, 4, 5]
EDGES = np.linspace(0, 1, 21); C = 0.5 * (EDGES[1:] + EDGES[:-1]); HI = 0.5
OUT = os.path.dirname(os.path.abspath(__file__))


def key(ev, hm, zm):
    return [(int(a), round(float(b), 3), round(float(c), 3)) for a, b, c in zip(ev, hm, zm)]


split_of = {}
f = uproot.open(ALLBKG)
for tr in ("train", "validation", "test"):
    a = f[tr].arrays(["event", "H_m", "Z_m"], library="np")
    for k in key(a["event"], a["H_m"], a["Z_m"]):
        split_of[k] = tr
print("split map:", {t: sum(1 for v in split_of.values() if v == t) for t in ("train", "validation", "test")}, flush=True)

brs = ["H_m", "Z_m", "event", "z_mumu", "z_ee", "weight_sideband_rwgt"] + ["MVA_Score_mA_M%d" % m for m in MASSES]
mc = []; data = []; unmatched = 0; total = 0
for e in ERAS:
    d = uproot.open("%s/Data/%s.root" % (B, e))["inclusive"].arrays(["H_m", "z_mumu", "z_ee"] + brs[6:], library="np")
    data.append({k: v[(d["H_m"] > 95) & (d["H_m"] < 180)] for k, v in d.items()})
    for s in BKG[e]:
        a = uproot.open("%s/%s/%s.root" % (B, s, e))["inclusive"].arrays(brs, library="np")
        m = (a["H_m"] > 95) & (a["H_m"] < 180)
        a = {k: v[m] for k, v in a.items()}
        sp = np.array([split_of.get(k, "none") for k in key(a["event"], a["H_m"], a["Z_m"])])
        unmatched += int((sp == "none").sum()); total += len(sp)
        a["split"] = sp
        mc.append(a)
data = {k: np.concatenate([x[k] for x in data]) for k in data[0]}
mc = {k: np.concatenate([x[k] for x in mc]) for k in mc[0]}
print("background events %d, unmatched %d" % (total, unmatched), flush=True)


def h(a, mA, w, sel):
    s = a["MVA_Score_mA_M%d" % mA]
    return (np.histogram(s[sel], EDGES, weights=w[sel])[0], np.histogram(s[sel], EDGES, weights=w[sel] ** 2)[0])


def chan(a, ch):
    if ch == "ee":
        return np.abs(a["z_mumu"]) != 1
    if ch == "mumu":
        return np.abs(a["z_ee"]) != 1
    return np.ones(len(a["H_m"]), bool)


res = {}; L = ["data/MC at score>0.5 after scaling MC to data at score<0.5; chi2 over score>0.5 bins",
               "MC = all (as in the AN), train split only, non-train (validation+test) only"]
for mA in MASSES:
    for region in ("incl", "sb"):
        for ch in ("all", "ee", "mumu"):
            dsel = chan(data, ch)
            if region == "sb":
                dsel &= (data["H_m"] < 115) | (data["H_m"] > 135)
            else:
                dsel &= data["MVA_Score_mA_M%d" % mA] < 0.9
            d, _ = h(data, mA, np.ones(len(data["H_m"])), dsel)
            top = 0.9 if region == "incl" else 1.0
            hi = (C > HI) & (C < top); lo = C < HI
            row = {}
            for tag, ss in (("all", None), ("train", ("train",)), ("nontrain", ("validation", "test"))):
                msel = chan(mc, ch)
                if region == "sb":
                    msel &= (mc["H_m"] < 115) | (mc["H_m"] > 135)
                if ss:
                    msel &= np.isin(mc["split"], ss)
                m, m2 = h(mc, mA, mc["weight_sideband_rwgt"].astype(float), msel)
                n = d[lo].sum() / m[lo].sum(); m = m * n; m2 = m2 * n * n
                ex = d[hi].sum() / m[hi].sum()
                err = ex * math.hypot(math.sqrt(m2[hi].sum()) / m[hi].sum(), 1 / math.sqrt(d[hi].sum()))
                v = np.where(d > 0, d, m) + m2; ok = hi & (v > 0)
                c2 = float(np.sum((d[ok] - m[ok]) ** 2 / v[ok])); nd = int(ok.sum())
                row[tag] = {"norm_lo": float(n), "excess_hi": float(ex), "err": float(err), "chi2": c2, "ndf": nd,
                            "p": float(stats.chi2.sf(c2, nd))}
            res["M%d/%s/%s" % (mA, region, ch)] = row
            L.append("mA%d %-4s %-4s | all %5.3f+-%.3f (p %.2g) | train %5.3f+-%.3f (p %.2g) | non-train %5.3f+-%.3f (p %.2g)" % (
                mA, region, ch, row["all"]["excess_hi"], row["all"]["err"], row["all"]["p"],
                row["train"]["excess_hi"], row["train"]["err"], row["train"]["p"],
                row["nontrain"]["excess_hi"], row["nontrain"]["err"], row["nontrain"]["p"]))
json.dump({"unmatched": unmatched, "total": total, "res": res}, open(os.path.join(OUT, "diag_trainsplit.json"), "w"), indent=1)
L.insert(0, "background events %d, unmatched to a split %d" % (total, unmatched))
open(os.path.join(OUT, "diag_trainsplit.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
