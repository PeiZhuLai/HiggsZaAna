#!/usr/bin/env python3
"""
Expected Asimov significance vs m_a for three strategies:
  RESOLVED         : two resolved photons, signal peak in ALP_mass (= m_a)
  MERGED           : merged-photon selection (pass_allcuts_merged_AN2020),
                     signal peak in H_mass (~125), m_a only enters via signal eff
  RESOLVED+MERGED  : quadrature sum in the m_a range where both exist (assumes
                     the two categories are orthogonal -> optimistic upper bound)

Asimov Z = sqrt( 2 * ( (S+B) ln(1+S/B) - S ) ).

Normalisation: weight_central already carries xs * lumi / sum(genW).
Signal uses the common reference cross section xs = 0.1 pb (same for every m_a),
so the three curves are directly comparable across m_a.  2024 only
(merged-category backgrounds exist only for 2024); lumi = 109.82 fb^-1.
"""
import glob
import json
import math
import numpy as np
import pyarrow.parquet as pq

EOS = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
H_LO, H_HI = 120.0, 130.0          # m_H window for the resolved baseline / merged peak search
NSIG = 2.0                          # +/- N sigma_eff mass window

RES_BKG = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]
MRG_BKG = ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"]

# m_a points
RES_MA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
MRG_MA = [("M0p%d" % i, i / 10.0) for i in range(1, 10)] + \
         [("M%d" % i, float(i)) for i in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]


def asimov(s, b):
    if s <= 0 or b <= 0:
        return 0.0
    return math.sqrt(max(0.0, 2 * ((s + b) * math.log(1 + s / b) - s)))


def load(path, cols):
    fs = glob.glob(path)
    if not fs:
        return None
    return pq.read_table(fs[0], columns=cols).to_pandas()


# ---- background caches (read each big file only once) -----------------------
_RES_BKG_CACHE = None   # list of (ALP_mass, H_mass, weight) arrays
_MRG_BKG_CACHE = None   # list of (H_mass, weight, sel) arrays


def res_bkg():
    global _RES_BKG_CACHE
    if _RES_BKG_CACHE is None:
        _RES_BKG_CACHE = []
        for bk in RES_BKG:
            d = load(f"{EOS}/parquet_DNA/Bkg_MC/{bk}_2024/merged_nominal.parquet",
                     ["ALP_mass", "H_mass", "weight_central"])
            if d is None:
                print(f"  [res bkg missing] {bk}", flush=True)
                continue
            _RES_BKG_CACHE.append((d.ALP_mass.to_numpy(), d.H_mass.to_numpy(),
                                   d.weight_central.to_numpy()))
        print(f"  [cached {len(_RES_BKG_CACHE)} resolved bkg]", flush=True)
    return _RES_BKG_CACHE


def mrg_bkg():
    global _MRG_BKG_CACHE
    if _MRG_BKG_CACHE is None:
        _MRG_BKG_CACHE = []
        for bk in MRG_BKG:
            d = load(f"{EOS}/parquet_friend/Bkg_{bk}_2024/{bk}_2024/merged_nominal.parquet",
                     ["H_mass", "weight_central", "pass_allcuts_merged_AN2020"])
            if d is None:
                print(f"  [mrg bkg missing] {bk}", flush=True)
                continue
            _MRG_BKG_CACHE.append((d.H_mass.to_numpy(), d.weight_central.to_numpy(),
                                   d.pass_allcuts_merged_AN2020.to_numpy().astype(bool)))
        print(f"  [cached {len(_MRG_BKG_CACHE)} merged bkg]", flush=True)
    return _MRG_BKG_CACHE


def eff_window(vals, w, lo, hi):
    """center=weighted-ish median and sigma_eff=(p84-p16)/2 inside [lo,hi]."""
    m = (vals > lo) & (vals < hi) & np.isfinite(vals)
    v = vals[m]
    if len(v) < 20:
        return None
    c = np.median(v)
    sig = (np.percentile(v, 84) - np.percentile(v, 16)) / 2.0
    return c, max(sig, 1e-3)


# ---------------------------------------------------------------- RESOLVED ----
def z_resolved(ma):
    sig = load(f"{EOS}/parquet_DNA/Sig_MC/mA_M{ma}_2024/merged_nominal.parquet",
               ["ALP_mass", "H_mass", "weight_central"])
    if sig is None:
        return None
    Hs = sig.H_mass.to_numpy(); As = sig.ALP_mass.to_numpy(); ws = sig.weight_central.to_numpy()
    inH = (Hs > H_LO) & (Hs < H_HI)
    win = eff_window(As[inH], ws[inH], 0.4 * ma, 1.6 * ma)
    if win is None:
        return None
    c, sg = win
    lo, hi = c - NSIG * sg, c + NSIG * sg
    S = ws[inH & (As > lo) & (As < hi)].sum()
    B = 0.0
    for Ab, Hb, wb in res_bkg():
        B += wb[(Hb > H_LO) & (Hb < H_HI) & (Ab > lo) & (Ab < hi)].sum()
    return dict(ma=ma, S=float(S), B=float(B), Z=asimov(S, B), win=[lo, hi], sigma=sg)


# ------------------------------------------------------------------ MERGED ----
def z_merged(tag, ma):
    sig = load(f"{EOS}/parquet_friend/mA_{tag}/mA_{tag}_2024/merged_nominal.parquet",
               ["H_mass", "weight_central", "pass_allcuts_merged_AN2020"])
    if sig is None:
        return None
    Hs = sig.H_mass.to_numpy(); ws = sig.weight_central.to_numpy()
    sel = sig.pass_allcuts_merged_AN2020.to_numpy().astype(bool)
    win = eff_window(Hs[sel], ws[sel], 100.0, 150.0)
    if win is None:
        return None
    c, sg = win
    lo, hi = c - NSIG * sg, c + NSIG * sg
    S = ws[sel & (Hs > lo) & (Hs < hi)].sum()
    B = 0.0
    for Hb, wb, selb in mrg_bkg():
        B += wb[selb & (Hb > lo) & (Hb < hi)].sum()
    return dict(ma=ma, S=float(S), B=float(B), Z=asimov(S, B), win=[lo, hi], sigma=sg)


def main():
    res = {r["ma"]: r for r in (z_resolved(ma) for ma in RES_MA) if r}
    mrg = {}
    for tag, ma in MRG_MA:
        r = z_merged(tag, ma)
        if r:
            mrg[ma] = r

    print(f"{'m_a':>6} {'Z_res':>8} {'Z_mrg':>8} {'Z_comb':>8}   "
          f"{'S_res':>8} {'B_res':>10} {'S_mrg':>8} {'B_mrg':>10}")
    print("-" * 80)
    out = []
    allma = sorted(set(list(res) + list(mrg)))
    for ma in allma:
        zr = res.get(ma, {}).get("Z")
        zm = mrg.get(ma, {}).get("Z")
        zc = math.sqrt(zr**2 + zm**2) if (zr and zm) else None
        sr = res.get(ma, {}).get("S"); br = res.get(ma, {}).get("B")
        sm = mrg.get(ma, {}).get("S"); bm = mrg.get(ma, {}).get("B")
        print(f"{ma:>6} {('%.3f'%zr) if zr is not None else '   -   ':>8} "
              f"{('%.3f'%zm) if zm is not None else '   -   ':>8} "
              f"{('%.3f'%zc) if zc is not None else '   -   ':>8}   "
              f"{('%.2f'%sr) if sr is not None else '-':>8} {('%.1f'%br) if br is not None else '-':>10} "
              f"{('%.2f'%sm) if sm is not None else '-':>8} {('%.1f'%bm) if bm is not None else '-':>10}")
        out.append(dict(ma=ma, Z_resolved=zr, Z_merged=zm, Z_combined=zc,
                        S_res=sr, B_res=br, S_mrg=sm, B_mrg=bm))

    js = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/significance_merged_vs_resolved.json"
    import os; os.makedirs(os.path.dirname(js), exist_ok=True)
    json.dump(out, open(js, "w"), indent=2)
    print("\nwrote", js)
    print("Signal normalised to reference xs = 0.1 pb; 2024, lumi 109.82 fb^-1.")
    print("Z_combined = sqrt(Z_res^2+Z_mrg^2) assumes orthogonal categories (upper bound).")


if __name__ == "__main__":
    main()
