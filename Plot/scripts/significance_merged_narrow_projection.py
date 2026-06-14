#!/usr/bin/env python3
"""
PROJECTED merged-category significance when a narrow window on the regressed
merged-photon mass (MergedML_mass) is added on top of the m_H window.

** This is a projection, not a measurement. **  No background sample carrying
MergedML_mass exists yet (the MLNANO merged-DNA production was run on signal
only; merged bkg/data are defined only for 2023postBPix and not produced).
We therefore estimate the background rejection of a +/-2 sigma_M window under a
stated assumption and propagate it onto the measured broad-window yields.

Signal side (rigorous, measured here):
  sigma_M(m_a) and the +/-2 sigma signal efficiency are taken from the signal
  MergedML_mass distribution (parquet_merged_DNA_tmp, sub-GeV + M1, 2024).

Background side (assumption, flagged):
  background MergedML_mass is taken uniform over the regression range [0, 1.2] GeV,
  so a +/-2 sigma_M window keeps f_B = min(1, 4 sigma_M / 1.2) of the background.
  -> optimistic gain is bounded because sigma_M is large vs the 1.2 GeV range.

Broad-window yields S_mrg, B_mrg come from
  Plot/output/significance_merged_vs_resolved.json  (m_H window only).
"""
import glob
import json
import math
import numpy as np
import pyarrow.parquet as pq

EOS = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
REG_RANGE = 1.2          # MergedML regression output range [0, 1.2] GeV
BASES = [f"{EOS}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_all",
         f"{EOS}/parquet_merged_DNA_tmp/Sig_MC_MLNANO_M1"]
BROAD = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/significance_merged_vs_resolved.json"


def asimov(s, b):
    if s <= 0 or b <= 0:
        return 0.0
    return math.sqrt(max(0.0, 2 * ((s + b) * math.log(1 + s / b) - s)))


def sig_path(tag):
    for b in BASES:
        p = f"{b}/mA_MLNANO_{tag}_2024/merged_nominal.parquet"
        if glob.glob(p):
            return p
    return None


def resolution(tag):
    p = sig_path(tag)
    if not p:
        return None
    m = pq.read_table(p, columns=["MergedML_mass"]).to_pandas()["MergedML_mass"].to_numpy()
    v = m[np.isfinite(m) & (m > 0) & (m < REG_RANGE)]
    if len(v) < 50:
        return None
    peak = np.median(v)
    sigma = (np.percentile(v, 84) - np.percentile(v, 16)) / 2.0
    lo, hi = peak - 2 * sigma, peak + 2 * sigma
    eff = float(np.mean((v > lo) & (v < hi)))      # signal kept by +/-2 sigma window
    return dict(peak=peak, sigma=sigma, effS=eff,
                window=[max(0.0, lo), min(REG_RANGE, hi)])


broad = {r["ma"]: r for r in json.load(open(BROAD))}
TAGS = [("M0p%d" % i, i / 10.0) for i in range(1, 10)] + [("M1", 1.0)]

print(f"{'m_a':>5} {'sigma_M':>8} {'sig/mA':>7} {'effS':>6} {'f_B':>6} "
      f"{'Z_broad':>8} {'Z_narrow':>9} {'gain':>6}")
print("-" * 70)
out = []
for tag, ma in TAGS:
    r = resolution(tag)
    b = broad.get(ma) or broad.get(float(ma))
    if r is None or b is None or b.get("Z_merged") is None:
        print(f"{ma:>5.1f}  (missing)")
        continue
    win_w = r["window"][1] - r["window"][0]
    f_B = min(1.0, win_w / REG_RANGE)              # uniform-background assumption
    S0, B0, Z0 = b["S_mrg"], b["B_mrg"], b["Z_merged"]
    S1 = S0 * r["effS"]
    B1 = B0 * f_B
    Z1 = asimov(S1, B1)
    gain = Z1 / Z0 if Z0 > 0 else float("nan")
    print(f"{ma:>5.1f} {r['sigma']:>8.3f} {r['sigma']/ma:>7.2f} {r['effS']:>6.2f} "
          f"{f_B:>6.2f} {Z0:>8.3f} {Z1:>9.3f} {gain:>6.2f}")
    out.append(dict(ma=ma, sigma_M=r["sigma"], effS=r["effS"], f_B=f_B,
                    Z_broad=Z0, Z_narrow=Z1, gain=gain,
                    S_narrow=S1, B_narrow=B1))

js = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/significance_merged_narrow_projection.json"
json.dump(out, open(js, "w"), indent=2)
print("\nwrote", js)
print("ASSUMPTION: background MergedML_mass uniform over [0,1.2] GeV (flagged projection).")
print("Real number requires producing the 2024 merged-DNA background / data with MergedML_mass.")
