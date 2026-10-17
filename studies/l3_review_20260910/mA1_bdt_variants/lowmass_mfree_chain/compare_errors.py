#!/usr/bin/env python3
"""Bootstrap errors (200 replicas) on R, R_lo, R_hi at the compared working points, and raw data
sideband counts per band, to check that the production vs win_mfree differences in
compare_models.txt are not statistical. Env: higgs-alp-ana."""
import json, os, numpy as np
from compare_models import load, ERAS, BKG, inr, is_sb, frac
HERE = os.path.dirname(os.path.abspath(__file__))
S = json.load(open(os.path.join(HERE, "compare_models.json")))
dy = load([(s, e) for e in ERAS for s in BKG[e]], "inclusive")
da = load([("Data", e) for e in ERAS], "inclusive")
rng = np.random.default_rng(7)
out = []
for key, rec in S.items():
    tag, mA = key.split("_")[0], int(key.split("_")[1][2:])
    sc = f"MVA_Score_mA_M{mA}"; c = rec["cut"]
    b = dy[tag]; m, w, s = b["H_mass"], b["factor"].astype(float), b[sc]
    ok = inr(m, 95, 180) & np.isfinite(w) & np.isfinite(s); m, w, s = m[ok], w[ok], s[ok]
    p = s > c; allm = np.ones_like(m, bool)
    def ratios(ww):
        return [frac(ww, m, lo, hi, p) / frac(ww, m, lo, hi, allm) for lo, hi in ((120, 130), (105, 115), (135, 145))]
    bt = np.array([ratios(w * rng.poisson(1.0, len(w))) for _ in range(200)])
    neff = [(w[p & inr(m, lo, hi)].sum() ** 2 / (w[p & inr(m, lo, hi)] ** 2).sum()) for lo, hi in ((120, 130), (105, 115), (135, 145))]
    d = da[tag]; dm, ds = d["H_mass"], d[sc]; okd = inr(dm, 95, 180) & is_sb(dm) & np.isfinite(ds); dm, ds = dm[okd], ds[okd]
    bands = [(95, 105), (105, 115), (135, 145), (145, 160), (160, 180)]
    n_pass = [int((inr(dm, lo, hi) & (ds > c)).sum()) for lo, hi in bands]
    n_all = [int(inr(dm, lo, hi).sum()) for lo, hi in bands]
    exp = [sum(n_pass) * a / sum(n_all) for a in n_all]
    line = (f"{key:18s} cut={c:.4f}  R={rec['R']:.2f}+-{np.nanstd(bt[:,0]):.2f}  R_lo={rec['R_lo']:.2f}+-{np.nanstd(bt[:,1]):.2f}  "
            f"R_hi={rec['R_hi']:.2f}+-{np.nanstd(bt[:,2]):.2f}  DY n_eff(120-130/105-115/135-145)={neff[0]:.0f}/{neff[1]:.0f}/{neff[2]:.0f}  "
            f"data SB pass/expected-if-unsculpted: " + " ".join(f"{a}/{e:.0f}" for a, e in zip(n_pass, exp)))
    print(line, flush=True); out.append(line)
open(os.path.join(HERE, "compare_errors.txt"), "w").write("\n".join(out) + "\n")
