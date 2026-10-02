"""V9 diagnostic (numbers only, no plots): fine scan 0.90-0.999 (step 0.001) on the stored scores,
same S/B/R definitions as scan_score_R_significance.py, to locate the true Z maximum (the plot grid
stops at 0.99) and to quote the MC statistical precision of R at the adopted cut (n_eff, bootstrap)."""
import sys, json, numpy as np
sys.argv = ["x"]
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts")
import scan_score_R_significance as S
rng = np.random.default_rng(12345)
cuts = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(S.JSON_PATH))["results"]}
out = {}
for m in [4,5,6,7,8,9,10,15,20,25,30]:
    bk = S.load_scored_bkg(m); sg = S.load_scored_sig(m)
    H, w, s = bk["H_mass"], bk["factor"], bk["s"]; pk = (H > 120) & (H < 130)
    rows = []
    for thr in np.round(np.arange(0.900, 0.9995, 0.001), 4):
        if ((s > thr) & pk).sum() < 20: continue
        rows.append((thr, S._R_direct(bk, thr), S._Z_direct(bk, sg, thr)))
    rows = np.array(rows); j = int(np.argmax(rows[:, 2]))
    c = cuts[m]; p = (s > c)
    wp_ = w[p & pk]; neff = wp_.sum()**2 / (wp_**2).sum()
    # bootstrap of R at the WP (resample all bkg events with Poisson(1) weights)
    idx = np.nonzero(p)[0]
    frac = w[pk].sum() / w.sum()
    Rb = []
    for _ in range(300):
        k = rng.poisson(1.0, size=idx.size)
        ww = w[idx] * k
        Rb.append(ww[pk[idx]].sum() / ww.sum() / frac)
    Rb = np.array(Rb)
    near = rows[np.abs(rows[:, 0] - c) <= 0.0101]
    out[m] = dict(cut=c, R=S._R_direct(bk, c), Z=S._Z_direct(bk, sg, c), nraw_pk=int((p & pk).sum()),
                  neff_pk=float(neff), negw_frac_pk=float((wp_ < 0).mean()), R_boot_sd=float(Rb.std()),
                  Zmax_cut=float(rows[j, 0]), Zmax=float(rows[j, 2]), R_at_Zmax=float(rows[j, 1]),
                  Rmin_near=float(near[:, 1].min()), Rmax_near=float(near[:, 1].max()),
                  n_near_gt1p1=int((near[:, 1] > 1.1).sum()), n_near=int(len(near)))
    o = out[m]
    print(f"mA{m:2d} cut {c:.3f} R={o['R']:.2f}±{o['R_boot_sd']:.2f} Z={o['Z']:.1f} nraw={o['nraw_pk']} neff={o['neff_pk']:.0f} negw={o['negw_frac_pk']:.2f} | "
          f"Zmax at {o['Zmax_cut']:.3f} Z={o['Zmax']:.1f} R={o['R_at_Zmax']:.2f} | R in ±0.01: {o['Rmin_near']:.2f}-{o['Rmax_near']:.2f} (>1.1: {o['n_near_gt1p1']}/{o['n_near']})")
json.dump(out, open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/v9v10_20260927/v9_fine_diag.json", "w"), indent=1)
