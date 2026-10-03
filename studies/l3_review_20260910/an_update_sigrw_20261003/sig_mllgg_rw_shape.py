#!/usr/bin/env python3
"""Signal m_llgg before/after the sideband reweight, after the BDT working point (nominal trees of the
MVA-cut outputs; 'before' = weight / weight_sideband_rwgt). Prints, per mass point and channel summed over
the five eras: ratio after/before in 5-GeV bins within 110-140 GeV (each normalized to its own total),
and the shift of the mean and of the 68% interval width."""
import uproot, numpy as np
B = "/eos/home-p/pelai/HZa/root_MVAcut/sig"
ERAS = ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")
edges = np.arange(110, 141, 5)
def width68(x, w):
    o = np.argsort(x); x, w = x[o], np.cumsum(w[o]) / w.sum()
    lo, hi = np.interp([0.16, 0.84], w, x); return hi - lo
print("mA ch   yield_ratio  mean_shift[GeV] width68_rel  ratio(5GeV bins 110-140, shape only)")
for m in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30):
    for ch in ("ele", "mu"):
        xs, wa, wb = [], [], []
        for y in ERAS:
            t = uproot.open("%s/mA_M%d/output_%s.root" % (B, m, y))["DiphotonTree/ggh_125_Za_%s_13p6TeV_cat0" % ch]
            a = t.arrays(["CMS_hza_mass", "weight", "weight_sideband_rwgt"], library="np")
            xs.append(a["CMS_hza_mass"]); wa.append(a["weight"]); wb.append(a["weight"] / a["weight_sideband_rwgt"])
        x, wa, wb = np.concatenate(xs), np.concatenate(wa), np.concatenate(wb)
        ha, _ = np.histogram(x, edges, weights=wa); hb, _ = np.histogram(x, edges, weights=wb)
        r = (ha / ha.sum()) / (hb / hb.sum())
        ma, mb = np.average(x, weights=wa), np.average(x, weights=wb)
        print("%2d %-4s %8.3f %12.3f %11.3f   %s" % (m, ch, wa.sum() / wb.sum(), ma - mb,
              width68(x, wa) / width68(x, wb) - 1, " ".join("%.3f" % v for v in r)))
