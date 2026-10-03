"""Per-era conversion-safe electron veto lnN range over the 30 mass points (Sec 9.2.1 text).
usage: csev_per_era.py <datacard_dir>"""
import sys, re
D = sys.argv[1]; R = {}
for m in range(1, 31):
    procs = None
    for line in open("%s/%d_pruned_datacard_leptons.txt" % (D, m)):
        w = line.split()
        if w and w[0] == "process" and procs is None: procs = w[1:]
        if len(w) > 2 and w[0].startswith("CMS_hza_photon_csev_") and w[1] == "lnN":
            era = w[0].replace("CMS_hza_photon_csev_", "")
            for p, t in zip(procs, w[2:]):
                if t == "-": continue
                if "/" in t: dn, up = map(float, t.split("/")); d = 0.5 * (abs(up - 1) + abs(dn - 1))
                else: d = abs(float(t) - 1)
                R.setdefault(era, []).append(100 * d)
for e, v in R.items(): print("%-14s %.2f -- %.2f %%" % (e, min(v), max(v)))
