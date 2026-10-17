#!/usr/bin/env python3
"""BDT cuts of a scored copy that give target data-sideband counts (95-115 or 135-180 GeV, the count definition of
the variable search): n_production(m_a) x factor. Prints one line per m_a: "<mA> <cut> <cut> ...".
Usage: cuts_for_counts.py <scored base> <mA,...> <factor,...>
"""
import json
import sys

import numpy as np
import uproot

base, masses, factors = sys.argv[1], [int(x) for x in sys.argv[2].split(",")], [float(x) for x in sys.argv[3].split(",")]
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
nprod = json.load(open("/eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples/meta.json"))["prod_data_sb_count"]
for m in masses:
    br = "MVA_Score_mA_M%d" % m
    s = []
    for e in ERAS:
        a = uproot.open("%s/Data/%s.root" % (base, e))["inclusive"].arrays(["H_mass", br], library="np")
        h = a["H_mass"]; sb = ((h > 95) & (h < 115)) | ((h > 135) & (h < 180))
        s.append(a[br][sb])
    s = np.sort(np.concatenate(s)); s = s[np.isfinite(s)][::-1]
    cuts = []
    for f in factors:
        n = int(round(nprod[str(m)] * f))
        cuts.append("%.5f" % (0.5 * (s[n - 1] + s[n])))
    print(m, " ".join(cuts))
