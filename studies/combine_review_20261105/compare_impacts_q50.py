#!/usr/bin/env python3
"""Expected impacts on the Asimov dataset with r = 1 versus r = q50 (median expected limit).

Input : Combine/output_impacts{,_xsec}/<m>_impacts.json      (r = 1)
        Combine/output_impacts{,_xsec}_q50/<m>_impacts.json  (r = q50)
Output: <W>/compare_impacts_q50.txt
The impact of a parameter is the larger of |r(+1 sigma) - r_hat| and |r(-1 sigma) - r_hat|. Since r_hat
differs by a factor ~50 between the two sets, impacts are quoted relative to the total uncertainty on r
(half the 68% interval), which makes them comparable. Per card and mass point: the five leading
parameters of each set, and the summed share of the signal-yield (lnN and signal-shape) nuisances
versus the background-model parameters among the ten leading ones.
"""
import json, os

W = os.path.dirname(os.path.abspath(__file__))
C = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/Combine"
BKG_PREFIX = ("Exp_", "Pow_", "Lau_", "Bern_", "shapeBkg_", "env_pdf")


def load(path):
    d = json.load(open(path))
    r = [p for p in d["POIs"] if p["name"] == "r"][0]["fit"]
    sig = 0.5 * (r[2] - r[0])
    imp = []
    for p in d["params"]:
        lo, nom, hi = p["r"]
        imp.append((max(abs(hi - nom), abs(lo - nom)) / sig, p["name"]))
    imp.sort(reverse=True)
    return r, imp


def is_bkg(n):
    return n.startswith(BKG_PREFIX)


def main():
    out = []
    for card, base in (("wilson", "output_impacts"), ("xsec", "output_impacts_xsec")):
        out.append("=== %s card: impact / sigma_r, five leading parameters" % card)
        for m in range(1, 31):
            a, b = "%s/%s/%d_impacts.json" % (C, base, m), "%s/%s_q50/%d_impacts.json" % (C, base, m)
            if not (os.path.exists(a) and os.path.exists(b)):
                out.append("mA%-2d missing" % m); continue
            ra, ia = load(a); rb, ib = load(b)
            fa = sum(x for x, n in ia[:10] if is_bkg(n)) / sum(x for x, _ in ia[:10])
            fb = sum(x for x, n in ib[:10] if is_bkg(n)) / sum(x for x, _ in ib[:10])
            out.append("mA%-2d r=1   (%.3f %+.3f/-%.3f) bkg share of top10 %.2f: %s" % (
                m, ra[1], ra[2] - ra[1], ra[1] - ra[0], fa, ", ".join("%s %.3f" % (n, x) for x, n in ia[:5])))
            out.append("     r=q50 (%.4f %+.4f/-%.4f) bkg share of top10 %.2f: %s" % (
                rb[1], rb[2] - rb[1], rb[1] - rb[0], fb, ", ".join("%s %.3f" % (n, x) for x, n in ib[:5])))
    open(W + "/compare_impacts_q50.txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
