# Expected median limit, wilson vs xsec card, all 30 points, with file mtimes
import ROOT, sys, os, time
C = sys.argv[1]
def med(p):
    f = ROOT.TFile.Open(p); t = f.Get("limit"); v = None
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for m in range(1, 31):
    p = {k: "%s/output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (C, k, m) for k in ("xsec", "wilson")}
    a, b = med(p["xsec"]), med(p["wilson"])
    print("mA%-2d xsec %.4f wilson %.4f ratio %.3f  mtime xsec %s wilson %s" % (m, a, b, b / a, time.strftime("%m-%d %H:%M", time.localtime(os.path.getmtime(p["xsec"]))), time.strftime("%m-%d %H:%M", time.localtime(os.path.getmtime(p["wilson"])))))
