# Expected median limits of m_a = 18-22 before (backup) and after the m_a = 20 Powell refit
import ROOT, sys
D, C = sys.argv[1], sys.argv[2]
def med(p):
    f = ROOT.TFile.Open(p); t = f.Get("limit")
    for e in t:
        if abs(e.quantileExpected - 0.5) < 1e-3: v = e.limit
    f.Close(); return v
for k in ("xsec", "wilson"):
    for m in (18, 19, 20, 21, 22):
        n = "output_combine_results_%s/higgsCombine%d.AsymptoticLimits.mH125.38.root" % (k, m)
        a, b = med("%s/old_extract/%s" % (D, n)), med("%s/%s" % (C, n))
        print("%-6s mA%d  before %.4f  after %.4f  (%+.1f%%)" % (k, m, a, b, 100 * (b / a - 1)))
