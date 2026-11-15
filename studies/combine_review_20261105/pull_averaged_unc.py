#!/usr/bin/env python3
"""Mean pull of every truth function of the bias scan with two uncertainty definitions.

Input : <bias>/mA_<m>/merged/BiasFits/biasStudy_<pdf>_fits.root (MultiDimFit --algo singles, three
        entries per toy: best fit, lower and upper edge of the 68% interval)
        the injected signal strength (mu_injected) and the function list from envelope_diagnostics.csv
Output: <W>/pull_averaged_unc.txt
For each function the pull distribution is filled as in RunBiasStudy.py, once with the uncertainty
on the side facing the injected value (the AN definition, to check against the CSV) and once with
the average of the two uncertainties, 0.5*(hi-lo); the mean is from a Gaussian fit in both cases.
Run inside cmsenv.
"""
import csv, os
import ROOT

ROOT.gROOT.SetBatch(True)
W = os.path.dirname(os.path.abspath(__file__))
BIAS = ("/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit/"
        "Combine/Checks/Bias_nominal/bias_outputs_highstat")
CSV = "/eos/home-p/pelai/HZa/output_tables/bkg_envelope/envelope_diagnostics.csv"


def gaus_mean(vals, name):
    h = ROOT.TH1F(name, "", 80, -4., 4.)
    for v in vals:
        h.Fill(v)
    h.Fit("gaus", "Q0")
    f = h.GetFunction("gaus")
    return f.GetParameter(1), f.GetParError(1)


def main():
    rows = [r for r in csv.DictReader(open(CSV)) if r["bias_mean"]]
    out = ["mA  pdf    facing(this)       facing(csv)   averaged"]
    avg_all, fac_all = [], []
    for r in rows:
        m, pdf, mu = int(r["mA"]), r["pdf"], float(r["mu_injected"])
        tf = ROOT.TFile.Open("%s/mA_%d/merged/BiasFits/biasStudy_%s_fits.root" % (BIAS, m, pdf))
        rv = [(e.r, e.quantileExpected) for e in tf.Get("limit")]
        tf.Close()
        fac, avg = [], []
        for i in range(len(rv) // 3):
            (bf, q0), (lo, q1), (hi, q2) = rv[3 * i:3 * i + 3]
            if q0 != -1 or abs(q1 + 0.32) > 0.001 or abs(q2 - 0.32) > 0.001:
                continue
            u = (hi - bf) if bf < mu else (bf - lo)
            if u > 0:
                fac.append((bf - mu) / u)
            a = 0.5 * (hi - lo)
            if a > 0:
                avg.append((bf - mu) / a)
        mf, ef = gaus_mean(fac, "f_%d_%s" % (m, pdf))
        ma, ea = gaus_mean(avg, "a_%d_%s" % (m, pdf))
        fac_all.append(mf); avg_all.append((ma, m, pdf))
        out.append("%-3d %-6s %+.3f+-%.3f     %+.3f        %+.3f+-%.3f" % (m, pdf, mf, ef, float(r["bias_mean"]), ma, ea))
    worst = max(avg_all, key=lambda x: abs(x[0]))
    out.append("functions: %d; mean over functions: facing %+.3f, averaged %+.3f; largest |averaged|: %+.3f (%s, mA%d)"
               % (len(rows), sum(fac_all) / len(fac_all), sum(a[0] for a in avg_all) / len(avg_all), worst[0], worst[2], worst[1]))
    open(W + "/pull_averaged_unc.txt", "w").write("\n".join(out) + "\n")
    print(out[-1])


if __name__ == "__main__":
    main()
