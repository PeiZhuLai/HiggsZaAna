#!/usr/bin/env python3
"""Build the background RooWorkspace (discrete-profiling envelope) for combine.

Reads the H_mass histogram from prep_bkg_data.py, fits a small family of smooth
functions (Exponential, power law, Bernstein 2 & 3), and packages them into a
RooMultiPdf + a discrete index nuisance `pdfindex_ch1`, the floating
normalisation `roomultipdf_norm`, and the binned `data_obs`.

This mirrors the flashgg background fTest / discrete-profiling step but for the
single inclusive merged category.

MUST run inside the flashgg combine environment (cmsenv), because RooMultiPdf
lives in libHiggsAnalysisCombinedLimit.
"""
import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "output"))
    args = ap.parse_args()

    hist_path = os.path.join(args.outdir, "bkg_hist.root")
    fin = ROOT.TFile.Open(hist_path)
    if not fin or fin.IsZombie():
        raise SystemExit(f"[bkg_ws] cannot open {hist_path}")
    hist = fin.Get("data_obs")
    lo = hist.GetXaxis().GetXmin()
    hi = hist.GetXaxis().GetXmax()
    nbkg = hist.Integral()
    print(f"[bkg_ws] data_obs integral={nbkg:.1f}  range=[{lo},{hi}]")

    m = ROOT.RooRealVar("H_mass", "m_{ll#gamma} [GeV]", lo, hi)
    data_obs = ROOT.RooDataHist("data_obs", "data_obs", ROOT.RooArgList(m), hist)

    pdfs = ROOT.RooArgList()
    keep = []

    # 1) exponential
    tau = ROOT.RooRealVar("bkg_exp_tau", "tau", -0.02, -1.0, 0.0)
    exp = ROOT.RooExponential("bkg_exp", "exp", m, tau)
    pdfs.add(exp)
    keep += [tau, exp]

    # 2) power law
    p1 = ROOT.RooRealVar("bkg_pow_p", "p", -3.0, -10.0, 0.0)
    pw = ROOT.RooGenericPdf("bkg_pow", "TMath::Power(@0,@1)",
                            ROOT.RooArgList(m, p1))
    pdfs.add(pw)
    keep += [p1, pw]

    # 3) Bernstein 2
    b20 = ROOT.RooRealVar("bkg_b2_0", "", 0.5, 0.0, 15.0)
    b21 = ROOT.RooRealVar("bkg_b2_1", "", 0.3, 0.0, 15.0)
    b22 = ROOT.RooRealVar("bkg_b2_2", "", 0.2, 0.0, 15.0)
    ber2 = ROOT.RooBernstein("bkg_bern2", "bern2", m, ROOT.RooArgList(b20, b21, b22))
    pdfs.add(ber2)
    keep += [b20, b21, b22, ber2]

    # 4) Bernstein 3
    b30 = ROOT.RooRealVar("bkg_b3_0", "", 0.5, 0.0, 15.0)
    b31 = ROOT.RooRealVar("bkg_b3_1", "", 0.3, 0.0, 15.0)
    b32 = ROOT.RooRealVar("bkg_b3_2", "", 0.2, 0.0, 15.0)
    b33 = ROOT.RooRealVar("bkg_b3_3", "", 0.1, 0.0, 15.0)
    ber3 = ROOT.RooBernstein("bkg_bern3", "bern3", m,
                             ROOT.RooArgList(b30, b31, b32, b33))
    pdfs.add(ber3)
    keep += [b30, b31, b32, b33, ber3]

    # pre-fit each so the multipdf starts from sane parameters
    for i in range(pdfs.getSize()):
        pdfs.at(i).fitTo(data_obs, ROOT.RooFit.PrintLevel(-1),
                         ROOT.RooFit.SumW2Error(True))

    cat = ROOT.RooCategory("pdfindex_ch1", "pdfindex_ch1")
    multipdf = ROOT.RooMultiPdf("roomultipdf", "roomultipdf", cat, pdfs)
    norm = ROOT.RooRealVar("roomultipdf_norm", "bkg norm",
                           nbkg, 0.0, 3.0 * max(nbkg, 1.0))

    ws = ROOT.RooWorkspace("wbkg", "wbkg")
    getattr(ws, "import")(data_obs)
    getattr(ws, "import")(cat)
    getattr(ws, "import")(norm)
    getattr(ws, "import")(multipdf)
    out = os.path.join(args.outdir, "bkg_ws.root")
    ws.writeToFile(out)
    print(f"[bkg_ws] wrote {out}  (RooMultiPdf with {pdfs.getSize()} families)")


if __name__ == "__main__":
    main()
