#!/usr/bin/env python3
"""Signal modeling for the low-mA merged H->Z a (a->merged photon) analysis.

Fit observable: H_mass = m(ll + merged-photon), which peaks ~125 GeV.
Merged photons carry NO m_a resolution (study doc Part 2-3), so this is a single
inclusive merged category with ONE H_mass signal peak per m_a hypothesis.

For each sub-GeV mass point this:
  - reads the merged-selection signal friend parquet,
  - fits a double-sided Crystal Ball to the H_mass peak,
  - writes a per-mA RooWorkspace (pdf + constant normalisation) for combine,
  - records mean / sigma_eff / yield (at the common reference xs=0.1 pb),
  - draws a CMS Preliminary plot.

Run in the `higgs-alp-ana` env (PyROOT/RooFit).  weight_central already carries
xs * lumi / sum(genW) with xs = 0.1 pb -> the yield is the expected S at 0.1 pb.
"""
import argparse
import glob
import json
import os

import numpy as np
import pyarrow.parquet as pq
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

EOS_DEFAULT = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"


def mass_points():
    # mA = 0.1 .. 0.9 GeV  (merged-only territory)
    return [("M0p%d" % i, i / 10.0) for i in range(1, 10)]


def load_signal(eos, tag, selection):
    fs = glob.glob(f"{eos}/parquet_friend/mA_{tag}/*/merged_nominal.parquet")
    if not fs:
        return None, None, 0
    cols = ["H_mass", "weight_central"]
    if selection:
        cols.append(selection)
    d = pq.read_table(fs[0], columns=cols).to_pandas()
    n_tot = len(d)
    if selection and selection in d.columns:
        d = d[d[selection] == 1]
    return d["H_mass"].to_numpy(), d["weight_central"].to_numpy(), n_tot


def cms_label(lumi):
    keep = []
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(42)
    l.SetTextSize(0.045)
    l.DrawLatex(0.13, 0.92, "#bf{CMS} #it{Preliminary}")
    l.SetTextSize(0.040)
    l.DrawLatex(0.62, 0.92, f"{lumi:.1f} fb^{{-1}} (13.6 TeV)")
    keep.append(l)
    return keep


def fit_one(tag, ma, h, w, args, ws_dir, plot_dir):
    lo, hi = args.fit_lo, args.fit_hi
    m = ROOT.RooRealVar("H_mass", "m_{ll#gamma} [GeV]", lo, hi)
    wv = ROOT.RooRealVar("w", "w", -1e9, 1e9)
    ds = ROOT.RooDataSet("sig_" + tag, "sig", ROOT.RooArgSet(m, wv), "w")
    sel = (h > lo) & (h < hi) & np.isfinite(h) & np.isfinite(w)
    for x, ww in zip(h[sel], w[sel]):
        m.setVal(float(x))
        ds.add(ROOT.RooArgSet(m, wv), float(ww))
    yield_ref = float(w[sel].sum())

    mean = ROOT.RooRealVar("mean_" + tag, "mean", 125.0, 118.0, 132.0)
    # m_llgamma resolution floor ~1 GeV; below that the DCB rails to an
    # unphysical spike and destabilises the downstream Asimov limit.
    sigma = ROOT.RooRealVar("sigma_" + tag, "sigma", 2.5, 1.0, 10.0)
    aL = ROOT.RooRealVar("aL_" + tag, "aL", 1.3, 0.1, 6.0)
    nL = ROOT.RooRealVar("nL_" + tag, "nL", 3.0, 0.1, 50.0)
    aR = ROOT.RooRealVar("aR_" + tag, "aR", 1.5, 0.1, 6.0)
    nR = ROOT.RooRealVar("nR_" + tag, "nR", 5.0, 0.1, 50.0)
    pdf = None
    try:
        # double-sided single-sigma Crystal Ball
        pdf = ROOT.RooCrystalBall("sig_model_" + tag, "dcb", m, mean, sigma,
                                  aL, nL, aR, nR)
    except Exception as exc:
        print(f"  [{tag}] RooCrystalBall unavailable ({exc}); fallback Gaussian")
        pdf = ROOT.RooGaussian("sig_model_" + tag, "g", m, mean, sigma)

    pdf.fitTo(ds, ROOT.RooFit.Save(True), ROOT.RooFit.SumW2Error(True),
              ROOT.RooFit.PrintLevel(-1))

    sig_eff = sigma.getVal()

    # ---- plot ----
    c = ROOT.TCanvas("c_" + tag, "", 800, 700)
    fr = m.frame(ROOT.RooFit.Bins(args.bins))
    ds.plotOn(fr)
    pdf.plotOn(fr, ROOT.RooFit.LineColor(ROOT.kRed + 1))
    fr.SetTitle("")
    fr.GetYaxis().SetTitle("Weighted events")
    fr.Draw()
    keep = cms_label(args.lumi)
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(42)
    t.SetTextSize(0.034)
    t.DrawLatex(0.60, 0.82, f"m_{{a}} = {ma:g} GeV  (merged)")
    t.DrawLatex(0.60, 0.77, f"#mu = {mean.getVal():.2f} GeV")
    t.DrawLatex(0.60, 0.72, f"#sigma_{{eff}} = {sig_eff:.2f} GeV")
    t.DrawLatex(0.60, 0.67, f"S(0.1pb) = {yield_ref:.1f}")
    for ext in ("png", "pdf"):
        c.SaveAs(f"{plot_dir}/sig_Hmass_{tag}.{ext}")

    # ---- workspace for combine ----
    ws = ROOT.RooWorkspace("wsig_" + tag)
    getattr(ws, "import")(pdf)
    norm = ROOT.RooRealVar("sig_model_" + tag + "_norm", "signal yield", yield_ref)
    norm.setConstant(True)
    getattr(ws, "import")(norm)
    ws.writeToFile(f"{ws_dir}/sig_ws_{tag}.root")

    return dict(tag=tag, ma=ma, mean=mean.getVal(), sigma_eff=sig_eff,
                yield_ref=yield_ref)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--eos", default=EOS_DEFAULT)
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "output"))
    ap.add_argument("--plotdir", default=os.path.join(os.path.dirname(__file__),
                                                      "../Plot/plots/merged_lowma"))
    ap.add_argument("--selection", default="pass_allcuts_merged_AN2020")
    ap.add_argument("--fit-lo", type=float, default=105.0)
    ap.add_argument("--fit-hi", type=float, default=170.0)
    ap.add_argument("--bins", type=int, default=65)
    ap.add_argument("--lumi", type=float, default=109.82)
    args = ap.parse_args()

    ws_dir = os.path.join(args.outdir, "workspaces")
    os.makedirs(ws_dir, exist_ok=True)
    os.makedirs(args.plotdir, exist_ok=True)

    results = []
    for tag, ma in mass_points():
        h, w, n_tot = load_signal(args.eos, tag, args.selection)
        if h is None or len(h) == 0:
            print(f"[{tag}] no signal parquet / empty after selection -> skip")
            continue
        n_sel = len(h)
        print(f"[{tag}] m_a={ma:g}  n_total={n_tot}  n_sel={n_sel}  "
              f"({100.0*n_sel/max(n_tot,1):.1f}% pass {args.selection})")
        res = fit_one(tag, ma, h, w, args, ws_dir, args.plotdir)
        res["n_sel"] = n_sel
        res["n_total"] = n_tot
        results.append(res)
        print(f"      mu={res['mean']:.2f}  sigma_eff={res['sigma_eff']:.2f}  "
              f"S(0.1pb)={res['yield_ref']:.2f}")

    out_json = os.path.join(args.outdir, "signal_yields.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n[signal] wrote {out_json}  ({len(results)} mass points)")
    print(f"[signal] workspaces in {ws_dir}")


if __name__ == "__main__":
    main()
