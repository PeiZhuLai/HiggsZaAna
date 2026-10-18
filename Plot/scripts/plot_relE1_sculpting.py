#!/usr/bin/env python3
"""Background m_llgg sculpting with and without the leading-photon sigma_E/E (AN appendix on the low-mass input choice).

Two low-mass BDTs are trained with the variable-search recipe (minscan4 search_job.train, same samples, both with the
signal of m_a = 1, 2, 3 trained inside 115-133 GeV):
  relE2         the production inputs (photon isolation, dR(g1,g2), R9(g1), sigma_E/E(g2), param = (m_gg - m_a)/m_llgg)
  relE2+relE1   the same plus sigma_E/E(g1)
At each m_a the cut keeps the number of data-sideband events that the deployed model (HZaMVA/using, low_relE2_w123)
keeps at its working point in the same sample (82 / 152 / 112, half of the previous production counts 165 / 305 / 224
that minscan4 used), so both models select the same number of background events as the analysis. Drawn, one figure per m_a: the DY simulation m_llgg shape after the
cut divided by the shape before it (each normalized to unit area); a flat ratio means no sculpting. Error bars: MC
statistics of the passing events, sqrt(sum w^2) (the inclusive shape has about 100 times more events).
R / R_lo / R_hi (search_job.evaluate, 100 Poisson bootstraps) are printed. With the minscan4 counts (N_TARGET=prev)
they reproduce minscan4/summary.txt, but only in the environment of the search (LCG_106: xgboost 2.0.3, one thread);
xgboost 2.1.1 with four threads moves R_lo at m_a = 1 from 0.76 / 0.23 to 0.96 / 0.56.
Numbers: relE1_sculpting_bkg.txt. Env: LCG_106 view, OMP_NUM_THREADS=1 NTHREADS=1.
Output: Plot/plots/var_mllgg_dependence/relE1_sculpting_bkg_mA{1,2,3}.pdf
"""
import ast
import json
import os
import pickle
import sys

import numpy as np
import ROOT

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/varsearch_validate/minscan4/code")
import search_job as sj  # noqa: E402

CONSTANTS = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py"
LUMI = next(ast.literal_eval(n.value) for n in ast.parse(open(CONSTANTS).read()).body
            if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "LUMI")

SAMPLES = "/eos/home-p/pelai/HZa/bdt_varsearch_20261013/samples"
OUT = os.path.join(os.path.dirname(HERE), "plots", "var_mllgg_dependence")
BASE = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho2_relEerr", "param"]
MODELS = [("relE2", BASE, ROOT.kBlue + 1, 20, "current inputs"),
          ("relE2+relE1", BASE + ["pho1_relEerr"], ROOT.kRed + 1, 25, "current inputs + #sigma_{E}/E(#gamma_{1})")]
USING = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using"
WP = {1: 0.98115, 2: 0.98837, 3: 0.99226}   # Plot/scripts/scan_score_R_significance.py FIXED_WP
EDGES = np.array([95.0, 105.0, 115.0, 125.0, 135.0, 145.0, 160.0, 180.0])


def ratio_hist(name, m, w, passed):
    hp, _ = np.histogram(m[passed], EDGES, weights=w[passed])
    e2, _ = np.histogram(m[passed], EDGES, weights=w[passed] ** 2)
    ht, _ = np.histogram(m, EDGES, weights=w)
    r = (hp / hp.sum()) / (ht / ht.sum())
    h = ROOT.TH1D(name, "", len(EDGES) - 1, EDGES)
    for i in range(len(r)):
        h.SetBinContent(i + 1, r[i]); h.SetBinError(i + 1, r[i] * np.sqrt(e2[i]) / hp[i] if hp[i] > 0 else 0.0)
    return h


def main():
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
    os.makedirs(OUT, exist_ok=True)
    meta = json.load(open(os.path.join(SAMPLES, "meta.json")))
    names = ["data_sb", "dy", "train_bkg"] + ["%s_sig_m%d" % (t, m) for t in ("train", "test") for m in sj.LOW]
    S = {n: sj.Sample(os.path.join(SAMPLES, n + ".npz")) for n in names}
    lumi = sum(LUMI[e] for e in meta["eras"])
    if os.environ.get("N_TARGET") == "prev":
        target = {m: int(meta["prod_data_sb_count"][str(m)]) for m in sj.LOW}
    else:
        dmeta = json.load(open(os.path.join(USING, "model_Za_BDT_lowmass_run3.meta.json")))
        dep = pickle.load(open(os.path.join(USING, "model_Za_BDT_lowmass_run3.pkl"), "rb"))
        target = {m: int((dep.predict_proba(sj.matrix(S["data_sb"], dmeta["features"], float(m)))[:, 1] > WP[m]).sum())
                  for m in sj.LOW}
    print("data-sideband targets", target, flush=True)

    res, lines = {}, []
    for tag, feats, _, _, _ in MODELS:
        cfg = {"feats_full": feats, "win_mA1": True, "win_masses": [2, 3]}
        clf = sj.train(S, cfg, sj.LOW)
        for m in sj.LOW:
            r = sj.evaluate(S, cfg, clf, m, target[m], np.random.default_rng(20261013))
            dy = S["dy"]
            sb = clf.predict_proba(sj.matrix(dy, feats, float(m)))[:, 1]
            bm, bw = dy["H_m"], dy["factor"]
            ok = np.isfinite(bw) & np.isfinite(sb)
            res[(tag, m)] = ratio_hist("r_%s_%d" % (tag.replace("+", "p"), m), bm[ok], bw[ok], sb[ok] > r["cut"])
            lines.append("%-12s m_a=%d cut %.5f n_data_sb %d  R %.2f+-%.2f  R_lo %.2f+-%.2f  R_hi %.2f+-%.2f  Z %.3f" % (
                tag, m, r["cut"], r["n_data_sb"], r["R"], r["R_err"], r["R_lo"], r["R_lo_err"], r["R_hi"],
                r["R_hi_err"], r["Z"]))
            print(lines[-1], flush=True)

    for m in sj.LOW:
        c = ROOT.TCanvas("c%d" % m, "", 800, 700)
        c.SetLeftMargin(0.15); c.SetRightMargin(0.04); c.SetTopMargin(0.08); c.SetBottomMargin(0.14)
        frame = ROOT.TH1D("f%d" % m, ";m_{ll#gamma#gamma} (GeV);shape after / before the cut", 1, 95, 180)
        frame.SetMinimum(0.0); frame.SetMaximum(3.0)
        for ax in (frame.GetXaxis(), frame.GetYaxis()):
            ax.SetTitleSize(0.055); ax.SetLabelSize(0.050)
        frame.GetYaxis().SetTitleOffset(1.25)
        frame.Draw("AXIS")
        box = ROOT.TBox(sj.WIN[0], 0.0, sj.WIN[1], 3.0); box.SetFillColorAlpha(ROOT.kGray, 0.4); box.Draw()
        one = ROOT.TLine(95, 1, 180, 1); one.SetLineStyle(2); one.Draw()
        leg = ROOT.TLegend(0.30, 0.62, 0.95, 0.79); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.035)
        for k, (tag, _, col, mk, lab) in enumerate(MODELS):
            h = res[(tag, m)]
            h.SetLineColor(col); h.SetMarkerColor(col); h.SetMarkerStyle(mk); h.SetMarkerSize(1.4); h.SetLineWidth(2)
            h.Draw("E1 X0 SAME")
            leg.AddEntry(h, lab, "lp")
        leg.AddEntry(box, "signal training window", "f")
        leg.Draw()
        frame.Draw("AXIS SAME")
        t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42)
        t.SetTextSize(0.040); t.DrawLatex(0.31, 0.82, "Z/#gamma* + jets simulation, m_{a} = %d GeV" % m)
        t.SetTextSize(0.050); t.DrawLatex(0.15, 0.935, "#bf{CMS} #it{Preliminary}")
        t.SetTextSize(0.040); t.SetTextAlign(31); t.DrawLatex(0.96, 0.935, "%.2f fb^{-1} (13.6 TeV)" % lumi)
        c.SaveAs(os.path.join(OUT, "relE1_sculpting_bkg_mA%d.pdf" % m))
    open(os.path.join(OUT, "relE1_sculpting_bkg.txt"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
