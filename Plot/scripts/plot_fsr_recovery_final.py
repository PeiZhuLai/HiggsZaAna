#!/usr/bin/env python3
"""FSR recovery on the final samples: the same events reconstructed without and with
the FSR photon added to the muon.

Answers the L3 request on Sec. 4.4 ("before/after plots demonstrating how adding the
FSR photon improves the mass resolution"). The final production stores both
reconstructions in every event (Z_mass / Z_noFSR_mass, H_mass / H_noFSR_mass), so
"before" is not a different production: it is the same selected event with the
photon left out. That isolates the effect of the photon from any change in the
selection.

Figures (all muon channel unless stated; Run 3 = 2022preEE..2024 combined):
  fsr_sig_<var>_mA<ma>[_dressed]   signal m_ll / m_llgg, all mumu events or the
                                   events where a photon was added
  fsr_sig_seff_vs_ma               relative change of sigma_eff vs m_a, muon channel
                                   per era, electron channel as the null test
  fsr_sig_dressfrac_vs_ma          fraction of mumu signal events with an FSR photon
  fsr_datamc_mll_{noFSR,FSR}       data vs Z+fake photon simulation, m_ll
  fsr_datamc_fsrpt                 pT of the added FSR photon, data vs simulation

Input : /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix/Sig_MC
        /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1/{Data,Bkg_MC}
        /eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024
          (2024 DY+jets re-produced with the overlap veto)
Output: Plot/plots/fsrRecovery_final/*.pdf|.png, fsr_final_summary.txt

Env   : LCG_104 (ROOT 6.28). ROOT 6.34 drops histograms from PDF output.
    set +u; source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh
    python3 plot_fsr_recovery_final.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import ROOT

if ROOT.gROOT.GetVersion().startswith("6.34"):
    sys.exit("ROOT 6.34 drops histograms from PDF output; use LCG_104")

# House style helpers (shared HZgamma / HZa module) and the authoritative lumi table.
HZG = "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3"
sys.path.insert(0, f"{HZG}/plot/scripts")
from hzg_style import (  # noqa: E402
    HEADROOM, LINE_WIDTH, MARKERS, PAIR_COLORS, PETROFF, below_frame,
    canvas_ratio, canvas_single, cms_header, color, house_legend, panel_label,
    save, style_axes, SIZE_LEGEND_SMALL, _keep,
)

# constants.py cannot be imported (set literals as dict keys); lumi_constants parses
# its LUMI assignment instead.
from lumi_constants import load_lumi  # noqa: E402
LUMI = load_lumi(path=f"{HZG}/higgs_dna/metaconditions/corrections/constants.py")
if not all(e in LUMI for e in ("2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024")):
    sys.exit("could not read LUMI from constants.py")

P = "/eos/project/h/htozg-dy-privatemc/pelai/HZa"
SIG = f"{P}/parquet_DNA_tmp_fsrfix/Sig_MC"
DATA = f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Data"
BKG = f"{P}/parquet_DNA_tmp_fsrfix_fpo1/Bkg_MC"
BKG24 = f"{P}/parquet_DNA_tmp_fsrfix_fpo1_dyveto/Bkg_MC_dyveto2024"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
MASSES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30]
SHOW_MASSES = [1, 5, 15, 30]
COLS = ["z_mumu", "z_ee", "Z_mass", "Z_noFSR_mass", "H_mass", "H_noFSR_mass",
        "gamma_fsr_pt", "weight_central"]
RUN3_LUMI = sum(LUMI[e] for e in ERAS)
SIM_EXTRA = "Simulation Preliminary"
OUT = Path("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/fsrRecovery_final")

VARS = {
    "mll": dict(before="Z_noFSR_mass", after="Z_mass", nbins=60, lo=60.0, hi=120.0,
                title="m_{#mu#mu} [GeV]"),
    "mllgg": dict(before="H_noFSR_mass", after="H_mass", nbins=60, lo=100.0, hi=150.0,
                  title="m_{#mu#mu#gamma#gamma} [GeV]"),
}


def read(paths):
    parts = []
    for p in ([paths] if isinstance(paths, str) else paths):
        t = pq.read_table(p, columns=COLS)
        parts.append({c: t[c].combine_chunks().to_numpy(zero_copy_only=False).astype(float)
                      for c in COLS})
    return {c: np.concatenate([d[c] for d in parts]) for c in COLS}


def concat(dicts):
    return {c: np.concatenate([d[c] for d in dicts]) for c in COLS}


def bkg_files(era):
    if era == "2024":
        dy = [f"{BKG24}/DYJetsTo2{f}_2024/merged_nominal.parquet" for f in ("E", "Mu", "Tau")]
    else:
        dy = [f"{BKG}/DYJetsToLL_{era}/merged_nominal.parquet"]
    return dy + [f"{BKG}/DYGto2LG_10to100_{era}/merged_nominal.parquet"]


def dressed_mask(d):
    return (d["z_mumu"] == 1) & (np.abs(d["Z_mass"] - d["Z_noFSR_mass"]) > 0.01)


def sigma_eff(v):
    v = np.sort(np.asarray(v, float))
    n = len(v)
    if n < 20:
        return float("nan")
    k = int(round(0.683 * n))
    return float(np.min(v[k:] - v[:n - k]) / 2.0)


def th1(name, values, weights, nbins, lo, hi):
    h = ROOT.TH1D(name, "", nbins, lo, hi)
    h.Sumw2()
    v = np.ascontiguousarray(values, dtype=np.float64)
    w = np.ascontiguousarray(weights, dtype=np.float64)
    if len(v):
        h.FillN(len(v), v, w)
    _keep.append(h)
    return h


def unit(h):
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    return h


# ───────────────────────────────────────────────────────────── signal shapes
def fig_signal_shape(d, ma, key, dressed_only, outdir, log):
    cfg = VARS[key]
    sel = dressed_mask(d) if dressed_only else (d["z_mumu"] == 1)
    b, a, w = d[cfg["before"]][sel], d[cfg["after"]][sel], d["weight_central"][sel]
    s0, s1 = sigma_eff(b), sigma_eff(a)
    tag = "_dressed" if dressed_only else ""
    h0 = unit(th1(f"h0_{key}{ma}{tag}", b, w, cfg["nbins"], cfg["lo"], cfg["hi"]))
    h1 = unit(th1(f"h1_{key}{ma}{tag}", a, w, cfg["nbins"], cfg["lo"], cfg["hi"]))
    c, up, dn = canvas_ratio(f"c_{key}{ma}{tag}")
    up.cd()
    for h, col in zip((h0, h1), PAIR_COLORS):
        h.SetLineColor(col)
        h.SetLineWidth(LINE_WIDTH)
        h.SetMarkerStyle(0)
    binw = (cfg["hi"] - cfg["lo"]) / cfg["nbins"]
    style_axes(h0, up, "", "A.U. / %.2g GeV" % binw)
    h0.GetXaxis().SetLabelSize(0.0)
    top = max(h0.GetMaximum(), h1.GetMaximum())
    h0.SetMinimum(-0.02 * top)
    h0.SetMaximum(HEADROOM * top)
    h0.Draw("hist")
    h1.Draw("hist same")
    labels = ["without FSR  #sigma_{eff} = %.2f GeV" % s0,
              "with FSR  #sigma_{eff} = %.2f GeV" % s1]
    leg = house_legend(up, labels, below_label=True, size=SIZE_LEGEND_SMALL)
    for h, t in zip((h0, h1), labels):
        leg.AddEntry(h, t, "l")
    leg.Draw()
    cms_header(up, RUN3_LUMI, "13.6", extra=SIM_EXTRA)
    panel_label(up, "Signal  m_{a} = %d GeV  #mu#mu%s" % (ma, ",  FSR-dressed events" if dressed_only else ""))
    dn.cd()
    r = h1.Clone(f"r_{key}{ma}{tag}")
    r.Divide(h0)
    r.SetLineColor(PAIR_COLORS[1])
    rlo, rhi = (0.0, 3.0) if dressed_only else (0.8, 1.2)
    r.SetMinimum(rlo)
    r.SetMaximum(rhi)
    style_axes(r, dn, cfg["title"], "with / without", ratio=True)
    r.GetXaxis().SetNdivisions(505)
    r.Draw("hist")
    one = ROOT.TLine(cfg["lo"], 1.0, cfg["hi"], 1.0)
    one.SetLineStyle(2)
    one.SetLineColor(ROOT.kGray + 2)
    one.Draw()
    _keep.extend([r, one])
    save(c, str(outdir / f"fsr_sig_{key}_mA{ma}{tag}"))
    log.append("%-6s ma=%-2d %-8s N=%7d  seff %.3f -> %.3f  (%+.2f%%)"
               % (key, ma, "dressed" if dressed_only else "all", int(sel.sum()), s0, s1,
                  100.0 * (s1 - s0) / s0))


# ───────────────────────────────────────────────────────── signal summaries
def fig_seff_vs_ma(table, outdir, log):
    """table[(ma, era)] = dict(mu=%, ee=%, mll=%, frac=...)"""
    c = canvas_single("c_seff_vs_ma")
    frame = ROOT.TH1D("frame_seff", "", 1, 0.0, 32.0)
    _keep.append(frame)
    vals = [v["mu"] for v in table.values()] + [v["mll"] for v in table.values()]
    lo = min(-6.0, 1.15 * min(vals))
    frame.SetMinimum(lo)
    frame.SetMaximum(-1.35 * lo)
    style_axes(frame, c, "m_{a} [GeV]", "#Delta#sigma_{eff} / #sigma_{eff} [%]")
    frame.Draw("axis")
    zero = ROOT.TLine(0.0, 0.0, 32.0, 0.0)
    zero.SetLineStyle(2)
    zero.SetLineColor(ROOT.kGray + 2)
    zero.Draw()
    _keep.append(zero)
    entries, graphs = [], []
    for i, era in enumerate(ERAS):
        g = ROOT.TGraph(len(MASSES), np.array([m - 0.6 + 0.3 * i for m in MASSES], "d"),
                        np.array([table[(m, era)]["mu"] for m in MASSES], "d"))
        g.SetMarkerStyle(MARKERS[i])
        g.SetMarkerSize(1.4)
        g.SetMarkerColor(color(PETROFF[i]))
        g.SetLineColor(color(PETROFF[i]))
        g.Draw("P same")
        graphs.append(g)
        entries.append(era)
    ge = ROOT.TGraph(len(MASSES), np.array(MASSES, "d"),
                     np.array([np.mean([table[(m, e)]["ee"] for e in ERAS]) for m in MASSES], "d"))
    ge.SetMarkerStyle(24)
    ge.SetMarkerSize(1.4)
    ge.SetMarkerColor(ROOT.kBlack)
    ge.Draw("P same")
    graphs.append(ge)
    entries.append("ee (no recovery)")
    _keep.extend(graphs)
    leg = house_legend(c, entries, ncols=2, size=0.036, below_label=True)
    for g, t in zip(graphs, entries):
        leg.AddEntry(g, t, "p")
    leg.Draw()
    cms_header(c, RUN3_LUMI, "13.6", extra=SIM_EXTRA)
    panel_label(c, "Signal  #sigma_{eff}(m_{#font[12]{ll}#gamma#gamma}), with vs without FSR")
    save(c, str(outdir / "fsr_sig_seff_vs_ma"))

    c2 = canvas_single("c_frac_vs_ma")
    frame2 = ROOT.TH1D("frame_frac", "", 1, 0.0, 32.0)
    _keep.append(frame2)
    fmax = max(v["frac"] for v in table.values())
    frame2.SetMinimum(0.0)
    frame2.SetMaximum(2.4 * fmax)
    style_axes(frame2, c2, "m_{a} [GeV]", "FSR-dressed fraction [%]")
    frame2.Draw("axis")
    graphs2 = []
    for i, era in enumerate(ERAS):
        g = ROOT.TGraph(len(MASSES), np.array(MASSES, "d"),
                        np.array([table[(m, era)]["frac"] for m in MASSES], "d"))
        g.SetMarkerStyle(MARKERS[i])
        g.SetMarkerSize(1.4)
        g.SetMarkerColor(color(PETROFF[i]))
        g.SetLineColor(color(PETROFF[i]))
        g.SetLineWidth(2)
        g.Draw("PL same")
        graphs2.append(g)
    _keep.extend(graphs2)
    leg2 = house_legend(c2, ERAS, ncols=2, size=0.036, below_label=True)
    for g, t in zip(graphs2, ERAS):
        leg2.AddEntry(g, t, "p")
    leg2.Draw()
    cms_header(c2, RUN3_LUMI, "13.6", extra=SIM_EXTRA)
    panel_label(c2, "Signal  #mu#mu")
    save(c2, str(outdir / "fsr_sig_dressfrac_vs_ma"))


# ─────────────────────────────────────────────────────────────── data vs MC
def fig_datamc(data, mc, key_before, key_after, name, xtitle, nbins, lo, hi, outdir,
               log, sel_fn, note, logy=False):
    sd, sm = sel_fn(data), sel_fn(mc)
    for tagname, col in (("noFSR", key_before), ("FSR", key_after)):
        if col is None:
            continue
        hd = th1(f"hd_{name}_{tagname}", data[col][sd], np.ones(int(sd.sum())), nbins, lo, hi)
        hm = th1(f"hm_{name}_{tagname}", mc[col][sm], mc["weight_central"][sm], nbins, lo, hi)
        if hm.Integral() > 0:
            hm.Scale(hd.Integral() / hm.Integral())
        c, up, dn = canvas_ratio(f"c_{name}_{tagname}")
        up.cd()
        if logy:
            up.SetLogy()
        hm.SetLineColor(ROOT.kAzure + 2)
        hm.SetFillColor(color(PETROFF[0]))
        hm.SetLineWidth(LINE_WIDTH)
        hd.SetMarkerStyle(20)
        hd.SetMarkerSize(1.2)
        hd.SetLineColor(ROOT.kBlack)
        top = max(hd.GetMaximum(), hm.GetMaximum())
        if logy:
            hm.SetMaximum(2000.0 * top)
            hm.SetMinimum(0.5)
        else:
            hm.SetMaximum(HEADROOM * top)
            hm.SetMinimum(-0.02 * top)
        binw = (hi - lo) / nbins
        style_axes(hm, up, "", "Events / %.2g GeV" % binw)
        hm.GetXaxis().SetLabelSize(0.0)
        hm.Draw("hist")
        hd.Draw("E1 same")
        sdat, smc = sigma_eff(data[col][sd]), sigma_eff(mc[col][sm])
        lbl = "without FSR" if tagname == "noFSR" else "with FSR"
        if key_before is None:
            lbl = "added FSR photon"
        texts = ["Data", "Z+fake #gamma sim."]
        leg = house_legend(up, texts, below_label=True)
        leg.AddEntry(hd, texts[0], "pe")
        leg.AddEntry(hm, texts[1], "f")
        leg.Draw()
        cms_header(up, RUN3_LUMI, "13.6", extra="Preliminary")
        panel_label(up, "#mu#mu  %s" % lbl, [note])
        dn.cd()
        r = hd.Clone(f"r_{name}_{tagname}")
        r.Divide(hm)
        r.SetMinimum(0.5)
        r.SetMaximum(1.5)
        style_axes(r, dn, xtitle, "Data / sim.", ratio=True)
        r.Draw("E1")
        one = ROOT.TLine(lo, 1.0, hi, 1.0)
        one.SetLineStyle(2)
        one.SetLineColor(ROOT.kGray + 2)
        one.Draw()
        _keep.extend([r, one])
        save(c, str(outdir / f"fsr_datamc_{name}_{tagname}" if key_before else outdir / f"fsr_datamc_{name}"))
        log.append("data/MC %-8s %-6s Ndata=%8d  sigma_eff data %.3f  sim %.3f"
                   % (name, tagname, int(sd.sum()), sdat, smc))


def fig_datamc_shift(data, mc, outdir, log, lo=70.0, hi=110.0, nbins=40):
    """What the recovery does to m_mumu, in data and in simulation side by side:
    each shape without and with the photon, and the with/without ratio of both."""
    sd, sm = data["z_mumu"] == 1, mc["z_mumu"] == 1
    hs = {}
    for src, d, s, w in (("data", data, sd, np.ones(int(sd.sum()))),
                         ("sim", mc, sm, mc["weight_central"][sm])):
        for tag, col in (("noFSR", "Z_noFSR_mass"), ("FSR", "Z_mass")):
            hs[(src, tag)] = unit(th1(f"hs_{src}_{tag}", d[col][s], w, nbins, lo, hi))
    c, up, dn = canvas_ratio("c_mll_shift")
    up.cd()
    style = {("sim", "noFSR"): (ROOT.kBlack, 2), ("sim", "FSR"): (ROOT.kRed + 1, 1)}
    first = hs[("sim", "noFSR")]
    top = max(h.GetMaximum() for h in hs.values())
    first.SetMaximum(2.4 * top)
    first.SetMinimum(-0.02 * top)
    style_axes(first, up, "", "A.U. / %.2g GeV" % ((hi - lo) / nbins))
    first.GetXaxis().SetLabelSize(0.0)
    for key, (col, ls) in style.items():
        h = hs[key]
        h.SetLineColor(col)
        h.SetLineStyle(ls)
        h.SetLineWidth(LINE_WIDTH)
        h.SetMarkerStyle(0)
        h.Draw("hist" if key == ("sim", "noFSR") else "hist same")
    for key, (col, mk) in ((("data", "noFSR"), (ROOT.kBlack, 24)), (("data", "FSR"), (ROOT.kRed + 1, 20))):
        h = hs[key]
        h.SetMarkerStyle(mk)
        h.SetMarkerSize(1.2)
        h.SetMarkerColor(col)
        h.SetLineColor(col)
        h.Draw("E1 same")
    texts = ["Data, without FSR", "Data, with FSR", "Sim., without FSR", "Sim., with FSR"]
    leg = house_legend(up, texts, ncols=1, size=0.038, below_label=True)
    leg.AddEntry(hs[("data", "noFSR")], texts[0], "pe")
    leg.AddEntry(hs[("data", "FSR")], texts[1], "pe")
    leg.AddEntry(hs[("sim", "noFSR")], texts[2], "l")
    leg.AddEntry(hs[("sim", "FSR")], texts[3], "l")
    leg.Draw()
    cms_header(up, RUN3_LUMI, "13.6", extra="Preliminary")
    panel_label(up, "#mu#mu  preselection")
    dn.cd()
    rd = hs[("data", "FSR")].Clone("r_shift_data")
    rd.Divide(hs[("data", "noFSR")])
    rm = hs[("sim", "FSR")].Clone("r_shift_sim")
    rm.Divide(hs[("sim", "noFSR")])
    rm.SetLineStyle(1)
    rm.SetMinimum(0.9)
    rm.SetMaximum(1.1)
    style_axes(rm, dn, "m_{#mu#mu} [GeV]", "with / without", ratio=True)
    rm.Draw("hist")
    rd.Draw("E1 same")
    one = ROOT.TLine(lo, 1.0, hi, 1.0)
    one.SetLineStyle(2)
    one.SetLineColor(ROOT.kGray + 2)
    one.Draw()
    _keep.extend([rd, rm, one])
    save(c, str(outdir / "fsr_datamc_mll_shift"))
    for src in ("data", "sim"):
        a, b = hs[(src, "FSR")], hs[(src, "noFSR")]
        i0, i1 = a.FindBin(86.0), a.FindBin(96.0 - 1e-6)
        log.append("%-4s fraction of mumu events with 86<m_mumu<96: without %.4f  with %.4f  (%+.2f%%)"
                   % (src, b.Integral(i0, i1), a.Integral(i0, i1),
                      100.0 * (a.Integral(i0, i1) / b.Integral(i0, i1) - 1.0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--skip-datamc", action="store_true")
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    log = ["Run 3 lumi (constants.py, %s) = %.2f /fb" % ("+".join(ERAS), RUN3_LUMI)]

    table, sig = {}, {}
    for ma in MASSES:
        per_era = []
        for era in ERAS:
            d = read(f"{SIG}/mA_M{ma}_{era}/merged_nominal.parquet")
            per_era.append(d)
            mu, ee = d["z_mumu"] == 1, d["z_ee"] == 1
            rel = lambda b, a, s: 100.0 * (sigma_eff(d[a][s]) - sigma_eff(d[b][s])) / sigma_eff(d[b][s])
            table[(ma, era)] = dict(mu=rel("H_noFSR_mass", "H_mass", mu),
                                    ee=rel("H_noFSR_mass", "H_mass", ee),
                                    mll=rel("Z_noFSR_mass", "Z_mass", mu),
                                    frac=100.0 * dressed_mask(d).sum() / mu.sum())
        sig[ma] = concat(per_era)
    for ma in SHOW_MASSES:
        for k in VARS:
            for dressed in (False, True):
                fig_signal_shape(sig[ma], ma, k, dressed, outdir, log)
    fig_seff_vs_ma(table, outdir, log)
    mu_all = np.array([v["mu"] for v in table.values()])
    ee_all = np.array([v["ee"] for v in table.values()])
    fr_all = np.array([v["frac"] for v in table.values()])
    log.append("per (ma, era), 70 points: mllgg mu mean %+.2f%% median %+.2f%% range %+.2f..%+.2f%%, improved %d/70;"
               " ee max |change| %.4f%%; dressed fraction %.2f-%.2f%%"
               % (mu_all.mean(), np.median(mu_all), mu_all.min(), mu_all.max(), (mu_all < 0).sum(),
                  np.abs(ee_all).max(), fr_all.min(), fr_all.max()))

    if not args.skip_datamc:
        data = concat([read(f"{DATA}/Data_{e}/merged_nominal.parquet") for e in ERAS])
        mc = concat([read(bkg_files(e)) for e in ERAS])
        mumu = lambda d: d["z_mumu"] == 1
        fig_datamc(data, mc, "Z_noFSR_mass", "Z_mass", "mll", "m_{#mu#mu} [GeV]",
                   70, 50.0, 120.0, outdir, log, mumu, "preselection, sim. normalized to data",
                   logy=True)
        fig_datamc(data, mc, None, "gamma_fsr_pt", "fsrpt", "p_{T}^{FSR #gamma} [GeV]",
                   30, 0.0, 30.0, outdir, log, dressed_mask, "FSR-dressed events, sim. norm. to data")
        fig_datamc_shift(data, mc, outdir, log)
        fd = 100.0 * dressed_mask(data).sum() / mumu(data).sum()
        w = mc["weight_central"]
        fm = 100.0 * w[dressed_mask(mc)].sum() / w[mumu(mc)].sum()
        log.append("dressed fraction mumu preselection: data %.2f%%  sim %.2f%%" % (fd, fm))
        log.append("FSR photon pT median: data %.2f  sim %.2f GeV"
                   % (np.median(data["gamma_fsr_pt"][dressed_mask(data)]),
                      np.median(mc["gamma_fsr_pt"][dressed_mask(mc)])))
    (outdir / "fsr_final_summary.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
