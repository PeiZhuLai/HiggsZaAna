#!/usr/bin/env python3
"""Combined-trigger-SF validation in the zee_zmmg control sample's clean dilepton CRs.

通用化 muon/electron 兩 channel（同一份程式，--channel muon|electron|both）。
在乾淨 CR 內（muon: pass_dimu_cr = 2mu,0gamma,80<m_mumu<100；electron: pass_diel_cr =
2e,0mu,0gamma,80<m_ee<100）算三條 ratio = N(pass OR trig) / N(pass di-lepton trig)：
  Data       : data(OR) / data(di-lep)
  MC nominal : mc(OR, lumi only) / mc(di-lep, di-lep-trig weight)
  MC combined: mc(OR, combined-trig weight) / mc(di-lep, di-lep-trig weight)
畫 ratio vs lead/sublead lepton pT，下方 Data/MC panel，含統計誤差。
目標：檢查 combined trigger SF 是否讓 MC 更貼 data（尤其 low-pT）。

注意：框架目前只有「combined(OR) HLT SF」= weight_hlt_sf_central，沒有專屬的
di-lepton-only trigger SF；與既有 plot_SF_validation.py 一樣，di-lep-trig 權重退回用
weight_hlt_sf_central（可用 --dilep-weight 改）。lumi 權重用 weight_central_initial(genweight)。
Env: higgs-alp-ana（ROOT 6.24 + pyarrow）。
輸出: Plot/plots/trigger_SF_validation/<era>/<channel>/trigger_ratio_{lead,sublead}_<lep>Pt.{pdf,png}
"""
import argparse, os, glob
from array import array
import numpy as np
import pyarrow.parquet as pq
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLegendBorderSize(0); ROOT.gStyle.SetErrorX(0.5); ROOT.gStyle.SetEndErrorSize(0)

DEFAULT_ROOT = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_zee_zmmg_ctrl"
LUMI_MAP = {"2022preEE":7.9804,"2022postEE":26.6717,"2023preBPix":18.063,"2023postBPix":9.693,"2024":108.83}
DYLL_ERAS = ["2022preEE","2022postEE","2023preBPix","2023postBPix"]

PT_MU_LEAD = [20,25,30,35,40,50,60,80,120]
# sublead 起點拉高：di-lepton 低 leg 門檻(μ~8, e~12)以下分母幾乎空 → ratio 爆掉，故 drop
PT_MU_SUB  = [10,12,15,20,25,30,40,60,120]
PT_E_LEAD  = [25,30,35,40,50,60,80,120]
PT_E_SUB   = [15,20,25,30,40,50,70,120]

# 每個 channel 的欄位/設定，muon 與 electron 完全共用後面流程
CHANNELS = {
    "muon": dict(
        data_sub="data_muon/Data_tnp_zmmg_%s", cr="pass_dimu_cr",
        num="pass_simu_or_dimu", den="pass_dimu",
        lead_pt="dimu_cr_lead_muon_pt", sub_pt="dimu_cr_sublead_muon_pt",
        lead_title="Lead muon p_{T} [GeV]", sub_title="Sublead muon p_{T} [GeV]",
        lead_bins=PT_MU_LEAD, sub_bins=PT_MU_SUB,
        ratio_ytitle="OR(single+di-#mu) / di-#mu trigger",
    ),
    "electron": dict(
        data_sub="data_egamma/Data_%s", cr="pass_diel_cr",
        num="pass_siel_or_diel", den="pass_diel",
        lead_pt="diel_cr_lead_electron_pt", sub_pt="diel_cr_sublead_electron_pt",
        lead_title="Lead electron p_{T} [GeV]", sub_title="Sublead electron p_{T} [GeV]",
        lead_bins=PT_E_LEAD, sub_bins=PT_E_SUB,
        ratio_ytitle="OR(single+di-e) / di-e trigger",
    ),
}


def mc_paths(base, eras):
    return [p for e in eras for p in glob.glob(f"{base}/mc/DYJetsToLL_{e}/merged_nominal.parquet")]

def data_paths(base, eras, sub):
    return [p for e in eras for p in glob.glob(f"{base}/{sub % e}/merged_nominal.parquet")]


def load(paths, cols):
    out = {c: [] for c in cols}
    for p in paths:
        avail = set(pq.ParquetFile(p).schema.names)
        rd = [c for c in cols if c in avail]
        t = pq.read_table(p, columns=rd)
        n = pq.ParquetFile(p).metadata.num_rows
        for c in cols:
            out[c].append(np.asarray(t[c]) if c in rd else np.zeros(n))
    return {c: (np.concatenate(v) if v else np.array([])) for c, v in out.items()}


def fill(name, vals, w, edges):
    edges_a = array("d", [float(x) for x in edges])
    h = ROOT.TH1D(name, name, len(edges) - 1, edges_a); h.Sumw2()
    vals = np.asarray(vals, float); w = np.asarray(w, float)
    g = np.isfinite(vals) & np.isfinite(w) & (vals > -900)
    vals, w = vals[g], w[g]
    be = np.asarray(edges, float)
    cnt, _ = np.histogram(vals, bins=be, weights=w)
    er2, _ = np.histogram(vals, bins=be, weights=w * w)
    for i in range(len(edges) - 1):
        h.SetBinContent(i + 1, cnt[i]); h.SetBinError(i + 1, er2[i] ** 0.5)
    h.SetDirectory(0)
    return h


def ratio(num, den, name):
    if den.Integral() <= 0:
        return None
    r = num.Clone(name); r.SetDirectory(0); r.Divide(den)
    return r


def cms_labels(pad, lumi):
    pad.cd(); t = ROOT.TLatex(); t.SetNDC(True); t.SetTextFont(42)
    t.SetTextAlign(13); t.SetTextSize(0.060); t.DrawLatex(0.17, 0.96, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.SetTextSize(0.055); t.DrawLatex(0.95, 0.915, f"{lumi:.2f} fb^{{-1}} (13.6 TeV)")


def draw(data_r, nom_r, comb_r, xtitle, ytitle, lumi, out, logx):
    blue = ROOT.TColor.GetColor("#1F78B4"); orange = ROOT.TColor.GetColor("#E66101")
    for h in (data_r, nom_r, comb_r):
        h.SetTitle(""); h.SetStats(False); h.SetFillStyle(0)
    data_r.SetMarkerStyle(20); data_r.SetMarkerSize(1.3); data_r.SetLineColor(ROOT.kBlack); data_r.SetLineWidth(3); data_r.SetMarkerColor(ROOT.kBlack)
    nom_r.SetMarkerStyle(24); nom_r.SetMarkerSize(1.0); nom_r.SetLineColor(blue); nom_r.SetLineWidth(3); nom_r.SetMarkerColor(blue)
    comb_r.SetMarkerStyle(25); comb_r.SetMarkerSize(1.0); comb_r.SetLineColor(orange); comb_r.SetLineWidth(3); comb_r.SetMarkerColor(orange)

    c = ROOT.TCanvas("c", "", 800, 800)
    up = ROOT.TPad("up", "", 0, 0.30, 1, 1); lo = ROOT.TPad("lo", "", 0, 0, 1, 0.30)
    up.SetBottomMargin(0.03); up.SetLeftMargin(0.15); up.SetRightMargin(0.05)
    lo.SetTopMargin(0.04); lo.SetBottomMargin(0.45); lo.SetLeftMargin(0.15); lo.SetRightMargin(0.05)
    if logx:
        up.SetLogx(); lo.SetLogx()
    up.Draw(); lo.Draw()

    up.cd()
    ymax = max(data_r.GetMaximum(), nom_r.GetMaximum(), comb_r.GetMaximum())
    nom_r.SetMinimum(0.75); nom_r.SetMaximum(min(max(1.4, ymax * 1.2), 2.5))  # cap 避免單一大 bin 壓扁
    nom_r.GetXaxis().SetLabelSize(0)
    nom_r.GetYaxis().SetTitle(ytitle); nom_r.GetYaxis().SetTitleSize(0.055); nom_r.GetYaxis().SetLabelSize(0.05); nom_r.GetYaxis().SetTitleOffset(1.25)
    nom_r.Draw("E1"); comb_r.Draw("E1 same"); data_r.Draw("E1 same")
    line = ROOT.TLine(nom_r.GetBinLowEdge(1), 1.0, nom_r.GetBinLowEdge(nom_r.GetNbinsX()+1), 1.0)
    line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(2); line.Draw("same")
    leg = ROOT.TLegend(0.20, 0.68, 0.60, 0.90); leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.04)
    leg.AddEntry(data_r, "Data", "lep"); leg.AddEntry(nom_r, "MC (Wo Trigger SF)", "lp"); leg.AddEntry(comb_r, "MC (Wi combined Trigger SF)", "lp"); leg.Draw()
    cms_labels(up, lumi)

    lo.cd()
    edges = array("d", [nom_r.GetBinLowEdge(i + 1) for i in range(nom_r.GetNbinsX() + 1)])
    frame = ROOT.TH1D("frame", "", nom_r.GetNbinsX(), edges); frame.SetStats(False); frame.SetDirectory(0)
    frame.GetYaxis().SetTitle("Data / MC"); frame.GetXaxis().SetTitle(xtitle)
    frame.GetYaxis().SetRangeUser(0.6, 1.4); frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.14); frame.GetYaxis().SetTitleOffset(0.45); frame.GetYaxis().SetLabelSize(0.12)
    frame.GetXaxis().SetTitleSize(0.14); frame.GetXaxis().SetTitleOffset(1.3); frame.GetXaxis().SetLabelSize(0.12)
    frame.Draw("axis")
    l2 = ROOT.TLine(edges[0], 1.0, edges[-1], 1.0); l2.SetLineColor(ROOT.kGray + 2); l2.SetLineStyle(2); l2.Draw("same")
    keep = []
    for mc_r, col, mk in ((nom_r, blue, 24), (comb_r, orange, 25)):
        dm = ratio(data_r, mc_r, "dm_" + mc_r.GetName())
        if dm:
            dm.SetMarkerStyle(mk); dm.SetMarkerSize(1.2); dm.SetMarkerColor(col); dm.SetLineColor(col); dm.SetLineWidth(3); dm.Draw("E1 same"); keep.append(dm)
    c.SaveAs(out + ".pdf"); c.SaveAs(out + ".png"); c.Close()


def run_channel(ch, base, eras, lumi, outdir, dilep_weight):
    cfg = CHANNELS[ch]
    cols = [cfg["cr"], cfg["num"], cfg["den"], cfg["lead_pt"], cfg["sub_pt"]]
    mc = load(mc_paths(base, eras), cols + ["weight_central_initial", "weight_hlt_sf_central"])
    da = load(data_paths(base, eras, cfg["data_sub"]), cols)
    if mc[cfg["cr"]].size == 0 or da[cfg["cr"]].size == 0:
        print(f"  [skip] {ch}: empty mc/data"); return
    # masks
    mc_cr = mc[cfg["cr"]].astype(bool); da_cr = da[cfg["cr"]].astype(bool)
    mc_num = mc_cr & mc[cfg["num"]].astype(bool); mc_den = mc_cr & mc[cfg["den"]].astype(bool)
    da_num = da_cr & da[cfg["num"]].astype(bool); da_den = da_cr & da[cfg["den"]].astype(bool)
    # weights
    wci = np.where(np.isfinite(mc["weight_central_initial"]), mc["weight_central_initial"], 0.0)
    whlt = np.where(np.isfinite(mc["weight_hlt_sf_central"]) & (mc["weight_hlt_sf_central"] > 0.01) & (mc["weight_hlt_sf_central"] < 100), mc["weight_hlt_sf_central"], 1.0)
    w_dilep = wci * whlt if dilep_weight == "hlt" else wci  # di-lep-trig 權重（預設退回 combined HLT SF）
    w_comb = wci * whlt                                     # combined numerator 權重
    os.makedirs(outdir, exist_ok=True)
    for which, ptcol, title, bins in (("lead", cfg["lead_pt"], cfg["lead_title"], cfg["lead_bins"]),
                                       ("sublead", cfg["sub_pt"], cfg["sub_title"], cfg["sub_bins"])):
        tag = f"{ch}_{which}"
        d_num = fill(f"dn_{tag}", da[ptcol][da_num], np.ones(da_num.sum()), bins)
        d_den = fill(f"dd_{tag}", da[ptcol][da_den], np.ones(da_den.sum()), bins)
        m_num_nom = fill(f"mnn_{tag}", mc[ptcol][mc_num], wci[mc_num], bins)
        m_num_comb = fill(f"mnc_{tag}", mc[ptcol][mc_num], w_comb[mc_num], bins)
        m_den = fill(f"md_{tag}", mc[ptcol][mc_den], w_dilep[mc_den], bins)
        data_r = ratio(d_num, d_den, f"r_data_{tag}")
        nom_r = ratio(m_num_nom, m_den, f"r_nom_{tag}")
        comb_r = ratio(m_num_comb, m_den, f"r_comb_{tag}")
        if not (data_r and nom_r and comb_r):
            print(f"  [skip] {tag}: empty denominator"); continue
        out = os.path.join(outdir, f"trigger_ratio_{which}_{ch}Pt")
        draw(data_r, nom_r, comb_r, title, cfg["ratio_ytitle"], lumi, out, logx=True)
        print(f"  [done] {tag} (data num/den={da_num.sum()}/{da_den.sum()} mc num/den={mc_num.sum()}/{mc_den.sum()})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=DEFAULT_ROOT)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "plots", "trigger_SF_validation"))
    ap.add_argument("--channel", default="both", choices=["muon", "electron", "both"])
    ap.add_argument("--eras", default="all")
    ap.add_argument("--dilep-weight", default="hlt", choices=["hlt", "lumi"], help="di-lep-trig 分母權重(預設hlt=退回combined HLT SF)")
    args = ap.parse_args()

    chans = ["muon", "electron"] if args.channel == "both" else [args.channel]
    groups = [(e, [e]) for e in DYLL_ERAS] + [("run3_all", DYLL_ERAS)] if args.eras == "all" \
        else [(e, [e]) for e in args.eras.split(",") if e.strip()]
    print(f"[cfg] base={args.base} out={args.out} channels={chans} dilep_weight={args.dilep_weight}")
    for label, eras in groups:
        if not mc_paths(args.base, eras):
            print(f"[skip] {label}: no DYJetsToLL MC"); continue
        lumi = sum(LUMI_MAP.get(e, 0.0) for e in eras)
        print(f"\n[era] {label} lumi={lumi:.2f}")
        for ch in chans:
            run_channel(ch, args.base, eras, lumi, os.path.join(args.out, label, ch), args.dilep_weight)
    print("\nAll done")


if __name__ == "__main__":
    main()
