#!/usr/bin/env python3
"""SF-validation plots for the zee_zmmg control sample.

讀 control-sample merged parquet（pyarrow，直接讀，不轉 root），畫 data/MC 的
「套 SF vs 不套 SF」比較（上圖分佈 + 下圖 Data/MC ratio），驗證 lepton / photon SF：
  - Z->ee   (data=EGamma) : electron reco / wplid(ID) SF
  - Z->mumugamma (data=Muon): muon reco / looseid SF、photon id / csev SF

MC 用 DYJetsToLL（與既有 plot_SF_validation.py 慣例一致；只有 2022/2023 4 個 era 有），
每張圖把 MC 縮放到 data（只比形狀，不需 lumi/xs）。
Env: higgs-alp-ana（ROOT 6.24 + pyarrow）。
輸出: Plot/plots/SF_validation_zee_zmmg/<era>/<sf>/<var>.{pdf,png}
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
# 所有會乘進 weight_central 的 SF central 分支。注意：electron SF 在 SelectedElectron
# 為空（即整個 z_mumu 區）時回傳 NaN，會毒化 weight_central。因此這裡用
# weight_central_initial（genweight）× 各 SF（NaN→1，因空集合的 SF 本應=1）自行重建權重。
SF_CENTRAL = [
    "weight_hlt_sf_central", "weight_pu_reweight_sf_central",
    "weight_electron_reco_sf_SelectedElectron_central", "weight_electron_wplid_sf_SelectedElectron_central",
    "weight_muon_reco_sf_SelectedMuon_central", "weight_muon_looseid_sf_SelectedMuon_central",
    "weight_photon_id_sf_SelectedPhoton_central", "weight_photon_csev_sf_SelectedPhoton_central",
]


def _sane_sf(a):
    """NaN/inf 或不合理 SF -> 1.0（空集合 SF 本應=1）。"""
    a = np.asarray(a, dtype=float)
    return np.where(np.isfinite(a) & (a > 0.01) & (a < 100.0), a, 1.0)
LUMI_MAP = {"2022preEE":7.9804,"2022postEE":26.6717,"2023preBPix":18.063,"2023postBPix":9.693,"2024":108.83}
DYLL_ERAS = ["2022preEE","2022postEE","2023preBPix","2023postBPix"]  # 有 DYJetsToLL inclusive 的 era
# 2024 沒有 DYJetsToLL inclusive，改用 split 樣本合併等同 inclusive（region mask 自然分流 ee/mumu）
DY_SPLIT  = ["DYJetsTo2E","DYJetsTo2Mu","DYJetsTo2Tau"]
ALL_ERAS  = DYLL_ERAS + ["2024"]

# ---- plot variables ----
PT_E   = [7,15,20,35,50,100,500]
PT_MU  = [5,10,15,20,25,30,40,50,60,120,500]
PT_MU2 = [10,15,20,25,30,40,50,60,120,500]
PT_PHO = [10,20,35,50,80]
ETA    = [round(-2.5+0.5*i,3) for i in range(11)]
MASSZ  = [80+i for i in range(0,21)]

# var: (key, column, title, bins, logx)
VARS_EE = [
    ("lead_electron_pt","Zee_lead_lepton_pt","Lead electron p_{T} [GeV]",PT_E,True),
    ("sublead_electron_pt","Zee_sublead_lepton_pt","Sublead electron p_{T} [GeV]",PT_E,True),
    ("lead_electron_eta","Zee_lead_lepton_eta","Lead electron #eta",ETA,False),
    ("mee","Zee_mass","m_{ee} [GeV]",MASSZ,False),
]
VARS_MU = [
    ("lead_muon_pt","zmmg_lead_muon_pt","Lead muon p_{T} [GeV]",PT_MU,True),
    ("sublead_muon_pt","zmmg_sublead_muon_pt","Sublead muon p_{T} [GeV]",PT_MU2,True),
    ("lead_muon_eta","zmmg_lead_muon_eta","Lead muon #eta",ETA,False),
    ("mmmg","Zmmg_mass","m_{#mu#mu#gamma} [GeV]",MASSZ,False),
]
VARS_PHO = [
    ("photon_pt","probe_photon_pt","Probe photon p_{T} [GeV]",PT_PHO,True),
    ("photon_eta","probe_photon_eta","Probe photon #eta",ETA,False),
]

# sf: (name, region, data_kind, weight_branch, after_label, before_label, vars)
SFGROUPS = [
    ("electron_reco_sf","z_ee","ee","weight_electron_reco_sf_SelectedElectron_central",
     "MC (Wi Electron reco SF)","MC (Wo Electron reco SF)",VARS_EE),
    ("electron_id_sf","z_ee","ee","weight_electron_wplid_sf_SelectedElectron_central",
     "MC (Wi Electron ID SF)","MC (Wo Electron ID SF)",VARS_EE),
    ("electron_trigger_sf","z_ee","ee","weight_hlt_sf_central",
     "MC (Wi Trigger SF)","MC (Wo Trigger SF)",VARS_EE),
    ("muon_reco_sf","z_mumu","mmg","weight_muon_reco_sf_SelectedMuon_central",
     "MC (Wi Muon reco SF)","MC (Wo Muon reco SF)",VARS_MU),
    ("muon_id_sf","z_mumu","mmg","weight_muon_looseid_sf_SelectedMuon_central",
     "MC (Wi Muon ID SF)","MC (Wo Muon ID SF)",VARS_MU),
    ("muon_trigger_sf","z_mumu","mmg","weight_hlt_sf_central",
     "MC (Wi Trigger SF)","MC (Wo Trigger SF)",VARS_MU),
    ("photon_id_sf","z_mumu","mmg","weight_photon_id_sf_SelectedPhoton_central",
     "MC (Wi Photon ID SF)","MC (Wo Photon ID SF)",VARS_PHO),
    ("photon_csev_sf","z_mumu","mmg","weight_photon_csev_sf_SelectedPhoton_central",
     "MC (Wi Photon CSEV SF)","MC (Wo Photon CSEV SF)",VARS_PHO),
]


def mc_paths(base, eras):
    """每個 era 優先用 DYJetsToLL inclusive；沒有（如 2024）則合併 split 樣本。"""
    out = []
    for e in eras:
        incl = glob.glob(f"{base}/mc/DYJetsToLL_{e}/merged_nominal.parquet")
        if incl:
            out += incl
        else:
            for s in DY_SPLIT:
                out += glob.glob(f"{base}/mc/{s}_{e}/merged_nominal.parquet")
    return out

def data_paths(base, eras, kind):
    sub = "data_egamma/Data_%s" if kind == "ee" else "data_muon/Data_tnp_zmmg_%s"
    return [p for e in eras for p in glob.glob(f"{base}/{sub % e}/merged_nominal.parquet")]


def stream_batches(paths, region, cols, batch_size=1000000):
    """逐檔 iter_batches 串流；每批套 region mask，yield (arrs, m)。

    只把單批（batch_size 列）載入記憶體，直方圖在呼叫端逐批累積，記憶體 O(bins)，
    與樣本大小無關。2024 split 樣本破億列且 z_ee/z_mumu 通過率極高（~1.4 億），
    若先把通過列全部累積成陣列會 OOM，故改為 streaming fill。
    """
    need = list(dict.fromkeys(cols + [region]))
    for p in paths:
        pf = pq.ParquetFile(p)
        avail = set(pf.schema_arrow.names)
        rd = [c for c in need if c in avail]
        if region not in rd:  # 沒有 region 欄無法遮罩，跳過該檔
            continue
        for batch in pf.iter_batches(columns=rd, batch_size=batch_size):
            arrs = {rd[j]: batch.column(j).to_numpy(zero_copy_only=False) for j in range(len(rd))}
            m = arrs[region].astype(bool)
            if m.any():
                yield arrs, m


def _acc_init(edges):
    nb = len(edges) - 1
    return [np.zeros(nb), np.zeros(nb), 0.0, 0.0]  # cnt, err2, overflow_lo, overflow_hi


def _acc_add(acc, vals, w, edges):
    """把一批 (vals,w) 累積進 acc；語意與 make_hist 相同（over/underflow 併入邊界 bin 內容）。"""
    vals = np.asarray(vals, dtype=float); w = np.asarray(w, dtype=float)
    good = np.isfinite(vals) & np.isfinite(w) & (vals > -900)
    vals, w = vals[good], w[good]
    e = np.asarray([float(x) for x in edges], dtype=float)
    inr = (vals >= e[0]) & (vals < e[-1])
    if inr.any():
        c, _ = np.histogram(vals[inr], bins=e, weights=w[inr]); acc[0] += c
        d, _ = np.histogram(vals[inr], bins=e, weights=w[inr] * w[inr]); acc[1] += d
    lo = vals < e[0]; hi = vals >= e[-1]
    if lo.any(): acc[2] += float(w[lo].sum())
    if hi.any(): acc[3] += float(w[hi].sum())


def _acc_hist(name, acc, edges):
    e = array("d", [float(x) for x in edges]); nb = len(e) - 1
    h = ROOT.TH1D(name, name, nb, e); h.Sumw2()
    for i in range(nb):
        h.SetBinContent(i + 1, acc[0][i]); h.SetBinError(i + 1, acc[1][i] ** 0.5)
    if acc[2]: h.SetBinContent(1, h.GetBinContent(1) + acc[2])
    if acc[3]: h.SetBinContent(nb, h.GetBinContent(nb) + acc[3])
    h.SetDirectory(0)
    return h


def make_hist(name, vals, w, edges):
    edges = array("d", [float(x) for x in edges])
    h = ROOT.TH1D(name, name, len(edges) - 1, edges)
    h.Sumw2()
    vals = np.asarray(vals, dtype=float); w = np.asarray(w, dtype=float)
    good = np.isfinite(vals) & np.isfinite(w) & (vals > -900)
    vals, w = vals[good], w[good]
    cnt, _ = np.histogram(vals, bins=np.asarray(edges), weights=w)
    er2, _ = np.histogram(vals, bins=np.asarray(edges), weights=w * w)
    # overflow/underflow into edge bins
    nb = len(edges) - 1
    for i in range(nb):
        h.SetBinContent(i + 1, cnt[i]); h.SetBinError(i + 1, er2[i] ** 0.5)
    lo = vals < edges[0]; hi = vals >= edges[-1]
    if lo.any():
        h.SetBinContent(1, h.GetBinContent(1) + w[lo].sum())
    if hi.any():
        h.SetBinContent(nb, h.GetBinContent(nb) + w[hi].sum())
    h.SetDirectory(0)
    return h


def scale_to_data(data, mc):
    di, mi = data.Integral(), mc.Integral()
    if di > 0 and mi > 0:
        mc.Scale(di / mi)


def cms_labels(pad, lumi):
    pad.cd()
    t = ROOT.TLatex(); t.SetNDC(True); t.SetTextFont(42)
    t.SetTextAlign(13); t.SetTextSize(0.060); t.DrawLatex(0.17, 0.96, "#bf{CMS} #it{Preliminary}")
    t.SetTextAlign(31); t.SetTextSize(0.055); t.DrawLatex(0.95, 0.915, f"{lumi:.2f} fb^{{-1}} (13.6 TeV)")


def draw(data, after, before, title, after_lab, before_lab, lumi, out, logx):
    for h in (data, after, before):
        h.SetTitle(""); h.SetStats(False)
    data.SetMarkerStyle(20); data.SetMarkerSize(1.3); data.SetLineColor(ROOT.kBlack); data.SetLineWidth(3); data.SetMarkerColor(ROOT.kBlack)
    after.SetLineColor(ROOT.kRed + 1); after.SetLineWidth(3); after.SetFillStyle(0)
    before.SetLineColor(ROOT.TColor.GetColor("#1F78B4")); before.SetLineWidth(3); before.SetLineStyle(2); before.SetFillStyle(0)

    c = ROOT.TCanvas("c", "", 800, 800)
    up = ROOT.TPad("up", "", 0, 0.30, 1, 1); lo = ROOT.TPad("lo", "", 0, 0, 1, 0.30)
    up.SetBottomMargin(0.03); up.SetLeftMargin(0.15); up.SetRightMargin(0.05)
    lo.SetTopMargin(0.04); lo.SetBottomMargin(0.45); lo.SetLeftMargin(0.15); lo.SetRightMargin(0.05)
    if logx:
        up.SetLogx(); lo.SetLogx()
    up.Draw(); lo.Draw()

    up.cd()
    ymax = max(data.GetMaximum(), after.GetMaximum(), before.GetMaximum())
    after.SetMinimum(0.0); after.SetMaximum(ymax * 1.5 if ymax > 0 else 1.0)
    after.GetXaxis().SetLabelSize(0); after.GetYaxis().SetTitle("Events")
    after.GetYaxis().SetTitleSize(0.065); after.GetYaxis().SetLabelSize(0.055); after.GetYaxis().SetTitleOffset(1.1)
    after.Draw("hist"); before.Draw("hist same"); data.Draw("E1 same")
    leg = ROOT.TLegend(0.18, 0.66, 0.55, 0.88); leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.045)
    leg.AddEntry(data, "Data", "lep"); leg.AddEntry(after, after_lab, "l"); leg.AddEntry(before, before_lab, "l"); leg.Draw()
    cms_labels(up, lumi)

    lo.cd()
    edges = array("d", [data.GetBinLowEdge(i + 1) for i in range(data.GetNbinsX() + 1)])
    frame = ROOT.TH1D("frame", "", data.GetNbinsX(), edges); frame.SetStats(False); frame.SetDirectory(0)
    frame.GetYaxis().SetTitle("Data / MC"); frame.GetXaxis().SetTitle(title)
    frame.GetYaxis().SetRangeUser(0.45, 1.55); frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.14); frame.GetYaxis().SetTitleOffset(0.45); frame.GetYaxis().SetLabelSize(0.12)
    frame.GetXaxis().SetTitleSize(0.14); frame.GetXaxis().SetTitleOffset(1.3); frame.GetXaxis().SetLabelSize(0.12)
    frame.Draw("axis")
    line = ROOT.TLine(edges[0], 1.0, edges[-1], 1.0); line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(2); line.Draw("same")
    keep = []
    for mc, col, mk in ((before, "#1F78B4", 21), (after, None, 20)):
        if mc.Integral() <= 0:
            continue
        r = data.Clone("r_" + mc.GetName()); r.SetDirectory(0); r.Divide(mc)
        r.SetMarkerStyle(mk); r.SetMarkerSize(1.2)
        c2 = ROOT.kRed + 1 if col is None else ROOT.TColor.GetColor(col)
        r.SetMarkerColor(c2); r.SetLineColor(c2); r.SetLineWidth(3); r.Draw("E1 same"); keep.append(r)

    c.SaveAs(out + ".pdf"); c.SaveAs(out + ".png"); c.Close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=DEFAULT_ROOT)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "plots", "SF_validation_zee_zmmg"))
    ap.add_argument("--eras", default="all", help="comma list or 'all' (+run3_all)")
    args = ap.parse_args()

    if args.eras.strip().lower() == "all":
        groups = [(e, [e]) for e in ALL_ERAS] + [("run3_all", ALL_ERAS)]
    else:
        groups = [(e, [e]) for e in args.eras.split(",") if e.strip()]

    os.makedirs(args.out, exist_ok=True)
    print(f"[cfg] base={args.base} out={args.out}")

    for label, eras in groups:
        lumi = sum(LUMI_MAP.get(e, 0.0) for e in eras)
        mcp = mc_paths(args.base, eras)
        if not mcp:
            print(f"[skip] {label}: no DYJetsToLL MC"); continue
        print(f"\n[era] {label} lumi={lumi:.2f} mc_files={len(mcp)}")
        for sfname, region, kind, wbranch, alab, blab, vlist in SFGROUPS:
            dpaths = data_paths(args.base, eras, kind)
            cols = [v[1] for v in vlist]
            # 逐批累積直方圖（記憶體 O(bins)）：acc_a=套SF, acc_b=不套此SF, acc_d=data
            acc_a = {v[0]: _acc_init(v[3]) for v in vlist}
            acc_b = {v[0]: _acc_init(v[3]) for v in vlist}
            acc_d = {v[0]: _acc_init(v[3]) for v in vlist}
            mc_n = da_n = 0
            # ---- MC：串流讀，逐批重建權重並填圖 ----
            for arrs, m in stream_batches(mcp, region, cols + ["weight_central_initial"] + SF_CENTRAL):
                mc_n += int(m.sum())
                # 自行重建 MC 權重：genweight × Π(各 SF，NaN→1)；繞過被 electron-SF NaN 毒化的 weight_central
                wci = arrs.get("weight_central_initial")
                wci = np.where(np.isfinite(wci[m]), wci[m], 0.0) if wci is not None else np.zeros(int(m.sum()))
                w_after = wci.copy()
                for sc in SF_CENTRAL:
                    v = arrs.get(sc)
                    w_after = w_after * _sane_sf(v[m] if v is not None else np.ones(int(m.sum())))
                vb = arrs.get(wbranch)
                w_before = w_after / _sane_sf(vb[m] if vb is not None else np.ones(int(m.sum())))
                for key, col, title, bins, logx in vlist:
                    v = arrs.get(col); v = v[m] if v is not None else np.zeros(int(m.sum()))
                    _acc_add(acc_a[key], v, w_after, bins)
                    _acc_add(acc_b[key], v, w_before, bins)
            # ---- Data：串流讀，單位權重填圖 ----
            for arrs, m in stream_batches(dpaths, region, cols + ["weight_central"]):
                n = int(m.sum()); da_n += n
                for key, col, title, bins, logx in vlist:
                    v = arrs.get(col); v = v[m] if v is not None else np.zeros(n)
                    _acc_add(acc_d[key], v, np.ones(n), bins)
            if mc_n == 0 or da_n == 0:
                print(f"  [skip] {sfname}: empty mc/data in {region}"); continue
            outdir = os.path.join(args.out, label, sfname); os.makedirs(outdir, exist_ok=True)
            for key, col, title, bins, logx in vlist:
                hd = _acc_hist(f"d_{label}_{sfname}_{key}", acc_d[key], bins)
                ha = _acc_hist(f"a_{label}_{sfname}_{key}", acc_a[key], bins)
                hb = _acc_hist(f"b_{label}_{sfname}_{key}", acc_b[key], bins)
                # Normalise MC-with-SF (after=ha) to data, then scale MC-without-SF
                # (before=hb) by the SAME factor. Normalising hb separately to data
                # (the old `scale_to_data(hd, hb)`) absorbs the SF's overall (few-%)
                # effect, leaving only the residual shape; sharing ha's factor keeps
                # the SF's full effect (normalisation + shape) visible on the plot.
                ai = ha.Integral()
                if ai > 0:
                    s = hd.Integral() / ai
                    ha.Scale(s); hb.Scale(s)
                draw(hd, ha, hb, title, alab, blab, lumi, os.path.join(outdir, key), logx)
            print(f"  [done] {sfname}: {len(vlist)} vars (mc={mc_n}, data={da_n})")
    print("\nAll done")


if __name__ == "__main__":
    main()
