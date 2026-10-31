#!/usr/bin/env python3
"""SF-validation plots for the zee_zmmg control sample.

讀 control-sample merged parquet（pyarrow，直接讀，不轉 root），畫 data/MC 的
「套 SF vs 不套 SF」比較（上圖分佈 + 下圖 Data/MC ratio），驗證 lepton / photon SF：
  - Z->ee   (data=EGamma) : electron reco / wplid(ID) SF
  - Z->mumugamma (data=Muon): muon reco / looseid SF、photon id / csev SF

MC 用 DYJetsToLL（與既有 plot_SF_validation.py 慣例一致；只有 2022/2023 4 個 era 有），
每張圖把 MC 縮放到 data（只比形狀，不需 lumi/xs）。
Env: higgs-alp-ana（ROOT 6.24 + pyarrow）。
輸出: /eos/home-p/pelai/HZa/output_plots/SF_validation_zee_zmmg/<era>/<sf>/<var>.{pdf,png}
每個 era 的直方圖存在 --cache-dir/<era>.pkl；只改版面時用 --from-cache 重畫（幾秒，不讀 parquet）。
各 era 可分開平行跑（--eras 2024），最後 --eras run3_all --from-cache 由各 era 相加。
"""
import argparse, os, glob, pickle, sys
from array import array
import numpy as np
import pyarrow.parquet as pq
# 直接呼叫 env 的 python（沒 activate）時 PATH 裡沒有 conda 的 c++，cling 啟動時會印
# "cannot extract standard library include paths"；把 python 所在的 bin 放到 PATH 最前面
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ.get("PATH", "")
import ROOT  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)
ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLegendBorderSize(0); ROOT.gStyle.SetErrorX(0.5); ROOT.gStyle.SetEndErrorSize(0)

DEFAULT_ROOT = "/eos/project/h/htozg-dy-privatemc/pelai/HZa/parquet_DNA_zee_zmmg_ctrl"
# 所有會乘進 weight_central 的 SF central 分支。
# 2026-10-30：control tagger 修好 SelectedElectron 為 None 的問題後，weight_central 不再是 NaN，
# 改用 weight_central（genweight × scale1fb × lumi × 全部 SF）。舊做法用 weight_central_initial（只有 genweight）自行乘 SF，
# 沒有截面正規化，無法把 DYJets 與 DYGto2LG 兩個樣本相加。
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
# lumi 從 HZgamma constants.py 讀（繪圖規範；HZa repo 那份會過期）
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/plot/scripts")
from lumi_constants import load_lumi  # noqa: E402
LUMI_MAP = load_lumi()
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


def mc_paths(base, eras, region):
    """每個 era 優先用 DYJetsToLL inclusive；沒有（如 2024）則合併 split 樣本。
    z_mumu（Z->mumugamma）區再加 DYGto2LG_10to100：2026-10-30 發現 DY+jets 在這一區幾乎沒有事件
    （2022postEE 一個 job 40898 列只有 2 列），FSR 光子事件在 DYGto2LG；control 沒有 overlap veto，重複計算可忽略。
    只用 DY+jets 時 MC 比 data 少約 20 倍，被「MC 縮放到 data」蓋掉，舊的 photon/muon SF 圖不可靠。"""
    out = []
    for e in eras:
        if region == "z_mumu":
            out += glob.glob(f"{base}/mc/DYGto2LG_10to100_{e}/merged_nominal.parquet")
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


def stream_batches(path, regions, cols, batch_size=1000000):
    """單檔 iter_batches 串流；每批對每個 region 做 mask，yield (arrs, {region: m})。

    只把單批（batch_size 列）載入記憶體，直方圖在呼叫端逐批累積，記憶體 O(bins)，
    與樣本大小無關。2024 split 樣本破億列且 z_ee/z_mumu 通過率極高（~1.4 億），
    若先把通過列全部累積成陣列會 OOM，故改為 streaming fill。
    一個檔同時屬於兩個 region 時（DYJets）只讀一遍。
    """
    pf = pq.ParquetFile(path)
    avail = set(pf.schema_arrow.names)
    regions = [r for r in regions if r in avail]  # 沒有 region 欄無法遮罩，跳過
    if not regions:
        return
    rd = [c for c in dict.fromkeys(cols + regions) if c in avail]
    for batch in pf.iter_batches(columns=rd, batch_size=batch_size):
        arrs = {rd[j]: batch.column(j).to_numpy(zero_copy_only=False) for j in range(len(rd))}
        masks = {r: arrs[r].astype(bool) for r in regions}
        masks = {r: m for r, m in masks.items() if m.any()}
        if masks:
            yield arrs, masks


REGION_KIND = {"z_ee": "ee", "z_mumu": "mmg"}  # region -> data stream（SFGROUPS 裡每個 region 只對應一種 data）


def region_vars(region):
    out = {}
    for g in SFGROUPS:
        if g[1] == region:
            for v in g[6]:
                out.setdefault(v[0], v)
    return list(out.values())


def era_inputs(base, era):
    """{path: [regions]}（MC）與 {path: [region]}（data）。"""
    mc, data = {}, {}
    for region, kind in REGION_KIND.items():
        for p in mc_paths(base, [era], region):
            mc.setdefault(p, []).append(region)
        for p in data_paths(base, [era], kind):
            data.setdefault(p, []).append(region)
    return mc, data


def fingerprint(paths):
    return {p: (os.path.getsize(p), int(os.path.getmtime(p))) for p in sorted(paths)}


def fill_era(base, era):
    """一個 era 的所有直方圖，每個輸入檔只讀一遍、8 個 SF 組同時填。

    舊做法每個 SF 組各讀一次 data 與 MC，run3_all 再把五個 era 全部重讀，同一個檔被讀 ~16 次（全量 3.7 小時）。
    acc["a"][region][var]：MC 套全部 SF；acc["d"][region][var]：data；acc["b"][sfname][var]：MC 拿掉該 SF。
    """
    mc, data = era_inputs(base, era)
    acc = {"a": {r: {v[0]: _acc_init(v[3]) for v in region_vars(r)} for r in REGION_KIND},
           "d": {r: {v[0]: _acc_init(v[3]) for v in region_vars(r)} for r in REGION_KIND},
           "b": {g[0]: {v[0]: _acc_init(v[3]) for v in g[6]} for g in SFGROUPS},
           "n_mc": {r: 0 for r in REGION_KIND}, "n_data": {r: 0 for r in REGION_KIND},
           "inputs": fingerprint(list(mc) + list(data))}
    for p, regions in mc.items():
        groups = [g for g in SFGROUPS if g[1] in regions]
        cols = [v[1] for r in regions for v in region_vars(r)] + ["weight_central"] + [g[3] for g in groups]
        for arrs, masks in stream_batches(p, regions, cols):
            for region, m in masks.items():
                n = int(m.sum()); acc["n_mc"][region] += n
                # weight_central = genweight × scale1fb × lumi × 全部 SF（樣本間可相加）
                wc = arrs.get("weight_central")
                w_after = np.where(np.isfinite(wc[m]), wc[m], 0.0) if wc is not None else np.zeros(n)
                vals = {}
                for key, col, title, bins, logx in region_vars(region):
                    v = arrs.get(col); vals[key] = v[m] if v is not None else np.zeros(n)
                    _acc_add(acc["a"][region][key], vals[key], w_after, bins)
                for g in groups:
                    if g[1] != region:
                        continue
                    vb = arrs.get(g[3])
                    w_before = w_after / _sane_sf(vb[m] if vb is not None else np.ones(n))
                    for key, col, title, bins, logx in g[6]:
                        _acc_add(acc["b"][g[0]][key], vals[key], w_before, bins)
        print(f"  [read] {era} {os.path.basename(os.path.dirname(p))} n_mc={acc['n_mc']}", flush=True)
    for p, regions in data.items():
        cols = [v[1] for r in regions for v in region_vars(r)]
        for arrs, masks in stream_batches(p, regions, cols):
            for region, m in masks.items():
                n = int(m.sum()); acc["n_data"][region] += n
                for key, col, title, bins, logx in region_vars(region):
                    v = arrs.get(col); v = v[m] if v is not None else np.zeros(n)
                    _acc_add(acc["d"][region][key], v, np.ones(n), bins)
        print(f"  [read] {era} {os.path.basename(os.path.dirname(p))} n_data={acc['n_data']}", flush=True)
    return acc


def sum_accs(accs):
    """直方圖可直接相加（MC 縮放到 data 在相加之後才做），run3_all 不必重讀檔案。"""
    def add(x, y):
        if isinstance(x, dict):
            return {k: add(x[k], y[k]) for k in x}
        if isinstance(x, list):
            return [add(a, b) for a, b in zip(x, y)]
        return x + y
    out = None
    for a in accs:
        a = {k: v for k, v in a.items() if k != "inputs"}
        out = a if out is None else add(out, a)
    return out


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

    # ratio 先算（未縮放的直方圖），y 範圍才能依資料點決定，不會把高 pT 的點畫出框外
    ratios = []
    for mc, col, mk in ((before, "#1F78B4", 21), (after, None, 20)):
        if mc.Integral() <= 0:
            continue
        r = data.Clone("r_" + mc.GetName()); r.SetDirectory(0); r.Divide(mc)
        r.SetMarkerStyle(mk); r.SetMarkerSize(1.2)
        c2 = ROOT.kRed + 1 if col is None else ROOT.TColor.GetColor(col)
        r.SetMarkerColor(c2); r.SetLineColor(c2); r.SetLineWidth(3); ratios.append(r)
    pts = [(r.GetBinContent(i) - r.GetBinError(i), r.GetBinContent(i) + r.GetBinError(i))
           for r in ratios for i in range(1, r.GetNbinsX() + 1) if r.GetBinContent(i) > 0]
    rlo = max(0.0, min([0.45] + [lo_ - 0.1 for lo_, _ in pts]))
    rhi = min(4.0, max([1.55] + [hi_ + 0.1 for _, hi_ in pts]))

    up.cd()
    ymax = max(data.GetMaximum(), after.GetMaximum(), before.GetMaximum())
    # 計數 >= 1e4 時改寫成 Events / 10^{k}：ROOT 自動加的「×10^{k}」落在框頂左角，會和 CMS 標籤疊字
    k = 3 * int(np.floor(np.log10(ymax * 1.5) / 3)) if ymax * 1.5 >= 1e4 else 0
    if k:
        for h in (data, after, before):
            h.Scale(10.0 ** -k)
        ymax *= 10.0 ** -k
    after.SetMinimum(0.0); after.SetMaximum(ymax * 1.5 if ymax > 0 else 1.0)
    after.GetXaxis().SetLabelSize(0); after.GetYaxis().SetTitle("Events / 10^{%d}" % k if k else "Events")
    after.GetYaxis().SetTitleSize(0.065); after.GetYaxis().SetLabelSize(0.055); after.GetYaxis().SetTitleOffset(1.1)
    after.Draw("hist"); before.Draw("hist same"); data.Draw("E1 same")
    leg = ROOT.TLegend(0.18, 0.66, 0.55, 0.88); leg.SetFillStyle(0); leg.SetBorderSize(0); leg.SetTextSize(0.045)
    leg.AddEntry(data, "Data", "lep"); leg.AddEntry(after, after_lab, "l"); leg.AddEntry(before, before_lab, "l"); leg.Draw()
    cms_labels(up, lumi)

    lo.cd()
    edges = array("d", [data.GetBinLowEdge(i + 1) for i in range(data.GetNbinsX() + 1)])
    frame = ROOT.TH1D("frame", "", data.GetNbinsX(), edges); frame.SetStats(False); frame.SetDirectory(0)
    frame.GetYaxis().SetTitle("Data / MC"); frame.GetXaxis().SetTitle(title)
    frame.GetYaxis().SetRangeUser(rlo, rhi); frame.GetYaxis().SetNdivisions(505)
    frame.GetYaxis().SetTitleSize(0.14); frame.GetYaxis().SetTitleOffset(0.45); frame.GetYaxis().SetLabelSize(0.12)
    frame.GetXaxis().SetTitleSize(0.14); frame.GetXaxis().SetTitleOffset(1.3); frame.GetXaxis().SetLabelSize(0.12)
    if logx:
        # 跨不到 20 倍（photon pT 10-80）時 log 軸一個十進位刻度都沒有，要補中間刻度的數字
        frame.GetXaxis().SetMoreLogLabels(edges[-1] / edges[0] < 20); frame.GetXaxis().SetNoExponent(True)
    frame.Draw("axis")
    line = ROOT.TLine(edges[0], 1.0, edges[-1], 1.0); line.SetLineColor(ROOT.kGray + 2); line.SetLineStyle(2); line.Draw("same")
    for r in ratios:
        r.Draw("E1 same")

    c.SaveAs(out + ".pdf"); c.SaveAs(out + ".png"); c.Close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=DEFAULT_ROOT)
    ap.add_argument("--out", default="/eos/home-p/pelai/HZa/output_plots/SF_validation_zee_zmmg")
    ap.add_argument("--eras", default="all",
                    help="comma list of eras and/or run3_all (= sum of the 5 eras); 'all' = 5 eras + run3_all")
    ap.add_argument("--cache-dir", default="/eos/home-p/pelai/HZa/output_tables/SF_validation_zee_zmmg/hist_cache",
                    help="per-era accumulated histograms (<era>.pkl), written after every fill")
    ap.add_argument("--from-cache", action="store_true",
                    help="do not read parquet; use --cache-dir (refused if an input file changed since the cache was written)")
    args = ap.parse_args()

    labels = ALL_ERAS + ["run3_all"] if args.eras.strip().lower() == "all" else [e.strip() for e in args.eras.split(",") if e.strip()]
    need = list(dict.fromkeys(e for l in labels for e in (ALL_ERAS if l == "run3_all" else [l])))

    os.makedirs(args.out, exist_ok=True); os.makedirs(args.cache_dir, exist_ok=True)
    print(f"[cfg] base={args.base} out={args.out} cache={args.cache_dir} from_cache={args.from_cache}", flush=True)

    accs = {}
    for era in need:
        cf = os.path.join(args.cache_dir, f"{era}.pkl")
        mc, data = era_inputs(args.base, era)
        if args.from_cache:
            if not os.path.exists(cf):
                sys.exit(f"[error] no cache for {era}: {cf}")
            with open(cf, "rb") as f:
                accs[era] = pickle.load(f)
            # 快取只存直方圖，不知道輸入是否換過：比對檔案大小與 mtime，換過就拒絕
            if accs[era]["inputs"] != fingerprint(list(mc) + list(data)):
                sys.exit(f"[error] inputs of {era} changed since {cf} was written; rerun without --from-cache")
            print(f"[cache] {era} <- {cf}", flush=True)
        else:
            accs[era] = fill_era(args.base, era)
            with open(cf + ".tmp", "wb") as f:
                pickle.dump(accs[era], f)
            os.replace(cf + ".tmp", cf)
            print(f"[cache] {era} -> {cf}", flush=True)

    for label in labels:
        eras = ALL_ERAS if label == "run3_all" else [label]
        lumi = sum(LUMI_MAP.get(e, 0.0) for e in eras)
        if not any(mc_paths(args.base, [e], "z_ee") for e in eras):
            print(f"[skip] {label}: no DYJetsToLL MC"); continue
        acc = sum_accs([accs[e] for e in eras])
        print(f"\n[era] {label} lumi={lumi:.2f}")
        for sfname, region, kind, wbranch, alab, blab, vlist in SFGROUPS:
            mc_n, da_n = acc["n_mc"][region], acc["n_data"][region]
            if mc_n == 0 or da_n == 0:
                print(f"  [skip] {sfname}: empty mc/data in {region}"); continue
            outdir = os.path.join(args.out, label, sfname); os.makedirs(outdir, exist_ok=True)
            for key, col, title, bins, logx in vlist:
                hd = _acc_hist(f"d_{label}_{sfname}_{key}", acc["d"][region][key], bins)
                ha = _acc_hist(f"a_{label}_{sfname}_{key}", acc["a"][region][key], bins)
                hb = _acc_hist(f"b_{label}_{sfname}_{key}", acc["b"][sfname][key], bins)
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
