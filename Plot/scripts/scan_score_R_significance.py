#!/usr/bin/env python3
"""Low-mass (mA=1,2,3) MVA working-point scan with the production 5-variable BDT.

For each mA, scan the BDT score cut and compute:
  - R   : 125-peak sculpting metric = (selected-bkg fraction in 120-130) / inclusive fraction.
          R=1 -> no sculpting (peak neither enriched nor depleted); R>1 -> fake-peak risk;
          R<1 -> peak depleted (conservative but wastes significance).
  - Z   : Asimov-like significance sqrt(2((s+b)ln(1+s/b))-2s) in the 120-130 window.

Selection philosophy (Pei-Zhu): pick the cut with R closest to 1 (no sculpting) rather than
maximum significance, to avoid reporting a sculpted/fake result. R=1 maximizes Z subject to
no peak enrichment.

Outputs a 2D plot (x=BDT score, y=R, colour=Z) per mA to
  Plot/plots/scan_score_R/lowmass_score_R_significance.{pdf,png}
and prints the R=1 working point (cut + Z) per mA.

Run in higgs-alp-ana (xgboost + uproot + ROOT 6.24).
"""
import os, math, json, shutil, argparse, datetime, numpy as np, pickle, uproot, ROOT


def _significance_asimov_like(s: float, b: float) -> float:
    """開根號版本：sqrt( (2*(s+b)*log(1+s/b)) - 2*s )"""
    if s <= 0.0:
        return 0.0
    if b <= 0.0:
        return 0.0
    val = (2.0 * (s + b) * math.log(1.0 + (s / b))) - (2.0 * s)
    # 數值誤差保護：避免極小負值導致 sqrt 出錯
    if val <= 0.0:
        return 0.0
    return math.sqrt(val)

ROOT_DIR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal"
MODEL    = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_lowmass_run3.pkl"
HIGHMASS_MODEL = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_highmass_run3.pkl"
FEATS    = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm"]   # + param
OUTDIR   = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/plots/scan_score_R"
JSON_PATH = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"
# Combined sculpting-ratio R(m_a) for ALL mass points (low-mass + high-mass BDT models),
# written for the AN table. Sourced from each model's .meta.json "R" field (R at the
# model's working point, evaluated on inclusive simulated background at training time).
SCULPT_R_JSON  = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/sculpt_R_run3.json"
LOWMASS_META   = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_lowmass_run3.meta.json"
HIGHMASS_META  = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/model_Za_BDT_highmass_run3.meta.json"
LOWM     = [1, 2, 3]
PEAK     = (120.0, 130.0)
SCORE_RANGE = (0.70, 0.99)   # extended so the chosen high-score working points (mA2 0.98, mA3 0.988) are on-plot

# --- stored-score inputs: regenerate the low-mass figure from the SAME BDT scores as the AN
# sculpt-R table (make_sculpt_R_table_from_wp.py), i.e. the analysis's stored MVA_Score_mA_M{n},
# so the figure R and the table R are computed identically (no re-scoring with the .pkl). ------
SCORED_BASE = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal"
_BKG_BY_YEAR = {
    "2022preEE":   ["DYGto2LG_10to50", "DYGto2LG_50to100", "DYJetsToLL"],
    "2022postEE":  ["DYGto2LG_10to50", "DYGto2LG_50to100", "DYJetsToLL"],
    "2023preBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
    "2023postBPix": ["DYGto2LG_10to100", "DYJetsToLL"],
    "2024":        ["DYGto2LG_10to100", "DYJetsTo2E", "DYJetsTo2Mu", "DYJetsTo2Tau"],
}
_ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]

def _load_scored(paths, mA):
    """Concatenate H_mass, factor and the stored score MVA_Score_mA_M{mA} over the given scored
    files; apply the 95<m<180 window. Uses the analysis's own BDT score (no re-scoring)."""
    sbr = f"MVA_Score_mA_M{mA}"
    H, W, Sc = [], [], []
    for p in paths:
        if not os.path.exists(p):
            continue
        a = uproot.open(p)["inclusive"].arrays(["H_mass", "factor", sbr], library="np")
        H.append(a["H_mass"]); W.append(a["factor"]); Sc.append(a[sbr])
    H = np.concatenate(H); W = np.concatenate(W); Sc = np.concatenate(Sc)
    m = (H > 95) & (H < 180) & np.isfinite(W) & np.isfinite(Sc)
    return {"H_mass": H[m], "factor": W[m], "s": Sc[m]}

def load_scored_bkg(mA):
    paths = [f"{SCORED_BASE}/{s}/{y}.root" for y, ss in _BKG_BY_YEAR.items() for s in ss]
    return _load_scored(paths, mA)

def load_scored_sig(mA):
    return _load_scored([f"{SCORED_BASE}/mA_M{mA}/{y}.root" for y in _ERAS], mA)

def _R_direct(bk, cut):
    """R at a fixed cut on the stored-score background (identical definition to the AN table:
    RAW weights, negatives kept, so the red star equals make_sculpt_R_table_from_wp.py)."""
    w = bk["factor"]; H = bk["H_mass"]; s = bk["s"]
    peak = (H > PEAK[0]) & (H < PEAK[1]); frac = w[peak].sum() / w.sum()
    pB = s > cut; denom = w[pB].sum()
    return float(w[pB & peak].sum() / denom / frac) if denom > 0 and frac > 0 else float("nan")
LUMI_LABEL  = "172.13 fb^{-1} (13.6 TeV)"

# ---------------------------------------------------------------------------------------------
# FIXED low-mass working points (HARDCODED -- override the bare R=1 crossing).  [Pei-Zhu 2026-06-17]
#   mA2 and mA3 do NOT use the plain R=1 crossing: at that (low-score) cut the selection keeps
#   too much data and the background-model goodness-of-fit FAILS (asymptotic-GOF regime). We
#   deliberately picked a slightly HIGHER BDT score near R=1, where N_data < ~425 so the toy-GOF
#   passes while R stays ~1.
#   mA2: cut 0.982, N_data~334, R~0.99  [Pei-Zhu 2026-06-19 -- moved RIGHT from 0.98 (N~394) to
#        suppress a 2.76 sigma excess at the old cut; 0.982 lands on the R=1 crossing, Z~18.4].
#   mA3: cut 0.988, N~364, R~1.05.
#   These are pinned HERE so that re-running this script -- in particular with --write-json --
#   can NOT silently revert mA2/mA3 back to the R=1 crossing. mA1 is NOT fixed: it keeps its
#   natural R=1 crossing (~0.963). See memory ref_hza_unblind_procedure / project_hza_lowma_pow1_R1cuts.
FIXED_WP = {2: 0.982, 3: 0.988}   # mA -> chosen BDT score cut (used for the red star AND --write-json)

import sys
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
from hza_features import BASE_VARS, FEATURES   # FEATURES == BASE_VARS + ['pho_pt_asym','param']

mdl      = pickle.load(open(MODEL, "rb"))            # low-mass model (mA1,2,3)
mdl_high = pickle.load(open(HIGHMASS_MODEL, "rb"))   # high-mass model (mA4-30; used here for mA30)
_HM_FEATS = json.load(open(HIGHMASS_META))["features"]   # high-mass model feature order (incl 'param')

def load(path, feats=FEATS):
    a = uproot.open(path)["inclusive"].arrays(list(feats) + ["H_mass", "ALP_mass", "factor"], library="np")
    m = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"])
    return {k: a[k][m] for k in a}

def score(a, mh):
    X = np.column_stack([a[c] for c in FEATS] + [(a["ALP_mass"] - mh) / a["H_mass"]])
    return mdl.predict_proba(X)[:, 1]

def score_high(a, mA):
    """BDT score from the high-mass parametric model. Builds the full feature matrix with the
    shared hza_features layer (BASE_VARS + derived pho_pt_asym + param), then selects the model's
    feature order -- same construction as _bkg_R_at_cuts(). 'a' must be loaded with BASE_VARS."""
    idx = [FEATURES.index(f) for f in _HM_FEATS]
    asym = (a["pho1Pt_oHm"] - a["pho2Pt_oHm"]) / (a["pho1Pt_oHm"] + a["pho2Pt_oHm"] + 1e-6)
    param = (a["ALP_mass"] - mA) / a["H_mass"]
    X = np.column_stack([a[b] for b in BASE_VARS] + [asym, param])[:, idx]
    return mdl_high.predict_proba(X)[:, 1]

def scan(mA):
    """R / Asimov-Z scan vs BDT cut for one low-mass hypothesis, using the STORED BDT scores
    (MVA_Score_mA_M{mA}) on the combined DY background and the mA signal. Returns (pts, bk)."""
    bk = load_scored_bkg(mA); sg = load_scored_sig(mA)
    wB = bk["factor"]; Hb = bk["H_mass"]; sb = bk["s"]   # RAW weights (match the AN sculpt-R table)
    peak = (Hb > PEAK[0]) & (Hb < PEAK[1]); frac = wB[peak].sum() / wB.sum()
    wS = np.clip(sg["factor"], 0, None); Hs = sg["H_mass"]; ss = sg["s"]
    pks = (Hs > PEAK[0]) & (Hs < PEAK[1])
    pts = []
    for thr in np.linspace(*SCORE_RANGE, 45):
        pB = sb > thr; nB = wB[pB & peak].sum()
        if nB <= 0 or (pB & peak).sum() < 20:   # too few bkg events -> R unreliable
            continue
        R = nB / wB[pB].sum() / frac
        S = wS[(ss > thr) & pks].sum()
        Z = _significance_asimov_like(S, nB)        # Asimov-like significance, S & B in 120-130
        pts.append((thr, R, Z))
    return np.array(pts), bk

def scan_high(mA):
    """Same R / Asimov-Z scan as scan() but with the HIGH-MASS model (for the mA30 panel).
    Self-contained: loads All_Bkg and the mA signal with BASE_VARS and scores via score_high()."""
    bk = load(f"{ROOT_DIR}/All_Bkg/run3.root", BASE_VARS)
    wB = np.clip(bk["factor"], 0, None); Hb = bk["H_mass"]
    peak = (Hb > PEAK[0]) & (Hb < PEAK[1]); frac = wB[peak].sum() / wB.sum()
    sg = load(f"{ROOT_DIR}/mA_M{mA}/run3.root", BASE_VARS)
    wS = np.clip(sg["factor"], 0, None); Hs = sg["H_mass"]; pks = (Hs > PEAK[0]) & (Hs < PEAK[1])
    sb, ss = score_high(bk, mA), score_high(sg, mA)
    pts = []
    for thr in np.linspace(*SCORE_RANGE, 45):
        pB = sb > thr; nB = wB[pB & peak].sum()
        if nB <= 0 or (pB & peak).sum() < 20:
            continue
        R = nB / wB[pB].sum() / frac
        S = wS[(ss > thr) & pks].sum()
        Z = _significance_asimov_like(S, nB)
        pts.append((thr, R, Z))
    return np.array(pts)

def r_equals_one(pts):
    """Linear-interpolate the (tightest) score where R crosses 1; return (score, R, Z)."""
    thr, R, Z = pts[:, 0], pts[:, 1], pts[:, 2]
    best = None
    for i in range(len(R) - 1):
        if (R[i] - 1) * (R[i + 1] - 1) <= 0 and R[i] != R[i + 1]:
            f = (1 - R[i]) / (R[i + 1] - R[i])
            best = (thr[i] + f * (thr[i + 1] - thr[i]), 1.0, Z[i] + f * (Z[i + 1] - Z[i]))
    if best is None:   # no crossing -> point closest to R=1
        j = int(np.argmin(np.abs(R - 1))); best = (thr[j], R[j], Z[j])
    return best

def chosen_wp(mA, pts, bk=None):
    """Final working point for this mA (red star). Priority for the cut: FIXED_WP (pinned
    mA2/mA3) > MVAcut in JSON_PATH > R=1 crossing. If a stored-score background `bk` is given
    (low-mass path) R is computed DIRECTLY at that cut so the red star equals the AN sculpt-R
    table value; otherwise (high-mass mA30 panel) R is interpolated from the scan curve. Z is
    always interpolated from the scan curve."""
    cut = None
    if mA in FIXED_WP:
        cut = float(FIXED_WP[mA])               # hardcoded: never reverts to the R=1 crossing
    else:
        try:
            with open(JSON_PATH) as f:
                doc = json.load(f)
            cut = next(float(e["MVAcut"]) for e in doc["results"]
                       if int(e["mA"]) == mA and e.get("MVAcut") is not None)
        except Exception:
            return r_equals_one(pts)
    R = _R_direct(bk, cut) if bk is not None else float(np.interp(cut, pts[:, 0], pts[:, 1]))
    Z = float(np.interp(cut, pts[:, 0], pts[:, 2]))
    return (cut, R, Z)

keep = []   # keep ROOT objects alive

def draw_one(mA, pts, r1, wp, big):
    """Draw one mA panel on the current pad. big=True -> larger axis/label/text for stand-alone."""
    Rlo, Rhi = pts[:, 1].min(), pts[:, 1].max(); pad = (Rhi - Rlo) * 0.12 + 1e-3
    xpad = (SCORE_RANGE[1] - SCORE_RANGE[0]) * 0.04   # small x-margin so edge points don't touch the frame
    xlo, xhi = SCORE_RANGE[0] - xpad, SCORE_RANGE[1] + xpad
    h = ROOT.TH2D(f"h{mA}_{int(big)}", ";BDT score cut;R (125-peak sculpting)",
                  45, xlo, xhi, 60, Rlo - pad, Rhi + pad)
    for bx in range(1, h.GetNbinsX() + 1):
        for by in range(1, h.GetNbinsY() + 1):
            h.SetBinContent(bx, by, -1e9)
    for thr, R, Z in pts:
        h.SetBinContent(h.GetXaxis().FindBin(thr), h.GetYaxis().FindBin(R), Z)
    zmin, zmax = pts[:, 2].min(), pts[:, 2].max()
    h.SetMinimum(zmin); h.SetMaximum(zmax); h.GetZaxis().SetTitle("Asimov significance Z")
    ts, ls = (0.055, 0.05) if big else (0.045, 0.04)
    for ax in (h.GetXaxis(), h.GetYaxis(), h.GetZaxis()):
        ax.SetTitleSize(ts); ax.SetLabelSize(ls)
    h.GetXaxis().SetTitleOffset(1.05 if big else 1.1); h.GetYaxis().SetTitleOffset(1.2 if big else 1.3)
    h.GetZaxis().SetTitleOffset(1.2 if big else 1.0)
    ROOT.gPad.SetRightMargin(0.19 if big else 0.16); ROOT.gPad.SetLeftMargin(0.16 if big else 0.13)
    ROOT.gPad.SetBottomMargin(0.14 if big else 0.12); ROOT.gPad.SetTopMargin(0.09)
    h.Draw("COLZ"); keep.append(h)
    mk = []
    for thr, R, Z in pts:
        col = ROOT.TColor.GetColorPalette(int(253 * (Z - zmin) / (zmax - zmin + 1e-9)))
        m = ROOT.TMarker(thr, R, 20); m.SetMarkerColor(col); m.SetMarkerSize(2.0 if big else 1.3); m.Draw(); mk.append(m)
    line = ROOT.TLine(xlo, 1.0, xhi, 1.0); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray + 2); line.SetLineWidth(2); line.Draw(); mk.append(line)
    wp_mk = ROOT.TMarker(wp[0], wp[1], 29); wp_mk.SetMarkerColor(ROOT.kRed); wp_mk.SetMarkerSize(3.3 if big else 2.6); wp_mk.Draw(); mk.append(wp_mk)
    keep.append(mk)
    tl = ROOT.TLatex(); tl.SetNDC(); tl.SetTextFont(42)
    tl.SetTextSize(0.042 if big else 0.038); tl.SetTextAlign(11); tl.DrawLatex(0.16 if big else 0.13, 0.935, "#bf{CMS} #it{Preliminary}")
    tl.SetTextSize(0.034 if big else 0.032); tl.SetTextAlign(31); tl.DrawLatex(0.81 if big else 0.84, 0.935, LUMI_LABEL)
    # m_a label inside the frame, top-left corner
    tl.SetTextAlign(13); tl.SetTextColor(ROOT.kBlack); tl.SetTextSize(0.05 if big else 0.04)
    tl.DrawLatex(0.185 if big else 0.155, 0.88, f"m_{{a}} = {mA} GeV")
    tl.SetTextSize(0.040 if big else 0.030); tl.SetTextColor(ROOT.kRed + 1)
    # Per-mA corner placement so the wrapped label sits in an empty region of each panel:
    #   mA1 bottom-left, mA2 bottom-right, mA3 top-left.  align = 10*h + v (h:1=L,3=R; v:1=bottom,3=top)
    #   mA30: R is always >1.2 (no R=1 crossing); the point band fills the top, so the label goes
    #   bottom-left where the panel is empty.
    wp_pos = {1: (0.20, 0.22, 11), 2: (0.78, 0.22, 31), 3: (0.20, 0.83, 13), 30: (0.20, 0.24, 11)}
    wx, wy, wa = wp_pos.get(mA, (0.20, 0.22, 11))
    tl.SetTextAlign(wa)
    tl.DrawLatex(wx, wy, f"#splitline{{Working point: cut>{wp[0]:.3f},}}{{R={wp[1]:.2f}, Z={wp[2]:.1f}}}")


def write_json(json_path, r1_by_mA):
    """Overwrite MVAcut (and matching Asimov Significance) for mA in r1_by_mA, keep the rest.

    r1_by_mA: {mA: (cut, R, Z)}. A timestamped .bak is written before overwriting.
    """
    with open(json_path) as f:
        doc = json.load(f)
    bak = json_path + ".bak_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(json_path, bak)
    changed = []
    for entry in doc["results"]:
        mA = int(entry["mA"])
        if mA in r1_by_mA:
            cut, _R, Z = r1_by_mA[mA]
            old = entry.get("MVAcut")
            entry["MVAcut"] = round(float(cut), 4)
            entry["Significance"] = round(float(Z), 3)
            changed.append((mA, old, entry["MVAcut"]))
    with open(json_path, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")
    print(f"[write-json] backup -> {bak}")
    for mA, old, new in changed:
        print(f"[write-json] mA={mA}: MVAcut {old} -> {new}")
    print(f"[write-json] updated {json_path}")


def _bkg_R_at_cuts(model_path, meta_path, cuts):
    """R(m_a) at a fixed BDT-score cut, per mA, for one parametric model, on the inclusive
    simulated background. R = f_peak^pass / f_peak^incl (same definition as scan()). Builds the
    full feature matrix with the shared hza_features layer (BASE_VARS + derived pho_pt_asym +
    param) so the high-mass derived feature is handled, then selects the model's feature order.
    NB: the 'inclusive' tree uses branch names H_mass/ALP_mass (not the train-tree H_m/ALP_m)."""
    import sys
    sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
    from hza_features import BASE_VARS, FEATURES   # FEATURES == BASE_VARS + ['pho_pt_asym','param']

    mdl = pickle.load(open(model_path, "rb"))
    feats = json.load(open(meta_path))["features"]          # model feature order (incl 'param')
    idx = [FEATURES.index(f) for f in feats]
    a = uproot.open(f"{ROOT_DIR}/All_Bkg/run3.root")["inclusive"].arrays(
        BASE_VARS + ["H_mass", "ALP_mass", "factor"], library="np")
    m = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"])
    a = {k: a[k][m] for k in a}
    wB = np.clip(a["factor"], 0, None); Hb = a["H_mass"]
    peak = (Hb > PEAK[0]) & (Hb < PEAK[1]); frac = wB[peak].sum() / wB.sum()
    asym = (a["pho1Pt_oHm"] - a["pho2Pt_oHm"]) / (a["pho1Pt_oHm"] + a["pho2Pt_oHm"] + 1e-6)
    base_cols = [a[b] for b in BASE_VARS]
    R = {}
    for mA, cut in cuts.items():
        param = (a["ALP_mass"] - mA) / a["H_mass"]
        X = np.column_stack(base_cols + [asym, param])[:, idx]  # FEATURES order -> model order
        s = mdl.predict_proba(X)[:, 1]
        pB = s > cut; denom = wB[pB].sum()
        R[mA] = float(wB[pB & peak].sum() / denom / frac) if denom > 0 else float("nan")
    return R


def write_sculpt_R_json():
    """Write SCULPT_R_JSON with R(m_a) for ALL mass points, RECOMPUTED at the FINAL working
    point (the MVAcut in JSON_PATH) on the inclusive simulated background: low-mass model for
    mA 1-3 (R=1-tuned cuts), high-mass model for mA 4-30. Consumed by make_sculpt_R_table.py."""
    cuts = {}
    for e in json.load(open(JSON_PATH))["results"]:
        if e.get("MVAcut") is not None:
            cuts[int(e["mA"])] = float(e["MVAcut"])
    low_masses  = [int(x) for x in json.load(open(LOWMASS_META))["masses"]]
    high_masses = [int(x) for x in json.load(open(HIGHMASS_META))["masses"]]
    low_R  = _bkg_R_at_cuts(MODEL,          LOWMASS_META,  {m: cuts[m] for m in low_masses  if m in cuts})
    high_R = _bkg_R_at_cuts(HIGHMASS_MODEL, HIGHMASS_META, {m: cuts[m] for m in high_masses if m in cuts})
    doc = {
        "low_mass":  {str(m): low_R[m]  for m in sorted(low_R)},
        "high_mass": {str(m): high_R[m] for m in sorted(high_R)},
    }
    with open(SCULPT_R_JSON, "w") as f:
        json.dump(doc, f, indent=2); f.write("\n")
    print(f"[sculpt-R] wrote {SCULPT_R_JSON} (R recomputed at WP cuts; low {sorted(low_R)}, high {sorted(high_R)})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--write-json", action="store_true",
                    help=f"overwrite mA1,2,3 MVAcut (R=1 working point) in {JSON_PATH}")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    # NB: the AN sculpt-R table (tab:sculpt_R) is now produced by make_sculpt_R_table_from_wp.py
    # from the STORED scores; write_sculpt_R_json() (which re-scores with the .pkl) is no longer
    # used for the table and is left disabled to avoid writing inconsistent values.
    # write_sculpt_R_json()
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.gStyle.SetPalette(ROOT.kBird)

    results = {}
    print(f"{'mA':>3}  {'R=1 cut':>9}  {'Z@R=1':>7}   {'R-min cut':>9}  {'R-min':>6}  {'Z@Rmin':>7}   {'WP cut':>8}  {'R@WP':>6}  {'Z@WP':>6}")
    for mA in LOWM:
        pts, bk = scan(mA)   # low-mass scan on the STORED BDT scores (matches the AN table)
        r1 = r_equals_one(pts); j = int(np.argmin(pts[:, 1]))
        wp = chosen_wp(mA, pts, bk)   # red star = R directly at the JSON MVAcut (matches the AN table)
        print(f"{mA:>3}  {r1[0]:>9.3f}  {r1[2]:>7.2f}   {pts[j,0]:>9.3f}  {pts[j,1]:>6.2f}  {pts[j,2]:>7.2f}   {wp[0]:>8.3f}  {wp[1]:>6.2f}  {wp[2]:>6.1f}")
        results[mA] = (pts, r1, wp)

    if args.write_json:
        # Write the CHOSEN working point (chosen_wp = FIXED_WP for mA2/3, R=1 crossing for mA1),
        # NOT the bare R=1 crossing -> --write-json can no longer revert mA2/mA3.
        write_json(JSON_PATH, {mA: results[mA][2] for mA in LOWM})
    else:
        print("[write-json] skipped (pass --write-json to overwrite mA1,2,3 MVAcut in the JSON)")

    # combined 3-pad overview
    c = ROOT.TCanvas("c", "", 1500, 520); c.Divide(3, 1)
    for i, mA in enumerate(LOWM):
        c.cd(i + 1); draw_one(mA, *results[mA], big=False)
    out = f"{OUTDIR}/lowmass_score_R_significance.pdf"
    c.SaveAs(out); c.SaveAs(out.replace(".pdf", ".png")); print("saved", out)

    # one stand-alone plot per mA with larger axis/label/text
    for mA in LOWM:
        cc = ROOT.TCanvas(f"c{mA}", "", 820, 740); cc.cd(); draw_one(mA, *results[mA], big=True)
        o = f"{OUTDIR}/lowmass_score_R_significance_mA{mA}.pdf"
        cc.SaveAs(o); cc.SaveAs(o.replace(".pdf", ".png")); print("saved", o)

    # stand-alone mA30 panel (HIGH-MASS model). mA30 has R~1.2-1.3 across the whole scan (sculpting
    # risk flagged by the AN R-table); its working point (red star) is the MVAcut from JSON_PATH.
    mA = 30
    pts30 = scan_high(mA); r1_30 = r_equals_one(pts30); wp30 = chosen_wp(mA, pts30)
    print(f"{mA:>3}  {r1_30[0]:>9.3f}  {r1_30[2]:>7.2f}   {'-':>9}  {'-':>6}  {'-':>7}   "
          f"{wp30[0]:>8.3f}  {wp30[1]:>6.2f}  {wp30[2]:>6.1f}   (high-mass model)")
    cc = ROOT.TCanvas("c30", "", 820, 740); cc.cd(); draw_one(mA, pts30, r1_30, wp30, big=True)
    o = f"{OUTDIR}/highmass_score_R_significance_mA{mA}.pdf"
    cc.SaveAs(o); cc.SaveAs(o.replace(".pdf", ".png")); print("saved", o)

if __name__ == "__main__":
    main()
