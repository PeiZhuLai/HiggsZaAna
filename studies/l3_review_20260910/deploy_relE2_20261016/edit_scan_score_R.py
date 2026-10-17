#!/usr/bin/env python3
"""One-off edit of Plot/scripts/scan_score_R_significance.py for the low_swapPt_relE2 low-mass BDT:
pin the m_a 1-3 working points to the relE2 cuts (so --write-json cannot revert them) and build the
low-mass feature matrix by name (hza_features.feature_matrix), which knows pho2_relEerr.
Backup: scan_score_R_significance.py.bak_preRelE2_20261016"""
P = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts/scan_score_R_significance.py"
s = open(P).read()
pairs = [
    ('FEATS    = ["ALP_calculatedPhotonIso", "var_dR_g1g2", "pho1R9", "pho1Pt_oHm"]   # + param\n',
     'FEATS    = [f for f in __import__("json").load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts/'
     'model_Za_BDT_lowmass_run3.meta.json"))["features"] if f != "param"]   # low-mass model inputs (+ param)\n'),
    ('FIXED_WP = {2: 0.975, 3: 0.988}   # mA3 0.992 -> 0.988 (2026-10-01): after the 2024 DY+jets overlap-veto\n',
     '# Low-mass BDT low_swapPt_relE2: the m_a 1-3 cuts come from the x0.5 rule of the background-only closure scan\n'
     '# (closure + data GOF pass), not from the R=1 crossing; pinned here so that --write-json keeps them.\n'
     'FIXED_WP = {1: 0.98036, 2: 0.98661, 3: 0.99194}\n'
     '# Previous low-mass model (pho1Pt_oHm input): {2: 0.975, 3: 0.988}, mA1 at the R=1 crossing 0.963.\n'
     '# mA3 0.992 -> 0.988: after the 2024 DY+jets overlap-veto\n'),
    ('from hza_features import BASE_VARS, FEATURES   # FEATURES == BASE_VARS + [\'pho_pt_asym\',\'param\']\n\nmdl      =',
     'from hza_features import BASE_VARS, FEATURES, feature_matrix, relerr_extra   # FEATURES == BASE_VARS + [\'pho_pt_asym\',\'param\']\n\nmdl      ='),
    ('''def load(path, feats=FEATS):
    a = uproot.open(path)["inclusive"].arrays(list(feats) + ["H_mass", "ALP_mass", "factor"], library="np")
    m = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"])
    return {k: a[k][m] for k in a}

def score(a, mh):
    X = np.column_stack([a[c] for c in FEATS] + [(a["ALP_mass"] - mh) / a["H_mass"]])
    return mdl.predict_proba(X)[:, 1]''',
     '''def load(path, feats=None):
    """Inclusive tree in 95-180 GeV; default columns: the low-mass model inputs (derived ones expanded)."""
    if feats is None:
        feats = [f for f in FEATS if f not in ("pho1_relEerr", "pho2_relEerr")] + relerr_extra(FEATS)
    a = uproot.open(path)["inclusive"].arrays(list(feats) + ["H_mass", "ALP_mass", "factor"], library="np")
    m = (a["H_mass"] > 95) & (a["H_mass"] < 180) & np.isfinite(a["factor"])
    return {k: a[k][m] for k in a}

def score(a, mh):
    b = dict(a, ALP_m=a["ALP_mass"], H_m=a["H_mass"])
    return mdl.predict_proba(feature_matrix(b, mh, FEATS + ["param"]))[:, 1]'''),
    ('''    from hza_features import BASE_VARS, FEATURES   # FEATURES == BASE_VARS + ['pho_pt_asym','param']

    mdl = pickle.load(open(model_path, "rb"))
    feats = json.load(open(meta_path))["features"]          # model feature order (incl 'param')
    idx = [FEATURES.index(f) for f in feats]
    a = uproot.open(f"{ROOT_DIR}/All_Bkg/run3.root")["inclusive"].arrays(
        BASE_VARS + ["H_mass", "ALP_mass", "factor"], library="np")''',
     '''    from hza_features import BASE_VARS, feature_matrix, relerr_extra

    mdl = pickle.load(open(model_path, "rb"))
    feats = json.load(open(meta_path))["features"]          # model feature order (incl 'param')
    a = uproot.open(f"{ROOT_DIR}/All_Bkg/run3.root")["inclusive"].arrays(
        BASE_VARS + relerr_extra(feats) + ["H_mass", "ALP_mass", "factor"], library="np")'''),
    ('''    asym = (a["pho1Pt_oHm"] - a["pho2Pt_oHm"]) / (a["pho1Pt_oHm"] + a["pho2Pt_oHm"] + 1e-6)
    base_cols = [a[b] for b in BASE_VARS]
    R = {}
    for mA, cut in cuts.items():
        param = (a["ALP_mass"] - mA) / a["H_mass"]
        X = np.column_stack(base_cols + [asym, param])[:, idx]  # FEATURES order -> model order''',
     '''    a["ALP_m"], a["H_m"] = a["ALP_mass"], a["H_mass"]
    R = {}
    for mA, cut in cuts.items():
        X = feature_matrix(a, float(mA), feats)                 # model feature order, by name'''),
]
for old, new in pairs:
    assert s.count(old) == 1, old[:70]
    s = s.replace(old, new)
open(P, "w").write(s)
print("edited", P)
