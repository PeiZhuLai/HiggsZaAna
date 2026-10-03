"""Add the normalized signal reweight (as in apply_bdt_sig.py / signal_eff_sumw.py) to the three
migration plotters. usage: patch_migration_sigrw.py <script.py> ..."""
import re, sys
HELPER = '''

# [PZ 2026-10-03] Signal reweight: the nominal sideband reweight with the true-mass param (ALP_m - m_a)/H_m,
# normalized per (m_a, era, channel) over the tree that is read, as in apply_bdt_sig.py and
# signal_eff_sumw.py. HZA_SIGNAL_REWEIGHT=0 restores the unreweighted signal.
_SIGNAL_RW_JSON = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/reweights/sideband_run3_iterative.json"
_SIGNAL_RW_OBJ = None
_SIGNAL_RW_CACHE = {}

def _signal_rw_weights(w, fp, tname, true_ma, t, wname):
    """Return w multiplied by the per-event normalized signal reweight (w must cover the whole tree)."""
    import os as _os
    if _os.environ.get("HZA_SIGNAL_REWEIGHT", "1") == "0" or true_ma is None or not wname:
        return w
    global _SIGNAL_RW_OBJ
    key = (str(fp), str(tname))
    if key not in _SIGNAL_RW_CACHE:
        if _SIGNAL_RW_OBJ is None:
            import sys as _sys
            _sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/scripts")
            from sideband_reweight import SidebandReweighter
            _SIGNAL_RW_OBJ = SidebandReweighter.from_json(_SIGNAL_RW_JSON)
        fr = t.arrays(library="pd")
        for logical, cands in (("pho1ECALIso", ("pho1PIso_noCorr",)), ("pho2ECALIso", ("pho2PIso_noCorr",)),
                               ("H_m", ("H_mass",)), ("ALP_m", ("ALP_mass",))):
            if logical not in fr.columns:
                for c in cands:
                    if c in fr.columns:
                        fr[logical] = fr[c]; break
        fr["param"] = (fr["ALP_m"].to_numpy(dtype=float) - float(true_ma)) / fr["H_m"].to_numpy(dtype=float)
        r = np.asarray(_SIGNAL_RW_OBJ.weights_for_dataframe(fr), dtype=float)
        ww = fr[wname].to_numpy(dtype=float)
        out = np.ones(len(fr), dtype=float)
        for col in ("n_electrons", "n_muons"):
            sel = fr[col].to_numpy() == 2
            den = float(np.sum(ww[sel] * r[sel]))
            out[sel] = r[sel] * (float(np.sum(ww[sel])) / den if den > 0 else 1.0)
        _SIGNAL_RW_CACHE[key] = out
    rw = _SIGNAL_RW_CACHE[key]
    if len(w) != len(rw):
        raise RuntimeError("signal reweight needs the whole tree in one chunk: %d vs %d (%s)" % (len(w), len(rw), fp))
    return w * rw
'''
for p in sys.argv[1:]:
    s = open(p).read()
    anchor = re.search(r"def _true_ma_from_dir\(ma_tag: str\) -> Optional\[int\]:\n(    .*\n)+", s)
    assert anchor, p
    s = s[:anchor.end()] + HELPER + s[anchor.end():]
    n = 0
    lines = s.split("\n")
    for i, l in enumerate(lines):
        if "ak.values_astype(arrs[wname], np.float64) if wname else" in l:
            # which true-mass variable is in scope: true_ma, else parse ma_tag
            ctx = "\n".join(lines[max(0, i - 80):i])
            tm = "true_ma" if re.search(r"\btrue_ma = _true_ma_from_dir", ctx) else "_true_ma_from_dir(ma_tag)"
            lines[i] = l.replace("ak.values_astype(arrs[wname], np.float64) if wname else",
                                 "_signal_rw_weights(ak.values_astype(arrs[wname], np.float64), fp, tname, %s, t, wname) if wname else" % tm)
            n += 1
    open(p, "w").write("\n".join(lines))
    print(p, "weight reads patched:", n)
