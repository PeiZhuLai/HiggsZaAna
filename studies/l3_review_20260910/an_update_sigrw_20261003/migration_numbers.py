"""Correct-m_a probability (diagonal of the migration, Run 3, test tree, as in plot_mAmigratedBar.py
_accumulate_pred_pass_sumw_by_year) per hypothesis, without and with the normalized signal reweight.
Also the fraction of each generated sample passing its own working point (Sec 6.5 '60% to 97%' candidate).
usage: migration_numbers.py   Output: stdout"""
import json, sys, numpy as np, uproot
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts")
import plot_mAmigratedBar as B
B_DIR = "/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"
MAS = [1,2,3,4,5,6,7,8,9,10,15,20,25,30]; ERAS = ["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
cuts = {e["mA"]: e["MVAcut"] for e in json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
S = {0: {}, 1: {}}; own = {0: {}, 1: {}}
for t_ma in MAS:
    for era in ERAS:
        fp = "%s/mA_M%d/%s.root" % (B_DIR, t_ma, era)
        t = uproot.open(fp)["test"]
        a = t.arrays(["weight"] + ["MVA_Score_mA_M%d" % p for p in MAS], library="np")
        w0 = a["weight"]
        w1 = np.asarray(B._signal_rw_weights(w0, fp, "test", t_ma, t, "weight"))
        for k, w in ((0, w0), (1, w1)):
            for p in MAS:
                S[k].setdefault(p, {}).setdefault(t_ma, 0.0)
                S[k][p][t_ma] += w[a["MVA_Score_mA_M%d" % p] >= cuts[p]].sum()
            o = own[k].setdefault(t_ma, [0.0, 0.0]); o[0] += w[a["MVA_Score_mA_M%d" % t_ma] >= cuts[t_ma]].sum(); o[1] += w.sum()
for k in (0, 1):
    diag = {p: S[k][p][p] / sum(S[k][p].values()) for p in MAS}
    eff = {m: own[k][m][0] / own[k][m][1] for m in MAS}
    print("rw=%d correct-mA prob:" % k, " ".join("%d:%.3f" % (p, diag[p]) for p in MAS), " range %.3f-%.3f" % (min(diag.values()), max(diag.values())))
    print("rw=%d own-WP eff (ee+mumu, test, all m_llgg):" % k, " ".join("%d:%.3f" % (m, eff[m]) for m in MAS))
