#!/usr/bin/env python3
"""Numbers for the AN/PAS update to low_relE2_w123 (read-only). Output: an_numbers/w123_numbers.txt
closure (pseudodata_closure/closure_results.json, 30 points), bias (Bias_nominal/bias_outputs_highstat mA_1-3),
effSigma (Signal/outdir_*/signalFit/Plots/effSigma_<m>_<era>_<ch>.json), m_a = 1 signal fits (signalFit log of the
deployment chain: in-window fraction, number of Gaussians, final chi2/ndf; order = ele eras then mu eras)."""
import glob, json, re
F = "/afs/cern.ch/work/p/pelai/HZa/flashgg_run3/CMSSW_14_1_0_pre4/src/flashggFinalFit"
D = "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/studies/l3_review_20260910/deploy_w123_20261017"
ERAS = ["2022preEE", "2022postEE", "2023preBPix", "2023postBPix", "2024"]
out = []
c = json.load(open(F + "/pseudodata_closure/closure_results.json"))
in1 = in2 = 0; up = []; dn = []
for m in sorted(c, key=int):
    l = c[m]["limit"]; o = l["obs"]
    i1 = l["exp_m1"] <= o <= l["exp_p1"]; i2 = l["exp_m2"] <= o <= l["exp_p2"]
    in1 += i1; in2 += i2
    if not i1: (up if o > l["exp_p1"] else dn).append(int(m))
sig = sorted(((float(v["signif"]), int(m), v["rhat"]["r"]) for m, v in c.items()), reverse=True)[:4]
out.append("closure: %d/30 within +-1sigma, %d/30 within +-2sigma; below -1sigma %s; above +1sigma %s" % (in1, in2, dn, up))
out.append("closure: largest local significance " + ", ".join("mA%d %.2f (rhat %+.3f)" % (m, s, r) for s, m, r in sig))
for m in ("1", "2", "3"):
    v = c[m]; l = v["limit"]
    out.append("closure mA%s: obs %.3f exp %.3f [-1s %.3f, +1s %.3f, -2s %.3f, +2s %.3f] Z %.2f rhat %+.3f -%.3f/+%.3f" % (
        m, l["obs"], l["exp_med"], l["exp_m1"], l["exp_p1"], l["exp_m2"], l["exp_p2"], v["signif"], v["rhat"]["r"], abs(v["rhat"]["rLoErr"]), v["rhat"]["rHiErr"]))
others = [abs(v["rhat"]["r"]) for m, v in c.items() if int(m) >= 3 and int(m) != 5]
out.append("closure: max |rhat| for mA>=3 except 5: %.3f" % max(others))
for m in (1, 2, 3):
    fr = json.load(open("%s/Combine/Checks/Bias_nominal/bias_outputs_highstat/mA_%d/merged/BiasJson/%d_gaussfit.json" % (F, m, m)))["fit_results"]
    out.append("bias mA%d: " % m + " ".join("%s %+.3f" % (k, v["mean"]) for k, v in sorted(fr.items(), key=lambda x: -abs(x[1]["mean"]))))
for m in (1, 2, 3):
    vals = []
    for ch in ("ele", "mu"):
        for e in ERAS:
            vals.append(json.load(open("%s/Signal/outdir_%s/signalFit/Plots/effSigma_%d_%s_%s.json" % (F, ch, m, e, ch)))[e] if False else
                        list(json.load(open("%s/Signal/outdir_%s/signalFit/Plots/effSigma_%d_%s_%s.json" % (F, ch, m, e, ch))).values())[0])
    comb = {ch: json.load(open("%s/Signal/outdir_%s/signalFit/Plots/effSigma_%d_run3_%s.json" % (F, ch, m, ch)))["combined"] for ch in ("ele", "mu")}
    out.append("effSigma mA%d: per era/channel %.2f-%.2f GeV; combined ele %.2f mu %.2f" % (m, min(vals), max(vals), comb["ele"], comb["mu"]))
log = open(D + "/logs/sig_3_runjob_sig_signalFit_m1.log").read().split(" --> Running fit for (proc,cat)")[1:]
keys = [(ch, e) for ch in ("ele", "mu") for e in ERAS]
assert len(log) == 10, len(log)
kept, chi = [], []
for k, b in zip(keys, log):
    kp = float(re.search(r"kept ([0-9.]+)%", b).group(1)); ng = int(re.search(r"convolution of nGaussians \((\d+)\)", b).group(1))
    c2 = float(re.findall(r"\(Total\) Running fit.*?chi2/n\(dof\) = ([0-9.]+)", b, re.S)[-1])
    kept.append(kp); chi.append(c2)
    out.append("mA1 sigfit %-4s %-13s kept %.1f%% nGauss %d chi2/ndf %.3f" % (k[0], k[1], kp, ng, c2))
out.append("mA1 sigfit: kept %.1f-%.1f%%, chi2/ndf %.2f-%.2f" % (min(kept), max(kept), min(chi), max(chi)))
meta = json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HZaMVA/using/model_Za_BDT_lowmass_run3.meta.json"))
out.append("low model: AUC test %.4f, KS sig %.3f bkg %.3f" % (meta["auc_test"], meta["ks_sig"], meta["ks_bkg"]))
open(D + "/an_numbers/w123_numbers.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
