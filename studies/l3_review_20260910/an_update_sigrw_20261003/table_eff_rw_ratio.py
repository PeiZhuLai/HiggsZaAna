"""tab:boundary-like efficiency (ee+mumu, 115<m_llgg<135, test tree, WP floored to the 0.005 score grid of the
optimization histograms) without and with the normalized signal reweight; indicates how the table would move."""
import json, math, sys, numpy as np, uproot
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts")
import plot_mAmigratedBar as B
D="/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_fsrfix"; ERAS=["2022preEE","2022postEE","2023preBPix","2023postBPix","2024"]
cuts={e["mA"]:e["MVAcut"] for e in json.load(open("/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/output/MVAcut_points_run3.json"))["results"]}
for m in [1,2,3,4,5,6,7,8,9,10,15,20,25,30]:
    c=math.floor(cuts[m]/0.005+1e-6)*0.005; n0=d0=n1=d1=0
    for e in ERAS:
        fp=f"{D}/mA_M{m}/{e}.root"; t=uproot.open(fp)["test"]
        a=t.arrays(["weight","H_mass",f"MVA_Score_mA_M{m}"],library="np")
        w0=a["weight"]; w1=np.asarray(B._signal_rw_weights(w0,fp,"test",m,t,"weight"))
        sr=(a["H_mass"]>115)&(a["H_mass"]<135); ps=sr&(a[f"MVA_Score_mA_M{m}"]>=c)
        n0+=w0[ps].sum(); d0+=w0[sr].sum(); n1+=w1[ps].sum(); d1+=w1[sr].sum()
    print("mA%-2d eff_noRW %.3f eff_RW %.3f ratio %.3f" % (m, n0/d0, n1/d1, (n1/d1)/(n0/d0)))
