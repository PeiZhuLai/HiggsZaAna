"""V9 prep: compare R at the adopted WP from (a) stored MVA_Score_mA_M{m} on the scored
combined DY background (AN Table 17 definition) vs (b) re-scoring All_Bkg with the high-mass pkl
(existing scan_high path used for the mA30 panel)."""
import sys, json, numpy as np
sys.argv = ["x"]
sys.path.insert(0, "/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/Plot/scripts")
import scan_score_R_significance as S
cuts = {int(e["mA"]): float(e["MVAcut"]) for e in json.load(open(S.JSON_PATH))["results"]}
bkA = S.load(f"{S.ROOT_DIR}/All_Bkg/run3.root", S.BASE_VARS)
print("All_Bkg entries", len(bkA["factor"]), "sumw", bkA["factor"].sum(), "neg-w frac", (bkA["factor"]<0).mean())
for m in [4,5,6,7,8,9,10,15,20,25,30]:
    bk = S.load_scored_bkg(m)
    Rst = S._R_direct(bk, cuts[m])
    s = S.score_high(bkA, m)
    w = np.clip(bkA["factor"], 0, None); H = bkA["H_mass"]
    pk = (H>120)&(H<130); frac = w[pk].sum()/w.sum(); p = s>cuts[m]
    Rpk = w[p&pk].sum()/w[p].sum()/frac
    print(f"mA{m:2d} cut {cuts[m]:.3f}  stored: N={len(bk['factor'])} sumw={bk['factor'].sum():.1f} R={Rst:.3f} npass={(bk['s']>cuts[m]).sum()} npk={((bk['s']>cuts[m])&(bk['H_mass']>120)&(bk['H_mass']<130)).sum()} | pkl: R={Rpk:.3f} npass={p.sum()} npk={(p&pk).sum()}")
