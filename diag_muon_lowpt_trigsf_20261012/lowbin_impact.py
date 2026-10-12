"""HZa signal mumu: effect of the sub-threshold (lowest-pT) muon trigger bins on the event trigger SF.
Recomputes P_data/P_MC (HLT_systematics.py formula) from the two Z muons with the current 2024 JSONs."""
import sys, numpy as np, uproot, correctionlib
D="/afs/cern.ch/work/p/pelai/HZa/HiggsZaAna/HiggsDNA/higgs_dna/systematics/data/2024_UL/"
cs={leg: correctionlib.CorrectionSet.from_file(D+f"hza_mutrig{leg}_2024_efficiencies.json") for leg in (8,17,24)}
LOW={8:8.0, 17:16.75, 24:23.0}          # upper edge of the lowest (sub-threshold) pT bin
def eff(leg, pt, eta, kind, mode):
    e = cs[leg][kind].evaluate(pt, np.abs(eta))
    low = pt < LOW[leg]
    if mode == "sf1" and kind == "effdata":     # lowest bin: data eff := MC eff (SF=1)
        e = np.where(low, cs[leg]["effmc"].evaluate(pt, np.abs(eta)), e)
    if mode == "gate":                          # offline gate: below threshold cannot fire the leg
        e = np.where(low, 0.0, e)
    return e
def prob(p1,e1,p2,e2,kind,mode):
    s1,s2=eff(24,p1,e1,kind,mode),eff(24,p2,e2,kind,mode)
    h1,h2=eff(17,p1,e1,kind,mode),eff(17,p2,e2,kind,mode)
    l1,l2=eff(8,p1,e1,kind,mode),eff(8,p2,e2,kind,mode)
    ps=1-(1-s1)*(1-s2); pd=h1*l2+l1*h2-h1*h2
    return ps+pd-ps*pd
br=["Z_lead_lepton_pt","Z_lead_lepton_eta","Z_sublead_lepton_pt","Z_sublead_lepton_eta","z_mumu","weight_central","weight_hlt_sf_central","n_muons"]
print(f"{'mA':>4} {'N':>6} {'sub<16.75':>9} {'sub<23':>7} {'lead<23':>7} | {'<SF>nom':>8} {'dY(sf1)%':>9} {'dY(gate)%':>9} | {'stored/mine':>11} {'n_mu>2':>6}")
for m in sys.argv[1:]:
    t=uproot.open(f"/eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_fsrfix/mA_M{m}/2024.root")["inclusive"]
    a=t.arrays(br,library="np"); k=a["z_mumu"]==1
    p1,e1,p2,e2,w=(a[x][k] for x in br[:4]+["weight_central"])
    w=w/a["weight_hlt_sf_central"][k]           # weight without trigger SF
    r={}
    for mode in ("nom","sf1","gate"):
        r[mode]=prob(p1,e1,p2,e2,"effdata",mode)/prob(p1,e1,p2,e2,"effmc",mode)
    Y={mo:(w*r[mo]).sum() for mo in r}
    st=a["weight_hlt_sf_central"][k]/r["nom"]
    print(f"{m:>4} {k.sum():>6} {np.average(p2<16.75,weights=w):9.3f} {np.average(p2<23,weights=w):7.3f} {np.average(p1<23,weights=w):7.3f} | {np.average(r['nom'],weights=w):8.4f} {100*(Y['sf1']/Y['nom']-1):9.3f} {100*(Y['gate']/Y['nom']-1):9.3f} | {np.median(st):5.3f}/{np.percentile(np.abs(st-1),99):.3f} {np.average(a['n_muons'][k]>2,weights=w):6.3f}")
