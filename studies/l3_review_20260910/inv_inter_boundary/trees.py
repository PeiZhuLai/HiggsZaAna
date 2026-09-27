import uproot, sys
P="/eos/home-p/pelai/HZa/root_P2Root"
for tag in ("run3_bdt_scored_fsrfix","run3_bdt_scored_nominal"):
    for s,y in (("DYGto2LG_10to100","2024"),("DYJetsTo2Mu","2024"),("DYJetsToLL","2023preBPix")):
        f=uproot.open(f"{P}/{tag}/{s}/{y}.root")
        ks=[(k.split(';')[0],f[k].classname) for k in f.keys()]
        tt=[k for k,c in ks if 'TTree' in c]
        print(tag,s,y,"trees:",[(k,f[k].num_entries) for k in tt][:8])
        t=f[tt[0]] if 'inclusive' not in [k for k,_ in ks] else f['inclusive']
        br=t.keys()
        print("   has weight:", 'weight' in br, "factor:", 'factor' in br, "H_m:", 'H_m' in br, "n MVA branches:", len([b for b in br if b.startswith('MVA_Score_mA_M')]))
