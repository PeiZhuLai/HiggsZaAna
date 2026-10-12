import uproot,sys
f=uproot.open(sys.argv[1]); print(list(f.keys())[:20])
t=f[[k for k in f.keys() if "inclusive" in k.lower() or "SR" in k][0]] if any("inclusive" in k.lower() for k in f.keys()) else f[f.keys()[0]]
print(t.name, t.num_entries)
print([b for b in t.keys() if any(s in b.lower() for s in ("mu","trig","lepton","weight","sf_"))][:120])
