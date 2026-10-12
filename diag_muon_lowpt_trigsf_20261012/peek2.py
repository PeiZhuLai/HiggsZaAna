import uproot,sys,re
t=uproot.open(sys.argv[1])["inclusive"]
print([b for b in t.keys() if re.match(r"(muon|electron)_\d", b) or b.startswith(("n_","pass_","z_"))])
