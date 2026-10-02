import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.gErrorIgnoreLevel = ROOT.kFatal
p = sys.argv[1]
f = ROOT.TFile.Open(p)
if not f or f.IsZombie(): print("BAD", p, "zombie"); sys.exit(1)
trees = []
def walk(d, pre=""):
    for k in d.GetListOfKeys():
        o = k.ReadObj()
        if o.InheritsFrom("TDirectory"): walk(o, pre + k.GetName() + "/")
        elif o.InheritsFrom("TTree"): trees.append((pre + k.GetName(), o))
walk(f)
if not trees: print("BAD", p, "no TTree"); sys.exit(1)
for name, t in trees:
    n = t.GetEntries()
    for i in range(n):
        if t.GetEntry(i) <= 0: print("BAD", p, "%s read failed at %d/%d" % (name, i, n)); sys.exit(1)
print("OK", p); sys.exit(0)
