import sys, ROOT
ROOT.gROOT.SetBatch(True); ROOT.TH1.AddDirectory(False)
src, dst = sys.argv[1], sys.argv[2]
fi = ROOT.TFile.Open(src); fo = ROOT.TFile.Open(dst, "RECREATE")
kept = dropped = 0
for dname in ("raw_plots", "sys_dir"):
    di = fi.Get(dname); do = fo.mkdir(dname); do.cd()
    for k in di.GetListOfKeys():
        n = k.GetName()
        if dname == "sys_dir" and "_DYJetsToLL_" not in n and "_DYGto2LG_" not in n:
            dropped += 1; continue
        k.ReadObj().Write(n); kept += 1
fo.Close(); print("kept", kept, "dropped", dropped)
