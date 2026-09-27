# Data sideband counts (95-115 & 135-180 GeV; SR 115-135 NOT looked at) in the MVA-cut data trees,
# to cross-check the MC ratios between neighbouring anchor / interpolated points.
import uproot, numpy as np
B="/eos/home-p/pelai/HZa/root_MVAcut/data"
res={}
for m in (10,11,12,15,16,19,20,21,22,24,25,26,29,30):
    t=uproot.open(f"{B}/mA_M{m}/run3.root")["DiphotonTree/Data_13p6TeV"]
    x=t["CMS_hza_mass"].array(library="np")
    sb=((x>95)&(x<115))|((x>135)&(x<180))
    res[m]=int(sb.sum())
    print(m, "sideband N =", res[m])
