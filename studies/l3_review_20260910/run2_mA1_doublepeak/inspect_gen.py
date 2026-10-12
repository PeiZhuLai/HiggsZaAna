"""Print the mothers of generator photons in one NanoAOD file (find the ALP pdgId)."""
import sys, collections
import numpy as np, uproot, awkward as ak
f = uproot.open(sys.argv[1], timeout=300)
a = f["Events"].arrays(["GenPart_pdgId", "GenPart_genPartIdxMother", "GenPart_status"], entry_stop=2000)
pid, mo = a.GenPart_pdgId, a.GenPart_genPartIdxMother
c = collections.Counter()
for p, m in zip(pid[:500], mo[:500]):
    p = np.asarray(p); m = np.asarray(m)
    for i in np.where(p == 22)[0]:
        if m[i] >= 0:
            c[int(p[m[i]])] += 1
print(c.most_common(8))
print([k for k in f["Events"].keys() if k.startswith("Photon_")][:60])
