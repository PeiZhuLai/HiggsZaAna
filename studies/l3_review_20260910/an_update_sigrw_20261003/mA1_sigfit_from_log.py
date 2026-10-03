"""m_a=1 signal fits (Sec 8.1): in-window yield fraction, final chi2/ndf and number of Gaussians, per
(channel, era), parsed from a signalFit driver log (shellScripts/sig_sys/3_runjob_sig_signalFit.sh output).
usage: mA1_sigfit_from_log.py <log>"""
import re, sys
cur = None; out = {}
for line in open(sys.argv[1]):
    m = re.search(r"root_MVAcut/sig\S*/mA_M(\d+)/ws_Tree2WS/ws_(ele|mu)_(\S+)\.root", line)
    if m:
        cur = (m.group(2), m.group(3)) if m.group(1) == "1" else None
        if cur: out[cur] = {}
        continue
    if cur is None: continue
    m = re.search(r"kept ([0-9.]+)%", line)
    if m: out[cur]["kept"] = float(m.group(1))
    m = re.search(r"nGaussians \((\d+)\)", line)
    if m: out[cur]["n"] = int(m.group(1))
    m = re.search(r"chi2/n\(dof\) = ([0-9.]+)", line)
    if m: out[cur]["chi2ndf"] = float(m.group(1))   # last one = post-fit
for k, v in out.items(): print(k, v)
ks = [v["kept"] for v in out.values()]; cs = [v["chi2ndf"] for v in out.values()]
print("kept %.1f-%.1f%%  chi2/ndf %.2f-%.2f  n_fits %d" % (min(ks), max(ks), min(cs), max(cs), len(out)))
