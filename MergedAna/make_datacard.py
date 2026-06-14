#!/usr/bin/env python3
"""Write one parametric combine datacard per m_a for the merged low-mA category.

Single channel `ch1`, observable H_mass.
  signal : double-CB from sig_ws_<tag>.root, normalisation = constant
           sig_model_<tag>_norm (yield at xs=0.1 pb) -> limit is in units of
           that reference signal strength r.
  bkg    : RooMultiPdf from bkg_ws.root, index nuisance pdfindex_ch1 (discrete).
"""
import argparse
import json
import os


CARD = """# merged low-mA  m_a={ma} GeV  ({src})
imax 1
jmax 1
kmax *
----------------------------------------------------------------------
shapes sig      ch1 {sig_ws}  wsig_{tag}:sig_model_{tag}
shapes bkg      ch1 {bkg_ws}  wbkg:roomultipdf
shapes data_obs ch1 {bkg_ws}  wbkg:data_obs
----------------------------------------------------------------------
bin          ch1
observation  -1
----------------------------------------------------------------------
bin              ch1     ch1
process          sig     bkg
process          0       1
rate             1       1
----------------------------------------------------------------------
lumi_13p6TeV  lnN  1.018   -
CMS_eff_merged lnN 1.05    -
pdfindex_ch1  discrete
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(os.path.dirname(__file__), "output"))
    args = ap.parse_args()

    ws_dir = os.path.join(args.outdir, "workspaces")
    card_dir = os.path.join(args.outdir, "datacards")
    os.makedirs(card_dir, exist_ok=True)

    yields = json.load(open(os.path.join(args.outdir, "signal_yields.json")))
    src_file = os.path.join(args.outdir, "bkg_source.txt")
    src = open(src_file).read().strip() if os.path.exists(src_file) else "unknown"
    bkg_ws = os.path.abspath(os.path.join(args.outdir, "bkg_ws.root"))

    n = 0
    for r in yields:
        tag = r["tag"]
        sig_ws = os.path.abspath(os.path.join(ws_dir, f"sig_ws_{tag}.root"))
        if not os.path.exists(sig_ws):
            print(f"[card] {tag}: missing {sig_ws} -> skip")
            continue
        txt = CARD.format(ma=r["ma"], src=src, tag=tag, sig_ws=sig_ws, bkg_ws=bkg_ws)
        path = os.path.join(card_dir, f"datacard_{tag}.txt")
        with open(path, "w") as f:
            f.write(txt)
        n += 1
        print(f"[card] wrote {path}")
    print(f"[card] {n} datacards in {card_dir}  (bkg source: {src})")


if __name__ == "__main__":
    main()
