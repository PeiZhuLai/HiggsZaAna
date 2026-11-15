#!/usr/bin/env python3
"""Re-align the columns of a Combine datacard (whitespace only).

The nuisance renaming replaced names in place without re-padding, so the type and value columns of the
systematics block no longer line up with the bin/process/rate header. This script rewrites every
table row (bin, process, rate, systematics, param) with common column widths and leaves all other lines
unchanged. It refuses to write if the whitespace-split tokens of any line change.
Usage: align_datacard.py <card.txt> [<card.txt> ...]   (in place)
"""
import sys

TYPES = {"lnN", "lnU", "gmN", "shape", "shape?", "shapeN", "param", "flatParam", "discrete", "rateParam", "extArg"}


def align(lines):
    rows = []  # (index, kind, tokens)
    for i, l in enumerate(lines):
        t = l.split()
        if not t:
            continue
        if t[0] in ("bin", "process", "rate") and len(t) > 2:
            rows.append((i, "head", t))
        elif len(t) >= 2 and t[1] in TYPES:
            rows.append((i, "syst", t))
    name_w = max(len(t[0]) for _, k, t in rows if k == "syst") + 2
    type_w = max(len(t[1]) for _, k, t in rows if k == "syst") + 2
    col_w = max(len(x) for _, k, t in rows if k in ("head", "syst") and (k == "head" or t[1] in ("lnN", "lnU"))
                for x in (t[1:] if k == "head" else t[2:])) + 2
    out = list(lines)
    for i, k, t in rows:
        if k == "head":
            s = t[0].ljust(name_w + type_w) + "".join(x.ljust(col_w) for x in t[1:])
        elif t[1] in ("lnN", "lnU"):
            s = t[0].ljust(name_w) + t[1].ljust(type_w) + "".join(x.ljust(col_w) for x in t[2:])
        else:  # param, discrete, rateParam, ...: name and type aligned, values separated by two spaces
            s = t[0].ljust(name_w) + t[1].ljust(type_w) + "  ".join(t[2:])
        out[i] = s.rstrip() + "\n"
    return out


for card in sys.argv[1:]:
    lines = open(card).readlines()
    new = align(lines)
    if [l.split() for l in lines] != [l.split() for l in new]:
        sys.exit("token mismatch in %s, not written" % card)
    open(card, "w").writelines(new)
    print("aligned", card)
