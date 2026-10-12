"""Add the Run 2 UL NanoAODv9 H->Za->2l2g datasets (from DAS) to the mA_M<mass> entries of the sample list."""
import json, re, sys
das_txt, sample_json = sys.argv[1], sys.argv[2]
CAMPAIGN_ERA = {"RunIISummer20UL16NanoAODAPVv9": "2016preVFP", "RunIISummer20UL16NanoAODv9": "2016postVFP",
                "RunIISummer20UL17NanoAODv9": "2017", "RunIISummer20UL18NanoAODv9": "2018"}
RUN2 = ["2016preVFP", "2016postVFP", "2017", "2018"]
found = {}
for line in open(das_txt):
    ds = line.strip()
    m = re.match(r"/HZaTo2l2g_M-?([0-9p]+)_[^/]*/(RunIISummer20UL\d+NanoAOD(?:APV)?v9)-[^/]*/NANOAODSIM$", ds)
    if not m:
        continue
    key, era = "mA_M" + m.group(1), CAMPAIGN_ERA[m.group(2)]
    assert (key, era) not in found, ("duplicate", key, era, ds)
    found[(key, era)] = ds
s = open(sample_json).read()
d = json.loads(s)
assert json.dumps(d, indent=4) == s, "formatting would change"
added = 0
for key in sorted({k for k, _ in found}):
    if key not in d:
        sys.exit("no entry %s in sample list" % key)
    files = d[key]["files"]
    eras = [e for e in RUN2 if (key, e) in found]
    assert len(eras) == 4, (key, eras)
    for e in eras:
        if e in files and files[e] != [found[(key, e)]]:
            sys.exit("conflict %s %s: %s" % (key, e, files[e]))
    new = {e: [found[(key, e)]] for e in eras}
    new.update({e: v for e, v in files.items() if e not in new})
    added += sum(e not in files for e in eras)
    d[key]["files"] = new
open(sample_json, "w").write(json.dumps(d, indent=4))
print("entries:", len({k for k, _ in found}), "era paths added:", added)
