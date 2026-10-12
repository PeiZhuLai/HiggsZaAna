"""Job config for one Run 2 signal file: the Run 3 job config of mA_M1_2024 job_1 with the physics
content (branches, branch_map, systematics, tag_sequence, variables_of_interest) taken from
metadata/za_signal_run2.json and the sample switched to mA_M1_<year>."""
import ast, json, sys
job3, run2cfg, year, out = sys.argv[1:5]
c = json.load(open(job3))
r = json.load(open(run2cfg))
for k in ("branches", "branch_map", "systematics", "tag_sequence", "variables_of_interest"):
    c[k] = r[k]
tree = ast.parse(open("/afs/cern.ch/work/p/pelai/HZgamma/higgsdna-hzg-run3/higgs_dna/metaconditions/corrections/constants.py").read())
lumi = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", None) == "LUMI")
c["sample"].update({"name": "mA_M1_%s" % year, "year": year, "lumi": lumi[year]})
c["name"] = "mA_M1_%s_1" % year
json.dump(c, open(out, "w"), indent=1)
print("year", year, "lumi", lumi[year])
