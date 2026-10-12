"""Build metadata/za_signal_run2.json from metadata/za_signal_run3.json (Run 2 UL NanoAODv9).

Changes:
- branches: drop the Run 3-only ones (pT/scale/smear branches written by the Run 3 corrections,
  bsConstrainedPt*, Rho_fixedGridRhoAll), add fixedGridRhoFastjetAll, L1PreFiringWeight_*,
  Electron/Photon dEsigma*/seedGain.
- branch_map: none (NanoAODv9 already carries the EGM-corrected pT; the Rochester correction
  writes Muon.pt directly; Muon ptErr stays the NanoAOD ptErr, see analysis.load_events).
- weights: add l1_prefiring_sf; drop photon_id_sf (no Run 2 custom-ID scale factors yet).
- independent collections: EGM scale from EGM_ScaleUnc.json, EGM smearing from dEsigma (GeV);
  FNUF and Material kept; no muon variations (Rochester replicas not used, as in HZgamma).
"""
import copy, json, sys
src, dst = sys.argv[1:3]
c = json.load(open(src))
r = copy.deepcopy(c)
drop = {"Rho_fixedGridRhoAll", "Muon_bsConstrainedPt", "Muon_bsConstrainedPtErr", "Photon_pt_ScaleUp",
        "Photon_pt_ScaleDown", "Muon_corrected_pt", "Muon_scaleUp_pt", "Muon_scaleDown_pt",
        "Muon_smearUp_pt", "Muon_smearDown_pt", "LHEScaleWeight_Zero"}
add = ["fixedGridRhoFastjetAll", "L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn",
       "Electron_dEsigmaUp", "Electron_dEsigmaDown", "Electron_seedGain", "Photon_seedGain",
       "Photon_dEsigmaUp", "Photon_dEsigmaDown"]
br = []
for b in c["branches"] + add:
    if b not in drop and b not in br:
        br.append(b)
r["branches"] = br
r["branch_map"] = {"mc": [], "data": []}
w = r["systematics"]["weights"]
w.pop("photon_id_sf")
w["l1_prefiring_sf"] = {"type": "event", "method": "from_function",
                        "function": {"module_name": "higgs_dna.systematics.run2_systematics", "name": "l1_prefiring_sf"},
                        "modify_central_weight": True,
                        "requires_branches": ["L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn"]}
ic = r["systematics"]["independent_collections"]
for k in ("Muon_scale", "Muon_smear"):
    ic.pop(k)
for obj in ("Photon", "Electron"):
    for kind in ("scale", "smear"):
        ic["%s_%s" % (obj, kind)] = {"method": "from_function", "branch_modified": [obj, "pt"],
                                     "function": {"module_name": "higgs_dna.systematics.run2_systematics",
                                                  "name": "egm_%s_unc_run2_%s" % (kind, obj.lower())}}
r["samples"]["years"] = ["2016preVFP", "2016postVFP", "2017", "2018"]
json.dump(r, open(dst, "w"), indent=4)
print("branches %d -> %d; weights %s; ICs %s" % (len(c["branches"]), len(br), sorted(w), sorted(ic)))
