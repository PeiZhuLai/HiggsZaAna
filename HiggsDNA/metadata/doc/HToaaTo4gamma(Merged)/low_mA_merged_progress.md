# Low-mA merged analysis — progress, conclusions, next steps

Reference: CMS AN-20-142 v6 (H → aa → 4γ, merged). Adaptation: **H → Zã, Z → ll, ã → γγ** with low mA (0.1–1 GeV); ã reconstructs as a **single merged photon Γ**.

Last updated: 2026-05-26.

---

## 1. What has been done

### 1.1 Terminology clarified
- `m_Γ ROI` is **not** a (0, 1.2) GeV mass window. It is the **regression residual / closure**: `m_Γ,regressed − m_a,gen` (AN §5.4, Fig 9).
- The (0, 1.2) GeV bound is the **m_a coverage** of the analysis, not the ROI cut.

### 1.2 Photon ID — AN-20-142 §4.2 added as 16th scenario
File: `HiggsDNA/higgs_dna/taggers/za_tagger_merged.py`

New constants in `DEFAULT_OPTIONS["photons"]` (block `AN2020_*`).
New photon-level flag `pass_phid_AN2020` computed from NanoAODv15 branches:

| AN cut | NanoAODv15 branch | Tier |
| --- | --- | --- |
| R9 > 0.5 | `Photon_r9` | always |
| H/E < 0.04596 (PU-corr) | `Photon_hoe_PUcorr` | always |
| no pixel seed | `Photon_pixelSeed` | always |
| EGM MVA > −0.96 | `Photon_mvaID` | always |
| I_tk / p_T,Γ < 0.07 | `Photon_trkSumPtHollowConeDR03 / Photon_pt` | always |
| barrel, \|η\| < 1.4 | `Photon_eta` | always |
| σ_iηiη < 0.015 | `Photon_sieie` | only R9 ≤ 0.85 |
| I_γ < 4 GeV (PU-corr) | `Photon_ecalPFClusterIso_PUcorr` | only R9 ≤ 0.85 |
| I_tk < 6 GeV | `Photon_trkSumPtHollowConeDR03` | only R9 ≤ 0.85 |

Nominal `id_cut` left unchanged (still `phid_custom_tight`) — AN2020 is study-only.

### 1.3 Gen-truth merged proxy (until regressor is integrated)
In `calculate_gen_info`:
- `GenALP_mass` — truth m_a (= ã mass).
- `GenALP_gammaL` — Lorentz boost `E_a / m_a`.
- `MergedGamma_*` — truth-matched single reco Γ (closest reco Photon to gen ã in ΔR):
  `MergedGamma_dR_to_genA, _pt, _eta, _phi, _r9, _mvaID, _hoe_PUcorr, _sieie, _pass_phid_AN2020`.

### 1.4 Merged category — nPhoton == 1 event flag and cutflow
In `produce_and_select_zgammas`:
- `n_photons_sel`, `has_1merged_cand` (any selected photon, ID-agnostic).
- `n_photons_AN2020`, `has_1merged_AN2020` (sum of `pass_phid_AN2020` == 1).
- Sequential cutflow `zgammas_merged_AN2020`(_w): `N_lep_sel → trig_cut → lep_pt_cut → has_z_cand → has_1merged`.
- Event flag `pass_allcuts_merged_AN2020`.

### 1.5 Submission scripts wired to merged configs
- `HiggsDNA/scripts/run_merged_signal.sh` → `metadata/za_merged_signal_run3.json`
- `HiggsDNA/scripts/run_merged_bkgmc.sh` → `metadata/za_merged_bkgmc_run3.json`
- `HiggsDNA/scripts/run_merged_data.sh` → `metadata/za_merged_data_run3.json`
- All three configs already point to `higgs_dna.taggers.za_tagger_merged` (`ZaTaggerRun3`).

### 1.6 RegressMergedPhoton **ported to Run3 / CMSSW_15_0_14**

- `cmsrel CMSSW_15_0_14` under `HiggsZaAna/`, SCRAM_ARCH `el9_amd64_gcc12`.
- Ported producer source from old Run2 UL package:
  - **C++ modernisation**: `iSetup.get<CaloGeometryRecord>().get(handle)` → `esConsumes<CaloGeometry, CaloGeometryRecord>` token + `iSetup.getData(token)`.
  - Added `<use name="PhysicsTools/NanoAOD"/>` to BuildFile.
  - **Registered** a `SimpleFlatTableProducer<reco::MLPhoton>` plugin → `SimpleMLPhotonFlatTableProducer` (needed because the generic `SimpleCandidateFlatTableProducer` only sees `reco::Candidate` methods and cannot expand `diphotonScore()`, `massEnergyRatio()`, etc.).
- New cfg `RecoEgamma/EgammaMLPhotonProducers/test/Prod_MLNanoAOD_run3_mc.py`:
  - era `Run3_2024`, GT `150X_mcRun3_2024_realistic_v2`.
  - Uses `nanoAOD_customizeCommon` (the v15 replacement for the deprecated `nanoAOD_customizeMC`).
  - InputTag process name left empty for MiniAOD auto-resolve.
- `scram b -j 4` clean, `edmPluginDump` shows both `MLPhotonProducer` and `SimpleMLPhotonFlatTableProducer` registered.
- **20-event end-to-end test on `mA_MINI_M0p4` MiniAOD: success.** Output NanoAOD carries the `MLPhoton_*` table:
  ```
  nMLPhoton, MLPhoton_pt, _eta, _phi, _mass, _massEnergyRatio,
  _diphotonScore, _monophotonScore, _hadronScore, _pfIsolation,
  _r1, _r2, _r3
  ```

### 1.7 HiggsDNA wiring of MLPhoton + ROI residual

In `za_tagger_merged.py`:

- **`produce_and_select_zgammas`** detects `events.MLPhoton` and adds:
  - `n_MLPhoton`, `n_MLPhoton_diphoton` (diphotonScore > 0.5), `has_1MLPhoton_diphoton`.
  - `MLPhoton_lead_{pt,eta,phi,mass,massEnergyRatio,diphotonScore,monophotonScore,hadronScore,pfIsolation,r1,r2,r3}` (leading by diphotonScore).
  - `zgammas_merged_ML`(_w) cutflow + `pass_allcuts_merged_ML`.

- **`calculate_gen_info`** truth-matches gen ã (pdgId=9000005) to the closest MLPhoton in ΔR and stores:
  - `MergedML_dR_to_genA, _pt, _eta, _phi, _mass, _massEnergyRatio, _diphotonScore, _monophotonScore, _hadronScore, _pfIsolation`.
  - **`MergedML_residual_mass = MLPhoton.mass − GenALP_mass`** ← this is the AN §5.4 ROI residual.

Three merged configs updated with new input branches (`MLPhoton_*`, `nMLPhoton`) and `variables_of_interest`.

**Logic verified** end-to-end on the 20-event test file: truth-matched events (ΔR<0.05) show residual ≈ ±0.03–0.08 GeV at gen mA = 0.4 GeV, consistent with AN Fig 9 (~200 MeV FWHM).

### 1.8 Condor production for mA = 0.4 GeV — dry-test submitted

- Submitter `RegressMergedPhoton/condor/submit_mass_point.sh`:
  - DAS-queries dataset, writes per-file args, writes JDL.
- Worker `RegressMergedPhoton/condor/run_one_file.sh`:
  - sources CMSSW_15_0_14, runs `cmsRun Prod_MLNanoAOD_run3_mc.py`, copies output to `/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/${mass_tag}/`.
- Generated JDL covers 33 files of `RunIII2024Summer24MiniAODv6 mA=0.4 GeV` (500k events).
- **First job submitted alone** (cluster `8309035`, `longlunch` flavour). Will widen to full 33 after it succeeds.

---

## 2. Conclusions

- HiggsDNA chain processes **NanoAODv15 mA_M0p\*** samples and produces:
  - the AN2020 photon-ID flag,
  - the nPhoton==1 merged-category flag and event-level cutflow,
  - truth-level m_a, γ_L, ΔR(γγ)_gen, plus a truth-matched reco-Γ proxy.
- RegressMergedPhoton end-to-end on Run3 2024 MiniAODv6 works in CMSSW_15_0_14, produces the `MLPhoton` table with the regressor mass.
- HiggsDNA reads that table and writes the ROI residual `MergedML_residual_mass` per event.
- Once the dry-test job lands, scale to full 33-file production for mA = 0.4 GeV, then the 8 other mass points.

---

## 3. Next steps

### NS-1 — finish the mA = 0.4 GeV production
- Wait for `8309035` to complete; verify EOS output `/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/mA_M0p4/mA_M0p4_0001.root` (~few hundred MB, NanoAOD-sized).
- If OK: `condor_submit` the full 33-job JDL.
- Plan-B: if held, inspect `condor_q -hold` and adjust memory/disk request.

### NS-2 — point HiggsDNA at the produced MLNanoAOD
- Edit `HiggsDNA/metadata/samples/zgamma_tutorial.json` so `mA_MINI_M0p4.files` lists the produced root files (e.g. `/eos/project/h/htozg-dy-privatemc/pelai/HZa/MLNanoAOD/mA_M0p4/mA_M0p4_*.root`).
- Run `HiggsDNA/scripts/run_merged_signal.sh` restricted to `mA_MINI_M0p4`.
- First plot: histogram of `MergedML_residual_mass` after `pass_allcuts_merged_ML`.

### NS-3 — sample-tag dispatch (Task #2 in tracker)
- In `za_tagger_merged.py`: when sample tag indicates low-mA merged (`mA_MINI_M0p*`), set nominal `id_cut = pass_phid_AN2020` and `has_2gamma_cand = has_1MLPhoton_diphoton`. Resolved samples (mA ≥ 2 GeV) keep current behaviour.

### NS-4 — Z→ee tag-and-probe ONNX validation (AN §5.2-5.3)
- Independent of the production. Use `DYJetsTo2E` MiniAOD + same Prod_MLNanoAOD_run3_mc.py.
- Extract probe-electron `MLPhoton_mass` distribution; fit data/MC scale + smearing.
- If scale offset > ~5% or smearing > 20 MeV, plan ONNX retraining on Run3 sim before trusting the residual closure.

### NS-5 — scale to 9 mass points
- For `mA_M0p1..mA_M0p9` (minus mA_M0p4 already running), call `submit_mass_point.sh` per mass point.
- Estimated time per job: ~3hr (cmsRun ~15k events/file).
- Total: 9 mass × 33 files × ~3hr = ~900 CPU-hr; ≲ 1 day wall-clock on condor.

### NS-6 — 1_grand.sh stage 0
- Add env-gated `RUN_MINI_TO_MLNANO=1` stage in `1_grand.sh` to drive `submit_mass_point.sh` for the mA list before HiggsDNA runs.

---

## 4. Files changed

```
A CMSSW_15_0_14/                                                       (new cmsrel workarea)
M CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/interface/MLPhotonProducer.h
M CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/src/MLPhotonProducer.cc
M CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/BuildFile.xml
A CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/src/MLPhotonTablePlugin.cc
A CMSSW_15_0_14/src/RecoEgamma/EgammaMLPhotonProducers/test/Prod_MLNanoAOD_run3_mc.py
A RegressMergedPhoton/condor/run_one_file.sh
A RegressMergedPhoton/condor/submit_mass_point.sh
A RegressMergedPhoton/condor/stage/mA_M0p4/{files.txt,args.txt,args_first1.txt,job.jdl,job_first1.jdl,logs/}
M HiggsDNA/higgs_dna/taggers/za_tagger_merged.py
M HiggsDNA/metadata/za_merged_signal_run3.json
M HiggsDNA/metadata/za_merged_bkgmc_run3.json
M HiggsDNA/metadata/za_merged_data_run3.json
M HiggsDNA/scripts/run_merged_signal.sh
M HiggsDNA/scripts/run_merged_bkgmc.sh
M HiggsDNA/scripts/run_merged_data.sh
A HiggsDNA/metadata/doc/HToaaTo4gamma(Merged)/low_mA_merged_progress.md   (this file)
```

---

## 5. Task tracker snapshot (2026-05-26)

| # | Status | Subject |
|---|--------|---------|
| 1 | done | Add phid_AN2020 photon ID scenario |
| 2 | pending | Add sample-tag-based ID dispatch (low_mass_merged vs nominal) |
| 3 | done | Truth-level merged ã path (gen info + acceptance) |
| 4 | done | RegressMergedPhoton Run3 upgrade (CMSSW_15_0_14) |
| 5 | pending | Integrate merged stage into 1_grand.sh |
| 6 | done | Add merged-category nPhoton==1 event flag + cutflow |
| 7 | done | Wire run_merged_*.sh to za_merged_*_run3.json |
| 8 | done | HiggsDNA wiring: read MLPhoton + compute m_Γ ROI residual |
| 9 | in_progress | Condor production for mA = 0.4 GeV MiniAOD → MLNanoAOD |
