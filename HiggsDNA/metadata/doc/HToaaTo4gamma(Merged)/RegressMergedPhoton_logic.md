# RegressMergedPhoton — code logic summary

End-to-end ML reconstruction of merged a → γγ candidates from CMS ECAL hits.
Original author: Steven Clark / Michael Andrews (CMU). Run3 port: 2026-05-26.

Source tree (after Run3 port into `CMSSW_15_0_14`):

```
CMSSW_15_0_14/src/
├── MLDataFormats/EgammaCandidates/                  ← reco::MLPhoton EDM class
│   ├── interface/MLPhoton.h
│   ├── interface/MLPhotonFwd.h
│   ├── src/MLPhoton.cc
│   ├── src/classes.h          src/classes_def.xml   ← ROOT dictionary
│   └── BuildFile.xml
└── RecoEgamma/EgammaMLPhotonProducers/              ← the EDProducer + plugins
    ├── interface/Cluster.h           src/Cluster.cc          ← clustering helper
    ├── interface/MLPhotonProducer.h  src/MLPhotonProducer.cc ← main producer
    ├── src/MLPhotonTablePlugin.cc                            ← NanoAOD flat-table plugin
    ├── data/classifier.onnx          data/regressor.onnx     ← Run2 UL18 trained models
    ├── test/Prod_MLNanoAOD_run3_mc.py                        ← cmsRun cfg for Run3 MC
    └── BuildFile.xml
```

---

## 1. The data class — `reco::MLPhoton`

Inherits `reco::RecoCandidate`. Per-candidate payload:

| Member | Setter / Getter | Meaning |
|---|---|---|
| `p4_` (PolarLorentzVector) | inherited from `RecoCandidate` | reconstructed 4-momentum |
| `moe_` | `setMassEnergyRatio` / `massEnergyRatio()` | **regressor output**: m_Γ / E |
| `diphoton_score_` | `setDiphotonScore` / `diphotonScore()` | classifier softmax for the "a → γγ" class |
| `monophoton_score_` | `setMonophotonScore` / `monophotonScore()` | classifier softmax for the "single photon" class |
| `hadron_score_` | `setHadronScore` / `hadronScore()` | classifier softmax for the "hadronic" class |
| `pfIso_` | `setPFIsolation` / `pfIsolation()` | E_Γ / (E_Γ + PF energy in ΔR<0.3) |
| `r1_, r2_, r3_` | `setR{1,2,3}` / `r{1,2,3}()` | shower-shape ratios computed from the crystal image |

Edm-collection alias: `reco::MLPhotonCollection = std::vector<reco::MLPhoton>`.

---

## 2. Inputs the producer consumes (set in the cfg)

From MiniAOD only:

| cfg key | InputTag (MiniAOD) | Why |
|---|---|---|
| `clusterInputTag` | `reducedEgamma:reducedEBEEClusters` | seed `reco::CaloCluster` list (not actually used in the current implementation — see note below) |
| `HEBInputTag` | `reducedEgamma:reducedEBRecHits` | ECAL barrel `EcalRecHit` collection — **the actual raw input** |
| `HEEInputTag` | `reducedEgamma:reducedEERecHits` | ECAL endcap (consumed but currently unused; analysis is barrel-only) |
| `vtxInputTag` | `offlineSlimmedPrimaryVertices` | primary vertex (`pvtx.z` and `pvtx.position`) |
| `pfCandInputTag` | `packedPFCandidates` (`pat::PackedCandidate`) | PF cone-isolation energy sum |
| `classifierPath` | `data/classifier.onnx` | 3-class softmax ONNX model |
| `regressorPath` | `data/regressor.onnx` | scalar regression ONNX model |
| ES (auto) | `CaloGeometryRecord → CaloGeometry` | crystal→(η,φ) geometry |

> Note: although `clustersToken` and `HEEToken` are declared, the actual algorithm builds clusters from scratch out of EB RecHits and only fills crystal images on the barrel. The cluster + EE tokens are kept for forward-compat / future endcap support.

---

## 3. Algorithm — `MLPhotonProducer::produce(Event&, EventSetup&)`

### Step 1. Read EB RecHits and turn each non-zero crystal into a `Cluster` seed

```cpp
for (auto ic = HEB->begin(); ic != HEB->end(); ++ic)
{
    if (ic->energy() > 0.) {
        EBDetId ebDID(ic->detid());
        auto cell = EBgeom->getGeometry(ebDID);
        Cluster tC(cell.eta, cell.phi,
                   ebDID.ieta(), ebDID.iphi(),
                   ic->energy(),
                   ebDID.isNextToBoundary(...));
        clusters.push_back(tC);
    }
}
```

Each `Cluster` starts as a single crystal carrying:
- `vec` — a (E/cosh η, η, φ) 3-vector.
- arrays `Es, ietas, iphis, ncracks` (length 1 initially).

### Step 2. Iteratively merge nearby clusters until fixed point

```cpp
while (sizeBefore != sizeAfter) {
    sort(clusters, by |vec|² descending);
    clusters = DoPairings(clusters, MATCH_DR=0.15);
}
```

`DoPairings(inC, R)`:
- pop the most-energetic remaining cluster `C1`,
- find the **nearest** other cluster `C2` (`FindNearest`),
- if `ΔR(C1, C2) ≤ R=0.15` → **`C1.combine(C2)`** (vector add `vec`, concatenate `Es, ietas, iphis, ncracks`),
  else keep both as separate output clusters.

This is a greedy ECAL "single-link" clustering tuned to keep merged-photon footprints together (cone matches AN-20-142 description of "shower-merged" topology).

### Step 3. For each surviving cluster — build an image and run the ML models

For each cluster `C` after merging:

#### 3a. `C.makeImage()` — produce a 30×30 (= `isize²`) crystal image
- Translate `ietas, iphis` so the seed lands at the centre of a `(isize/2, isize/2)` grid;
- Handle φ-wrap when crystals straddle iφ=360;
- Clip out-of-bounds crystals to the image edge;
- Replicate the flattened 900-element image `isize²` times in `image` (the model is fed 900 copies — quirk of the original ONNX export, preserved as-is for compatibility).

#### 3b. Run the **classifier**
```cpp
class_outputs = ortClassifier.run(class_input_names, img, {}, class_output_names, 1)[0];
```
- input: 30×30 image (replicated 900×),
- output: 3 raw logits.

After softmax:
```cpp
mlpho.setMonophotonScore( exp(class_outputs[0]) / Σ exp );
mlpho.setDiphotonScore  ( exp(class_outputs[1]) / Σ exp );
mlpho.setHadronScore    ( exp(class_outputs[2]) / Σ exp );
```

#### 3c. Run the **regressor**
```cpp
regress_data = { normalised_image_30x30,  {η} };
regress_outputs = ortRegressor.run(regress_input_names, regress_data, {}, regress_output_names, 1)[0];
float moe = regress_outputs[0];          // m_Γ / E_Γ
mlpho.setMassEnergyRatio(moe);
```

The image is L1-normalised before being fed to the regressor (sum to unity).
The regressor also takes the cluster η as an extra scalar input.

#### 3d. Build the 4-vector from `moe` + raw `energy`

`calculateLorentzVector(moe, energy, eta, phi, zpv)`:
- starts from the calorimeter η of the cluster (`eta`),
- reprojects onto a ray from the primary vertex z-coordinate `zpv` to the ECAL barrel face (R=129 cm), so the η used in the 4-vector is **vertex-corrected** (`etaprime`).
- `m = moe * energy`, `p = √(E² − m²)`, then build `PtEtaPhiM`.

#### 3e. PF cone isolation (ΔR=0.3)

```cpp
for (auto pf : *pfCand)
    if (ΔR(pf, Γ) < 0.3) pfCandE += pf.energy();
if (pfCandE < energy) pfCandE = energy;       // floor to Γ energy
mlpho.setPFIsolation(energy / pfCandE);
```
Higher = more isolated; equals 1 when the cluster alone explains the cone energy.

#### 3f. Shower-shape ratios r1, r2, r3

```cpp
r_n = Σ_i (xᵢ² + yᵢ²)^(n/2)   /   Σ_i Eᵢ
```
(`x, y` = crystal coordinates in the 30×30 image, `n=1,2,3`.) Cheap proxies of the spatial moments of the energy distribution that the ONNX classifier already saw internally — exposed for offline use.

### Step 4. Push back & write

Each cluster yields one `reco::MLPhoton`; they are stored as `reco::MLPhotonCollection` via `iEvent.put`.

---

## 4. Producing a NanoAOD-style flat table — `MLPhotonTablePlugin.cc`

```cpp
typedef SimpleFlatTableProducer<reco::MLPhoton> SimpleMLPhotonFlatTableProducer;
DEFINE_FWK_MODULE(SimpleMLPhotonFlatTableProducer);
```

This registers a `SimpleFlatTableProducer` specialised on `reco::MLPhoton`, so NanoAOD `Var("massEnergyRatio()", float)` expressions resolve against the **derived** class instead of the `reco::Candidate` base. Without it the table producer fails with “no method named diphotonScore found for type reco::Candidate”.

The cfg uses it to expose the table to NanoAOD:
```python
process.mlphotonsTable = cms.EDProducer("SimpleMLPhotonFlatTableProducer",
    src   = cms.InputTag("mlphotons", "mlphotons"),
    name  = cms.string("MLPhoton"),
    cut   = cms.string(""),
    variables = cms.PSet(P4Vars,
        massEnergyRatio = Var("massEnergyRatio()", float),
        diphotonScore   = Var("diphotonScore()",   float),
        monophotonScore = Var("monophotonScore()", float),
        hadronScore     = Var("hadronScore()",     float),
        pfIsolation     = Var("pfIsolation()",     float),
        r1 = Var("r1()", float),
        r2 = Var("r2()", float),
        r3 = Var("r3()", float),
    ))
```

Output NanoAOD branches:
```
nMLPhoton
MLPhoton_pt, _eta, _phi, _mass               (from P4Vars)
MLPhoton_massEnergyRatio
MLPhoton_diphotonScore, _monophotonScore, _hadronScore
MLPhoton_pfIsolation
MLPhoton_r1, _r2, _r3
```

---

## 5. cmsRun cfg — `test/Prod_MLNanoAOD_run3_mc.py`

| Knob | Value |
|---|---|
| era | `Run3_2024` (was `Run2_2018`) |
| GlobalTag | `150X_mcRun3_2024_realistic_v2` (matches `RunIII2024Summer24MiniAODv6`) |
| NanoAOD customise | `nanoAOD_customizeCommon` (CMSSW_15 replacement for the removed `nanoAOD_customizeMC`) |
| InputTag process | left empty, e.g. `cms.InputTag('reducedEgamma', 'reducedEBEEClusters')`, so CMSSW resolves against whatever PAT process wrote the MiniAOD |
| Schedule | `mlphotons → mlphotonsTable → nanoSequenceMC → endOfProcess → NANOAODSIM output` |

The result is a NanoAODv15-shaped file with the standard `Photon`, `Electron`, `Muon`, `GenPart`, … tables **plus** the new `MLPhoton` table — exactly the input HiggsDNA expects.

---

## 6. Differences from the original Run2 UL18 package

- **C++ Event-Setup API modernised**: legacy `iSetup.get<CaloGeometryRecord>().get(handle)` is removed in CMSSW ≥ 13; replaced by `esConsumes<CaloGeometry, CaloGeometryRecord>` token + `iSetup.getData(token)`.
- **Added `<use name="PhysicsTools/NanoAOD"/>`** to `BuildFile.xml` for the flat-table specialisation.
- **New plugin `MLPhotonTablePlugin.cc`** registering `SimpleFlatTableProducer<reco::MLPhoton>` — necessary because the generic `SimpleCandidateFlatTableProducer` cannot expand `diphotonScore()` / `massEnergyRatio()` member-functions on the base type.
- **New cfg** with Run3_2024 era, 150X GT, and `nanoAOD_customizeCommon` (`nanoAOD_customizeMC` was removed in CMSSW_15).
- ONNX models (`classifier.onnx`, `regressor.onnx`) are still the **Run2 UL18-trained ones** — not yet retrained on Run3 sim. They will need Z→ee tag-and-probe scale/smearing validation (AN §5.2-5.3) before being trusted for Run3.

---

## 7. Output file size and timing (20-event Run3 test, 26 May 2026)

- Input file: one MiniAODv6 file (~15k events; the test ran on the first 20).
- Wall time: ~26 s end-to-end (xrootd open + cmsRun + LZMA compression).
- Output `MLPhoton` table populated for all 20 events; `nMLPhoton` ranges 0–8 per event; both classifier and regressor pass through their full softmax / regression paths.

---

## 8. Where this fits in the analysis chain

```
MiniAODv6
  │
  │  cmsRun Prod_MLNanoAOD_run3_mc.py        ←  this package
  ▼
custom NanoAOD (Photon_* + GenPart_* + …  +  MLPhoton_*)
  │
  │  HiggsDNA  ZaTaggerRun3   (za_tagger_merged.py)
  ▼
parquet with merged-category branches + MergedML_residual_mass
                                          (ROI = m_Γ,regressed − m_a,gen)
  │
  │  Parquet2Root + plotting
  ▼
ROI residual histograms per mA hypothesis (AN-20-142 Fig 9 analogue)
```
