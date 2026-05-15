# H -> Za Angular Distributions

Use `plot_hza_angular_distributions.py` to make unit-area angular-shape overlays for
`H -> Z a -> l+ l- gamma gamma` signal samples.

Example:

```bash
python3 Plot/scripts/plot_hza_angular_distributions.py \
  --input /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_inputs_nominal \
  --masses 1,5,10,20,30 \
  --years 2022preEE,2022postEE,2023preBPix,2023postBPix,2024
```

If final BDT-score branches are available, add a score cut, for example:

```bash
python3 Plot/scripts/plot_hza_angular_distributions.py \
  --input /eos/home-p/pelai/HZa/root_P2Root/run3_bdt_scored_nominal \
  --masses 1,5,10,20,30 \
  --final-score-cut 0.95
```

The script first tries to build generator-level vectors from full nano-like
`GenPart_*` branches.  This is the preferred input for "gen level before
selection".  If full `GenPart` is not present, it falls back to the flattened
HiggsDNA branches written by `za_tagger_resolved.py`, such as
`GenHzaZLepMinus_*`, `GenHzaZLepPlus_*`, `GenALPPhoton1_*`, and
`GenALPPhoton2_*`.  Older pT-ordered aliases like `GenHzaZLeadLep_*`,
`GenHzaZSubleadLep_*`, `GenALPLeadPho_*`, and `GenALPSubleadPho_*` are still
supported as fallbacks.  Reconstructed-level plots use
`Z_lead_lepton_*`, `Z_sublead_lepton_*`, `ALP_lead_photon_*`, and
`ALP_sublead_photon_*`.

Physics interpretation:

* `cos(theta_Z)` probes the spin-1 Z helicity/polarization.  A longitudinally
  polarized Z gives an approximate `1 - cos^2(theta_Z)` shape before detector
  and selection effects.
* `cos(theta_a)` checks the spin-0 ALP decay.  At generator level,
  `a -> gamma gamma` should be isotropic, so `cos(theta_a)` should be
  approximately flat.  The generator-level convention uses the higher-pT ALP
  photon when flattened branches are used; for a spin-0 decay either photon is
  equivalent up to `cos(theta_a) -> -cos(theta_a)`, so also compare
  `abs(cos(theta_a))`.
* The `phi` or signed `phi` angle between the `Z -> ll` and `a -> gamma gamma`
  decay planes is the more CP-sensitive observable.
* Reconstructed-level non-flatness should not be immediately interpreted as
  intrinsic spin or CP structure.  Compare with generator level first, because
  lepton/photon acceptance, photon pT and eta thresholds, photon ID, and
  leading/subleading photon ordering can sculpt the shapes.

Main outputs include:

* `hza_costhetaZ_gen.pdf/png`
* `hza_costhetaA_gen.pdf/png`
* `hza_abscosthetaA_gen.pdf/png`
* `hza_phi_decayplanes_gen.pdf/png`
* `hza_signedphi_decayplanes_gen.pdf/png`
* `hza_cosThetaZ_gen.pdf/png`

Corresponding `*_reco` and `*_reco_final` plots are written when the necessary
branches and final-selection information are available.
