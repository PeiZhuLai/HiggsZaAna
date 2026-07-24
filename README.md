# HiggsZaAna

This is the full analysis chain of the **Run3 (2022 / 2023 / 2024)** **H &rarr; Z(&rarr;&#8467;&#8467;) + a** search, where `a` is a light axion-like particle (ALP) decaying to two photons (`a` &rarr; &gamma;&gamma;). At low ALP mass the two photons are reconstructed either as a **resolved** di-photon pair or as a single **merged** photon, so both regimes are treated in this repository. It is a combination of HiggsDNA (tagging), machine learning (categorization), spurious signal test and final fit.

> This repository is a fork of [HiggsZGammaAna](https://github.com/Vvvvvvvictor/HiggsZGammaAna) adapted to the H &rarr; Z + ALP analysis.

```
git clone https://github.com/PeiZhuLai/HiggsZaAna.git --recursive
```

## HiggsDNA for tagging
**1. Setup environment**

The conda env would occupy a large quota. Make sure that you have a few GB left in your path (e.g. `/eos/user`). This step would take a few hours, please keep patient.
```
cd HiggsZaAna/HiggsDNA/
conda env create -f environment.yml -p <some_path_where_you_have_more_disk_space>/.conda/envs/
```
Then, activate the conda environment (the env name is defined in `environment.yml`)

```
cd HiggsZaAna/HiggsDNA
conda deactivate
conda activate higgs-zg-ana
```

You may also want to increase your disk quota at [this link](https://resources.web.cern.ch/resources/Manage/EOS/Default.aspx), otherwise you may run out of space while installing your `conda` environment.

One additional package, `correctionlib`, must be installed via `pip`, rather than `conda`. Run
```
source setup.sh
```
to install this script.

**2. Add HiggsDNA package**

Install it by:
```
pip install -e .
```
If you notice issues with the `conda pack` command for creating the tarball, try updating and cleaning your environment with (after activating the analysis env):
```
conda env update --file environment.yml --prune
```

**3. Run the tagger**

The samples are documented in the `metadata/za_*_run3.json` config files (e.g. `metadata/za_data_run3.json`, `metadata/za_signal_run3.json`, `metadata/za_bkgmc_run3.json`). If you want to analyse other samples, please modify these files.

The event selection and the cutflow are defined in the Za taggers:

- `higgs_dna/taggers/za_tagger_resolved.py` &mdash; resolved di-photon regime.
- `higgs_dna/taggers/za_tagger_merged.py` &mdash; merged single-photon regime (very low ALP mass).

Control-region / tag-and-probe taggers used for scale-factor validation live alongside them (`control_zee_zmmg_tagger.py`, `tnp_zee_tagger.py`, `tnp_zmmg_tagger.py`).

The core runner is `scripts/run_analysis.py`. Ready-to-use wrapper scripts are provided per task, for example:
```
bash scripts/run_ana_data.sh      # data
bash scripts/run_ana_signal.sh    # signal MC
bash scripts/run_ana_bkgmc.sh     # background MC
```
The merged-photon regime has its own wrappers (`scripts/run_merged_data.sh`, `scripts/run_merged_signal.sh`, `scripts/run_merged_bkgmc.sh`).

The output (parquet) is stored in the path documented at the top of each wrapper script. Please check it before running the jobs. The meaning of every option can be found in `scripts/run_analysis.py`.

**4. Some useful notes**

The logical operator precedence is: 1. `&`, 2. `|`. It is better to wrap each `&` in `()`.

If you have questions about usage or code structure, please look at the [HiggsDNA contents](https://sam-may.github.io/higgs_dna_tutorial.github.io/).

The dataset names can be found in the DAS system [(example)](https://cmsweb.cern.ch/das/request?instance=prod/global&input=file+dataset%3D%2FVBFHToZG_M-125_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer19UL17NanoAODv2-106X_mc2017_realistic_v8-v1%2FNANOAODSIM). Then you can put them into the `metadata/za_*_run3.json` config files.

If you are confused about the cross-section of a MC dataset, please get it through [this link](https://xsdb-temp.app.cern.ch/?columns=38289920&currentPage=0&ordDirection=1&ordFieldName=process_name&pageSize=50).

If you want to find the path of a dataset, please get it through [this link](https://cms-pdmv.cern.ch/mcm/requests?page=0&mcdb_id=102X_dataRun2_v11&shown=127).

If you want to get the golden json, these files, `/eos/user/c/cmsdqm/www/CAF/certification/Collisions*/Cert_Collisions*_*_*_Golden.json`, are recommended. Please put them into `metadata/golden_json/` and also change the golden json choice in the tagger you use.

You can find the name of variables in NanoAOD in [this link](https://cms-nanoaod-integration.web.cern.ch/autoDoc/NanoAODv11/2022postEE/doc_Muon_Run2022G-PromptNanoAODv11_v1-v2.html)

## Processed lumi check

The processed lumi information will be printed in the log file with the `[[INFO]] Processed Lumi:` prefix. The code `scripts/CheckProcessedLumi.py` can grab this information from the log files and combine it into a final file named `processedLumi.json`.

Then, run the lumi check command:

1. Load the CMS environment using `cmsenv`.
2. Set the brilcalc environment:
```
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda3/bin:$PATH
```
3. Install the latest version of brilcalc:
```
pip install --user --upgrade brilws
```
Sometimes you may get a strange error. Then try to uninstall brilcalc and reinstall it. Somehow it works many times.
```
pip uninstall -y brilws
```
4. Get the luminosity:
```
brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb  -i <JSON FILE PATH>/processedLumis.json
```

## Z Constrain Refit

The di-lepton kinematics are refit with a Z-mass constraint to improve the ALP mass resolution. See details in `Plotter`.

## Parquet &rarr; ROOT conversion

Before the ML / fit steps the HiggsDNA parquet output is converted to flat ROOT trees (adding the training inputs / BDT-score branches). This is handled in `Parquet2Rootfile/`:
```
bash Parquet2Rootfile/1_run_P2Root.sh
```
The input/output paths are set at the top of the script (`Parque2Root_za.py` is the core converter). Dedicated variants exist for the control samples and the normalizing-flow inputs (`Parque2Root_BDT_control.py`, `Parque2Root_za_NFlow.py`).

## Machine learning for categorization

The BDT training and application live in `HZaMVA/`. Two XGBoost models are trained to avoid the mass-sculpting effect at very low ALP mass:

- `model_Za_BDT_lowmass_run3.pkl` &mdash; low-mass region (mA = 1&ndash;3 GeV).
- `model_Za_BDT_highmass_run3.pkl` &mdash; high-mass region (mA = 4&ndash;30 GeV).

Setup the ML/plotting environment (XGBoost stack) as described in `HZaMVA/`.

**1. Sideband reweighting**

Data-driven sideband reweighting factors (applied consistently to signal and background before training) are produced by:
```
python HZaMVA/scripts/1_make_sideband_reweight.py
```

**2. Training a model**

The training driver reads the feature list and hyper-parameters from `HZaMVA/scripts/hza_features.py` and the `run3_Za_BDT.py` configuration, trains in four folds, and transforms the scores such that the unweighted signal distribution is flat:
```
cd HZaMVA/scripts
bash 2_train.sh            # runs run3_Za_BDT.py, writes model_Za_BDT_run3.pkl
```
Once you are happy with the model, copy it into the `using/` folder that the downstream steps read from:
```
bash 3_save_model.sh
```

**3. Applying the scores**

The BDT scores are written back into the ROOT trees during the Parquet &rarr; ROOT step (`Parquet2Rootfile/Parque2Root_za.py`), which loads the model from `HZaMVA/using/`.

**4. Optimizing the categories**

The category boundaries and the number-counting significances are optimized as part of the grand pipeline `1_grand.sh` (see the `MVA_CUT` / `DETERMINE_MVA_CUT` stages). For a single interactive scan you can use the boundary-scan helpers under `HZaMVA/scripts/`.

> The end-to-end resolved-photon workflow (parquet &rarr; ROOT &rarr; MVA cut &rarr; Tree2WS &rarr; background/signal fit &rarr; datacard &rarr; combine limits) is driven by `1_grand.sh`. The merged-photon low-mass workflow is driven by `1_grand_merged.sh`.

## Spurious signal test

For this step, there is no need to activate the conda analysis environment.

### Setup environment

Get the CMS env first. Please make sure you have entered the `HiggsZaAna` directory.
```
cd HiggsZaAna/SSTest/
scram project CMSSW CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
```
Two more packages are needed and must be cloned from github.
```
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit -b 112x-comb2021
git clone https://github.com/cms-analysis/CombineHarvester.git HiggsAnalysis/CombineHarvester -b v2.0.0
scram b -j 9
```
Please note that you need to initialize the environment **each time you set up a terminal**:
```
cd CMSSW_11_3_4/src
cmsenv
```

### Scripts to run a task
**1. Prepare the background and signal template**

We can get a file (`bkg_sig_template.root`) with the bkg and sig shape in each channel and category by running:
```
python SSTest/Generate_template.py
```
You need to modify the range of the histogram and make sure it is the same as in `SSTest.cpp` or `SSTest_core_function.cpp`.

**2. Run the spurious signal test and find the best bkg function**

Run this code:
```
root -l -q SSTest/SSTest.cpp
```
Some options in `SSTest.cpp`: `cat` = which category of this channel to test, `channel` = one_jet and similar, `sig` = how many times the signal is injected, `bkg_fun` = function used to describe the background.

A spurious signal test based on the core-function method is also provided:
```
root -l -q SSTest/SSTest_core_function.cpp
```
Options in `SSTest_core_function.cpp`: `cat` = which category to test, `channel`, `sig`. The background function for the fit can be changed inside the code.

**Take care that the observed variable in this code should be `CMS_hza_mass`, matching the name used in the input workspace. We must use the same name, or it will fail.**

We select the best background function via the spurious signal, the chi-square and the F-test. The fitting result and the name of the best function are written to the log file.

## Merged-photon low-mass analysis

The dedicated merged-photon (very low mA) statistical analysis &mdash; per-mA background/signal workspaces, datacards and combine limits &mdash; lives in `MergedAna/`. It must be run inside the flashgg combine `cmsenv`:
```
bash MergedAna/run_limits_merged.sh
```
The merged-photon energy regression used upstream is in `RegressMergedPhoton/`.

## Final Fit
### Setup environment
**1. Setup the CMSSW environment**
```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
```

**2. Install the dependency packages**
```
git clone https://github.com/jonathon-langford/HiggsAnalysis.git
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.2.0
```

**3. Compile the dependency packages**
```
cd ../
scramv1 b clean; scramv1 b
```

**4. Install the Flashgg Final Fit package**
```
cd final_fit/
mv flashggFinalFit CMSSW_10_2_13/src
```
Setup the CMSSW environment:
```
cd CMSSW_10_2_13/src
cmsenv
cd flashggFinalFit/
```

```
cd Signal
make clean
make
```

```
cd Background
make clean
make
```
