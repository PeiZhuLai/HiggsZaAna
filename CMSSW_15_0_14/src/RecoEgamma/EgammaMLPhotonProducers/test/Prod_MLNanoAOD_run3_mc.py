"""
MLPhoton + NanoAODv15-style production for Run3 2024 MC MiniAODv6.

Inputs : RunIII2024Summer24MiniAODv6, GlobalTag 150X_mcRun3_2024_realistic_v2
Outputs: standard NanoAODSIM tables + an additional MLPhoton table populated
         from the MLPhotonProducer (regressor + classifier).

Usage:
    cmsRun Prod_MLNanoAOD_run3_mc.py \\
        inputFiles=root://cms-xrd-global.cern.ch//store/...MINIAODSIM.root \\
        outputFile=mlnano_mA0p4.root \\
        maxEvents=100
"""

import os
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing("analysis")
options.register("year", "2024",
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Data-taking year (2022preEE/2022postEE/2023preBPix/2023postBPix/2024)")
options.parseArguments()

# ---- era + GT selection based on year ----
if options.year.startswith("2024"):
    from Configuration.Eras.Era_Run3_2024_cff import Run3_2024 as ERA
    _GT = "150X_mcRun3_2024_realistic_v2"
elif options.year.startswith("2023"):
    from Configuration.Eras.Era_Run3_2023_cff import Run3_2023 as ERA
    _GT = "130X_mcRun3_2023_realistic_v15" if "preBPix" in options.year else "130X_mcRun3_2023_realistic_postBPix_v6"
elif options.year.startswith("2022"):
    from Configuration.Eras.Era_Run3_cff import Run3 as ERA
    _GT = "130X_mcRun3_2022_realistic_v5" if "preEE" in options.year else "130X_mcRun3_2022_realistic_postEE_v6"
else:
    raise RuntimeError(f"Unsupported year: {options.year}")

process = cms.Process("MLNANO", ERA)

process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.NanoAOD.nano_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(*options.inputFiles),
    secondaryFileNames=cms.untracked.vstring(),
)
process.options = cms.untracked.PSet()
process.configurationMetadata = cms.untracked.PSet(
    annotation=cms.untracked.string("MLNANO nevts:{}".format(options.maxEvents)),
    name=cms.untracked.string("Applications"),
    version=cms.untracked.string("$Revision: 1.0 $"),
)

# ----------------------------------------------------------------------
# MLPhotonProducer
# Note: InputTag process name is left empty so CMSSW auto-resolves
#       it against whatever PAT process the MiniAOD was written with.
# ----------------------------------------------------------------------
CMSSW_BASE = os.environ["CMSSW_BASE"]
process.mlphotons = cms.EDProducer(
    "MLPhotonProducer",
    collectionLabel=cms.string("mlphotons"),
    classifierPath=cms.string(CMSSW_BASE + "/src/RecoEgamma/EgammaMLPhotonProducers/data/classifier.onnx"),
    regressorPath=cms.string(CMSSW_BASE + "/src/RecoEgamma/EgammaMLPhotonProducers/data/regressor.onnx"),
    clusterInputTag=cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
    HEEInputTag=cms.InputTag("reducedEgamma", "reducedEERecHits"),
    HEBInputTag=cms.InputTag("reducedEgamma", "reducedEBRecHits"),
    pfcandInputTag=cms.InputTag("packedPFCandidates"),
    vtxInputTag=cms.InputTag("offlineSlimmedPrimaryVertices"),
    pfCandInputTag=cms.InputTag("packedPFCandidates"),
)

from PhysicsTools.NanoAOD.common_cff import P4Vars, Var
process.mlphotonsTable = cms.EDProducer(
    "SimpleMLPhotonFlatTableProducer",
    src=cms.InputTag("mlphotons", "mlphotons"),
    name=cms.string("MLPhoton"),
    doc=cms.string("MLPhoton candidates: regressed m_Gamma + classifier scores"),
    singleton=cms.bool(False),
    cut=cms.string(""),
    variables=cms.PSet(
        P4Vars,
        massEnergyRatio=Var("massEnergyRatio()", float, doc="regressed m_Gamma / energy"),
        diphotonScore=Var("diphotonScore()", float, doc="Diphoton classifier score"),
        monophotonScore=Var("monophotonScore()", float, doc="Single-photon classifier score"),
        hadronScore=Var("hadronScore()", float, doc="Hadronic classifier score"),
        pfIsolation=Var("pfIsolation()", float, doc="PF cone-0.3 isolation ratio"),
        r1=Var("r1()", float, doc="shower-shape ratio r1"),
        r2=Var("r2()", float, doc="shower-shape ratio r2"),
        r3=Var("r3()", float, doc="shower-shape ratio r3"),
    ),
)

process.NANOAODSIMoutput = cms.OutputModule(
    "NanoAODOutputModule",
    compressionAlgorithm=cms.untracked.string("LZMA"),
    compressionLevel=cms.untracked.int32(9),
    dataset=cms.untracked.PSet(
        dataTier=cms.untracked.string("NANOAODSIM"),
        filterName=cms.untracked.string(""),
    ),
    fileName=cms.untracked.string(options.outputFile or "mlnano_run3.root"),
    # Slim output: only keep MLPhoton flat table; run/luminosityBlock/event are
    # part of the Events tree skeleton and are written automatically.
    # Standard NanoAOD content (Photon, Electron, Muon, Jet, HLT, L1, ...) is
    # NOT duplicated — downstream HiggsDNA reads those from the existing
    # NanoAODv15 files via friend tree on (run, luminosityBlock, event).
    outputCommands=cms.untracked.vstring(
        "drop *",
        "keep nanoaodFlatTable_mlphotonsTable_*_*",
    ),
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, _GT, "")

process.mlphotons_step = cms.Path(process.mlphotons)
process.mlphotonsTable_step = cms.Path(process.mlphotonsTable)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Slim production: skip nanoSequenceMC entirely. We only need MLPhoton; the
# standard NanoAOD tables (Photon/Electron/Muon/Jet/HLT/L1/...) live in the
# upstream NanoAODv15 file that HiggsDNA reads as the main input.
process.schedule = cms.Schedule(
    process.mlphotons_step,
    process.mlphotonsTable_step,
    process.endjob_step,
    process.NANOAODSIMoutput_step,
)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeCommon
process = nanoAOD_customizeCommon(process)

process.add_(cms.Service("InitRootHandlers", EnableIMT=cms.untracked.bool(False)))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
