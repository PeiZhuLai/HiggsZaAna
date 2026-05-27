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

from Configuration.Eras.Era_Run3_2024_cff import Run3_2024

options = VarParsing("analysis")
options.register("year", "2024",
                 VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Data-taking year")
options.parseArguments()

process = cms.Process("MLNANO", Run3_2024)

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
    outputCommands=process.NANOAODSIMEventContent.outputCommands,
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "150X_mcRun3_2024_realistic_v2", "")

process.mlphotons_step = cms.Path(process.mlphotons)
process.mlphotonsTable_step = cms.Path(process.mlphotonsTable)
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

process.schedule = cms.Schedule(
    process.mlphotons_step,
    process.mlphotonsTable_step,
    process.nanoAOD_step,
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
