# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --no_exec --mc --python_filename run_crab.py --fileout NanoAODv2.root --eventcontent NANOAODSIM --datatier NANOAODSIM --step NANO -n 6284 --conditions 106X_upgrade2018_realistic_v15_L1v1 --era Run2_2018,run2_nanoAOD_106Xv1
import os

import FWCore.ParameterSet.Config as cms

# from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
from PhysicsTools.NanoAOD.common_cff import Var, ExtVar, CandVars

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv1_cff import run2_nanoAOD_106Xv1

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register('year',
                 '2018',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Year to process')

options.parseArguments()

inputFiles = options.inputFiles
outputFile = options.outputFile
maxEvents = options.maxEvents


process = cms.Process('NANO',Run2_2018,run2_nanoAOD_106Xv1)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(*inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:6284'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


# ML photons
process.mlphotons = cms.EDProducer("MLPhotonProducer",
    collectionLabel = cms.string("mlphotons"),
    classifierPath = cms.string("RecoEgamma/EgammaMLPhotonProducers/data/classifier.onnx"),
    regressorPath = cms.string("RecoEgamma/EgammaMLPhotonProducers/data/regressor.onnx"),
    clusterInputTag = cms.InputTag('reducedEgamma', 'reducedEBEEClusters', 'PAT'),
    HEEInputTag = cms.InputTag('reducedEgamma', 'reducedEERecHits', 'PAT'),
    HEBInputTag = cms.InputTag('reducedEgamma', 'reducedEBRecHits', 'PAT'),
    pfcandInputTag = cms.InputTag('packedPFCandidates', '', 'PAT'),
    vtxInputTag = cms.InputTag('offlineSlimmedPrimaryVertices', '', 'PAT'),
    pfCandInputTag = cms.InputTag('packedPFCandidates', '', 'PAT')
    # PhoInputTag = cms.InputTag('slimmedPhotons', '', 'PAT'),
)

# Define the mlphotonsTable module
# from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
process.mlphotonsTable = cms.EDProducer(
    'SimpleCandidateFlatTableProducer',
    src = cms.InputTag('mlphotons', 'mlphotons'),
    name = cms.string('MLPhoton'),
    doc = cms.string('Diphoton Objects and Tagging Variables'),
    singleton = cms.bool(False), # the number of entries is variable
    cut = cms.string(''),
    variables = cms.PSet(
        pt = Var("pt", float, precision=-1),
        eta = Var("eta", float, precision=12),
        phi = Var("phi", float, precision=12),
        massEnergyRatio = Var("massEnergyRatio()", float, doc="Regressed mass/energy"),
        diphotonScore = Var("diphotonScore()", float, doc="Diphoton Classifier score"),
        monophotonScore = Var("monophotonScore()", float, doc="Single Photon Classifier score"),
        hadronScore = Var("hadronScore()", float, doc="Hadronic Classifier score"),
        pfIsolation = Var("pfIsolation()", float, doc="Ratio of mlphoton energy to sum of PF candidates energy in a cone of 0.3 around the mlphoton)"),
        r1 = Var("r1()", float, doc="IDK"),
        r2 = Var("r2()", float, doc="IDK"),
        r3 = Var("r3()", float, doc="IDK"),
    ),
)

# process.dump = cms.EDAnalyzer('EventContentAnalyzer')

# Output definition
process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAODSIM'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(outputFile),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# process.NANOAODSIMoutput.outputCommands.extend([
#     'keep *_MLPhotons_*_*',
# ])


# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1', '')

# Path and EndPath definitions
process.mlphotons_step = cms.Path(process.mlphotons)
process.mlphotonsTable_step = cms.Path(process.mlphotonsTable)
# process.dump_step = cms.Path(process.dump)
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.mlphotons_step,
                                # process.dump_step,
                                process.mlphotonsTable_step,
                                process.nanoAOD_step,
                                process.endjob_step,
                                process.NANOAODSIMoutput_step)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
