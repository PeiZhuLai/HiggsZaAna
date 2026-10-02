import FWCore.ParameterSet.Config as cms
import sys
files = [l.strip() for l in open(sys.argv[-1]) if l.strip()]
process = cms.Process("GenXSec")
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring(["root://cms-xrd-global.cern.ch/"+f for f in files]))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))
process.genxsec = cms.EDAnalyzer("GenXSecAnalyzer")
process.p = cms.Path(process.genxsec)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000
