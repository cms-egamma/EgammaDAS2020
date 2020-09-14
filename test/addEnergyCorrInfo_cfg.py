import FWCore.ParameterSet.Config as cms

process = cms.Process("EGUSER")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis') 
options.register('nrThreads',8,options.multiplicity.singleton,options.varType.int,"number of threads to use")
options.register('isData',True,options.multiplicity.singleton,options.varType.bool,"is this data (sets GT to data and lumi sections)")
options.parseArguments()

process.source = cms.Source("PoolSource",
                           fileNames = cms.untracked.vstring(options.inputFiles),  
                          )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
if options.isData:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v20', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v6', '')


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(5000),
    limit = cms.untracked.int32(10000000)
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.options = cms.untracked.PSet(
    numberOfStreams = cms.untracked.uint32(options.nrThreads),
    numberOfThreads = cms.untracked.uint32(options.nrThreads),
    wantSummary = cms.untracked.bool(True)
)

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=False, 
                       era='2017-UL')    

process.path = cms.Path(process.egammaPostRecoSeq)

process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = cms.untracked.vstring(
                                      "keep *",
                                      #we drop all the intermediate products we produce
                                      #keeping just resulting slimmedElectrons and slimmedPhotons
                                      "drop *_*_*_EGUSER",
                                      "keep *_slimmedElectrons_*_*",
                                      "keep *_slimmedPhotons_*_*"
                                  )
                                  fileName = cms.untracked.string(options.outputFile),
                                  SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('path')),                                     
)

process.outPath = cms.EndPath(process.output)
