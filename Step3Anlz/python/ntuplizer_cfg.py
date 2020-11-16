import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.convertHGCalDigisSim_cff import convertHGCalDigisSim
from Configuration.Eras.Era_Phase2_cff import Phase2

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing('analysis')
options.parseArguments()

process = cms.Process('ntuplizer')

from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('PROD',Phase2C9)
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.TFileService = cms.Service( "TFileService",
                                    fileName = cms.string(options.outputFile),
                                    closeFileFast = cms.untracked.bool(True) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles[0])
)

process.ntuplizer = cms.EDAnalyzer( 'SinglePhotonSpatialResolutionNtuplizer',
                                    hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
                                    hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
                                    hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
                                    simVertexes    = cms.InputTag("g4SimHits", "") )

process.p = cms.Path(process.ntuplizer)
