# Modification of PhysicsTools.PatAlgos.patTemplate_cfg

# IMPORTS
import FWCore.ParameterSet.Config as cms			# Standard config file information.
#from PhysicsTools.PatAlgos.tools.pfTools import *		# Deep in here (PhysicsTools.PatAlgos.tools.jetTools) is addJetCollection.
from Configuration.AlCa.autoCond import autoCond		# For automatically determining global tags

# CONSTRUCT
process = cms.Process("fatjets")

# DEFINE CONFIGURATION
# Reports
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options = cms.untracked.PSet(
	wantSummary = cms.untracked.bool(True)			# Long summary after job.
)
# Input
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring()
)
# Output
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
	fileName = cms.untracked.string('out.root'),
	SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),	# Save only events passing the full pathway
	outputCommands = cms.untracked.vstring('drop *', "keep *_*_*_fatjets")	# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideSelectingBranchesForOutput
#	outputCommands = cms.untracked.vstring('drop *', *patEventContent )
)
process.outpath = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService",
	fileName = cms.string('out2.root')
)
# Events
process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(100)
)

# PRODUCTION
## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string(autoCond['startup'])		# Set global tags automatically.
#process.GlobalTag.globaltag = cms.string('START53_V27::All')		# Set global tags by hand.
##process.load('RecoJets.JetProducers.caSubjetFilterGenJets_cfi')
##process.caSubjetFilterGenJets.verbose = cms.untracked.bool(False)
##process.caSubjetFilterGenJets.src     = 'genParticles'	# Input module label
##process.caSubjetFilterGenJets.Rho_EtaMax = cms.double(5.0)
##process.caSubjetFilterGenJets.rParam = cms.double(0.5)

process.analyzer = cms.EDAnalyzer('Fatjets',
	v = cms.bool(False),		# verbose mode
	R = cms.double(0.7),
	make_gen = cms.bool(True),
	make_pf = cms.bool(True),
	do_true = cms.bool(False),
	do_prune = cms.bool(False)
)

## ------------------------------------------------------
#  NOTE: you can use a bunch of core tools of PAT to
#  taylor your PAT configuration; for a few examples
#  uncomment the lines below
## ------------------------------------------------------
#from PhysicsTools.PatAlgos.tools.coreTools import *

## remove MC matching from the default sequence
# removeMCMatching(process, ['Muons'])
# runOnData(process)

## remove certain objects from the default sequence
# removeAllPATObjectsBut(process, ['Muons'])
# removeSpecificPATObjects(process, ['Electrons', 'Muons', 'Taus'])


## let it run

process.p = cms.Path(
	  process.analyzer
)
## ------------------------------------------------------
# Commonly Changed User Parameters:
## ------------------------------------------------------

# PARAMETERS
msg = 850
msq = 250
R = 1.4
#process.analyzer.v = True
process.analyzer.v = False
process.analyzer.R = R
process.analyzer.make_gen = True
process.analyzer.make_pf = False
process.analyzer.do_true = True
process.analyzer.do_prune = False
# INPUT FILES
process.source.fileNames = cms.untracked.vstring(
	'file:/cms/tote/jet/reco/mgo{0}_msq{1}.root'.format(msg, msq),
)
# OUTPUT FILES
process.out.fileName = 'test/fatjets_test.root'
output_string = ""
if (process.analyzer.do_true == True):
	output_string += "_true"
if (process.analyzer.do_prune == True):
	output_string += "_pruned"
process.TFileService.fileName = 'test/test.root'
#process.TFileService.fileName = 'test/fatjets_sgtosq_msg{0}_msq{1}_R{2}{3}.root'.format(msg, msq, int(R*10), output_string)
# NUMBER OF EVENTS
process.maxEvents.input = 1000
# REPORTS
process.options.wantSummary = cms.untracked.bool(False)		# Turn off the long job summary.

