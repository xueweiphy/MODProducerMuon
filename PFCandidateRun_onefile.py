import sys
import os
from subprocess import call
from utilities.Analyzer_Onefile import *

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList

from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets



input_dir = sys.argv[2]
output_dir = sys.argv[3]

input_fileList = sys.argv[4]


map_file_path = sys.argv[5]
process_from_the_beginning = (sys.argv[6] == "1")

files_to_process = file_Names( input_fileList)
print files_to_process

dataFlag= sys.argv[7]

segments = files_to_process[0].split("/")
input_file = segments[len(segments) - 1:len(segments)][0]


readFiles = cms.untracked.vstring()
readFiles.extend(files_to_process)


process = cms.Process("MITCMSOpenData")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')

process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'


#process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db')
#process.GlobalTag.globaltag = cms.string('START53_LV6A1::All')




process.source = cms.Source("PoolSource", fileNames=readFiles)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

if ( dataFlag == "data" ) :
   goodJSON = "file_paths/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt"
   myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
   process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
   process.source.lumisToProcess.extend(myLumis)


process.ak5PFJets = ak5PFJets.clone(doAreaFastjet = cms.bool(True))
		    	
process.kt6PFJetsForIsolation = kt4PFJets.clone(rParam = 0.6, doRhoFastjet = True)

process.PFCandidateProducer = cms.EDProducer("PFCandidateProducer",
					rho = cms.InputTag("kt6PFJets","rho"),
					PFCandidateInputTag = cms.InputTag("particleFlow"),
					AK5PFInputTag = cms.InputTag("ak5PFJets"),
					mapFilename = cms.string(map_file_path),
					processFromTheBeginning = cms.bool(process_from_the_beginning),
					inputFile = cms.string(input_file),
					outputDir = cms.string(output_dir), 
					primaryVertices = cms.InputTag("offlinePrimaryVertices"),
					dataVersion = cms.string("5"),
               isMCarlo        = cms.untracked.bool(False),
					dataSet = cms.string("DoubleMu"), 
               skim = cms.untracked.bool(True)
				)
				
process.producer = cms.Path( process.ak5PFJets * process.kt6PFJetsForIsolation * process.PFCandidateProducer)
process.schedule = cms.Schedule( process.producer )
