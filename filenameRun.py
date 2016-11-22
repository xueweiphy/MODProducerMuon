import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.PythonUtilities.LumiList as LumiList
import sys

data_file_link = sys.argv[2]
registry_file_path = sys.argv[3]
file_directory = data_file_link[7:len(data_file_link) - 41]
ifMC = sys.argv[4]

file_name = data_file_link[len(data_file_link) - 41:len(data_file_link)]

process = cms.Process("filenameMapProducer")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source ("PoolSource", fileNames=cms.untracked.vstring( data_file_link ) )


#goodJSON = "file_paths/JSON_2011.txt"
if (ifMC == "data") :
   goodJSON = "file_paths/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON.txt"
   myLumis = LumiList.LumiList(filename = goodJSON).getCMSSWString().split(',')
   process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
   process.source.lumisToProcess.extend(myLumis)


process.filenameMapProducer = cms.EDProducer("filenameMapProducer", 
						filename = cms.string(file_name), 
						file_dir = cms.string(file_directory), 
						outputFile = cms.string(registry_file_path) 
						)


process.MessageLogger = cms.Service("MessageLogger",
        				default   =  cms.untracked.PSet(
                                                     timespan = cms.untracked.int32(60)
       )                                                  
)       
						
process.producer = cms.Path(process.filenameMapProducer)
process.schedule = cms.Schedule( process.producer )
