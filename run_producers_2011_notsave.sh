#!/bin/bash

# Clean up everything and recompile.
#cmsenv
#scram b clean
#scram b



sharedF=/media/sf_weixue/
stepN=2
fileNum=-1
files0=list1
files1=./file_paths/DiMuon/${files0}.txt
#register=registry_Jet_${files0}.txt



python main_analyzer.py ${files1} ${sharedF} ${fileNum} ${stepN} 1 test0 data true






