import os
from subprocess import call
import sys
from time import time
import copy

#list_of_downloads_file = sys.argv[1]
#destination_dir = sys.argv[2]
#fileNum = int(sys.argv[3] )

def download_root_files(list_of_file, destination, runFlag = 1):
   try :
      f = open(list_of_file, 'r')
      files_to_download = f.read().split("\n")
   except :
      files_to_download = list_of_file

   #destination_paths = []
   rootfileNs =[]
   for ii, root_link in enumerate(files_to_download):
      print "Downloading file #" + str(ii)
      http_link = root_link.replace("root://eospublic.cern.ch//", "http://opendata.cern.ch/")
      index0= http_link.rindex('/') + 1
      root_file_name = http_link[index0:len(http_link)]
      rootfileNs.append( root_file_name)
      destination_dir_path = destination + root_link[26:index0+1]
      if ( ii != 0 and destination_dir_path !=  destination_paths   ) :
         print "Warning! : the destination_dir_paths are different."
         
      #destination_paths = copy.deepcopy ( destination_dir_path ) 
      destination_paths = destination_dir_path +'/'
      #print destination_paths
      

      if runFlag == 1 :
         call(["wget", "--continue", http_link, "-N", "-P", destination_dir_path])

      print "\n"*4
   #print destination_paths
   return [destination_paths,  rootfileNs ]

def download_one_root_file(file1, destination):
   print "Downloading file #" 
   http_link = file1.replace("root://eospublic.cern.ch//", "http://opendata.cern.ch/")
   index0= http_link.rindex('/') + 1
   root_file_name = http_link[index0:len(http_link)]
   destination_dir_path = destination + file1[26:index0+1]
   #print "path ?" ,destination_dir_path 
   #print root_file_name
   call(["wget", "--continue", http_link, "-N", "-P", destination_dir_path])
   return [destination_dir_path , root_file_name]



def file_Names ( list_of_downloads_file ):
   f = open(list_of_downloads_file, 'r')
   files_to_download = f.read().split("\n")
      
   return files_to_download[:-1]




#fileName = file_Names( list_of_downloads_file ) 
#if ( fileNum != -1 ) :
   #fileName = fileName[:fileNum]
#
#for ii, ff in enumerate ( fileName ) :
   #start = time()
   #if ii == 0  :
      #start0 = start * 1 
   #print ff
   ### download file
   #AODpath ,fn0 = download_one_root_file(ff, destination_dir)
  # 
   ### create registry
   #registry_fileName = AODpath +'/' + fn0[:-5] +'_reg.txt'
  # 
   #reg_file_exist= os.path.isfile ( registry_fileName)  and os.path.getsize( registry_fileName) > 0
   #registry_file =  'file://'+ AODpath+'/' + fn0
   #registry_file0 =  AODpath+'/' + fn0
   ##print registry_fileName
  # 
   #if ( not reg_file_exist ) :
      #call(["cmsRun", "filenameRun.py", registry_file, registry_fileName])
   #else :
      #print registry_fileName
      #print "EXIST!! we will not overwrite it"
#
#
   ### run PFCandiate 
   #print "Running PFCandidate Producers now."
   #MODpath0  = destination_dir
   #call(["cmsRun", "PFCandidateRun.py", AODpath, MODpath0, registry_file0,  "1", "1"])
  # 
#
#
#
#
#
   #end = time()
   #print "Analyz file ",ii,  " done in " + str(end - start) + " seconds!"
#
#print "Analyz All files  done in " + str(end - start0) + " seconds!"
#
#
#
#
#
#
#
#
#






