from utilities.Analyzer_Onefile import *
import sys
import numpy as np
import json
#from skim_2010_Jet import *
#from save_jetpt_2010 import *

if __name__ =='__main__':

   list_of_downloads_file = sys.argv[1]
   destination_dir = sys.argv[2]
   fileNum =  int(sys.argv[3] )
   stepN = int(sys.argv[4] )
   post= sys.argv[6] 
   ifMC = sys.argv[7]
   GetRegisterFile = json.loads(sys.argv[8])

   rmAODMOD = True


   fileName = file_Names( list_of_downloads_file ) 
   if ( fileNum != -1 ) :
      fileName = fileName[:fileNum]

   fileName2 = []
   for i in range ( 0, len ( fileName) , stepN) :
      fileName2.append ( fileName[i:i+stepN])


   fileName= fileName2

   for ii, ff in enumerate ( fileName ) :
      start = time()
      if ii == 0  :
         start0 = start * 1 
      print ff



      AODpath, rootfilesName = download_root_files(ff, destination_dir)



      ## create registry
      print AODpath
      print rootfilesName
      registry_fileName = AODpath +'/' + rootfilesName[0][:-5] +'_to'+ str( stepN) + '_reg'\
            + post + '.txt'
      completed_fileName = AODpath +'/' + rootfilesName[0][:-5] +'_to'+ str( stepN) + '_completed'\
            + post + '.txt'
      print registry_fileName
      reg_file_exist= os.path.isfile ( registry_fileName)  and os.path.getsize( registry_fileName) > 0

      files_to_process = []
      for f1 in rootfilesName :
         root_file = "file://" + AODpath + '/'+f1
         files_to_process.append(root_file)

      #print files_to_process  
      files_to_process = sorted( files_to_process )
      #print files_to_process  
      files_to_process_file = AODpath +'/' + rootfilesName[0][:-5] +'_to'+ str( stepN) + \
            '_List'+ post + '.txt'
      print "files to process are ", files_to_process
      with open ( files_to_process_file, 'w') as f  :
         for s in files_to_process :
            f.write( "%s\n" % s)
      #np.savetxt(files_to_process_file , files_to_process , )
      
      
      #call(["rm",  registry_fileName ])
      #for root_file in files_to_process:              
         #call(["cmsRun", "filenameRun.py", root_file, registry_fileName ])

      #print "rootfile", files_to_process
      #break
      
      if ( GetRegisterFile)  :
         if ( reg_file_exist ) :
            print registry_fileName
            print "EXIST!! we will overwrite it"
            call(["rm",  registry_fileName ])

         for root_file in files_to_process:              
            call(["cmsRun", "filenameRun.py", root_file, registry_fileName , ifMC])

      #if ( ifMC == "MC") : 
      #   MODpath = AODpath + "/MOD/"
      #else :
      MODpath = AODpath.replace("AOD","MOD" )
      skimmedpath = AODpath.replace("AOD","skimmed" )
      print MODpath
      if not os.path.exists(MODpath):
         os.makedirs(MODpath)
      ## run PFCandiate 


      ### check MOD or plain ROOT file exist. If yes, delete them


      generateMODandROOT = True

      print "Running PFCandidate Producers now."


      if generateMODandROOT :
         for fn in rootfilesName : 
            fmod = fn[:-5]+".mod"
            call(["rm", MODpath + fmod])
            call(["rm", MODpath + fn])
         call(["cmsRun", "PFCandidateRun_onefile.py", AODpath, MODpath, files_to_process_file, 
               registry_fileName, "1", ifMC ])
               #registry_fileName, completed_fileName ,"1"])



      if ii < 1 :
         pass 
      elif rmAODMOD :
         for fn in rootfilesName : 
            fmod = fn[:-5]+".mod"
            call(["rm", MODpath + fmod])
            #call(["rm", MODpath + fn])
            call(["rm", AODpath + fn])
         
      end = time()
      print "Analyz file ",ii,  " done in " + str(end - start) + " seconds!"
      continue





