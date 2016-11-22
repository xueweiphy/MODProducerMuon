# MODProducerMuon

This package downloads AOD files from the 
[CERN Open Data Portal](http://opendata.cern.ch "CERN Open Data Portal")
release and 
converts them into a human-readable file format called MOD root(MIT Open Data). 
Currently, the following information are stored:

- 4-momentum, Jet Area and Jet Energy Correction factors for AK5 Jets.
- Trigger Names, prescale values and whether or not that trigger fired for each event.
- Run Number and Event Number.
- Luminosity Block, Average Instantaneous Luminosity and Number of Primary Vertices.
- Muon information

## Usage Instruction

-  Create the working area:

   ```
   cmsrel CMSSW_5_3_32

   cd ./CMSSW_5_3_32/src

   cmsenv
   
   git clone https://github.com/xueweiphy/MODProducerMuon.git
   
   cd MODProducerMuon

   scram b
   ```

-  Create links to the condition databases

   ```
   ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA FT_53_LV5_AN1 
   ln -sf /cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db FT_53_LV5_AN1_RUNA.db
   ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1 START53_LV6A1
   ln -sf /cvmfs/cms-opendata-conddb.cern.ch/START53_LV6A1.db START53_LV6A1.db
   ```

- Run 

   ```
   ./run_producers_2011_notsave.sh
   ```

### Workflow

1.  Download some of the ROOT files and arrange them in the same directory structure
   as they live on the CMS servers

2. Create a registry that maps each event and run number to a certain ROOT file.

3. Run the Producer on those AOD files. This reads the download directory and 
   processes only the files in there. This produces N MOD files.










