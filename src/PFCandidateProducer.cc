#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <limits>
#include <cmath>
#include <algorithm>
#include <sys/time.h>
#include <sys/stat.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"


#include "RecoJets/JetProducers/interface/BackgroundEstimator.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HLTrigger/HLTcore/interface/HLTConfigData.h"


// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"




#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// Muon
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Provenance/interface/Timestamp.h"

// ROOT headers.
#include "TLeaf.h"
#include "TTree.h"
#include "TFile.h"




using namespace std;
using namespace edm;
using namespace trigger;
using namespace reco;
using namespace fastjet;


class Writer {

public:

  // Methods.
  void    init (TTree *tree);
  void    add  (string branch);
  void    addint  (string branch);
  void    addstr  (string branch);
  int     size (string branch);
  double &var  (string branch, int idx = -2);
  int &varint  (string branch, int idx = -2);
  string &varstr  (string branch, int idx = -2);
  void    clear();

protected:

  // Members.
  double null;
  int nullint;
  string nullstr;
  vector< pair<string, vector<double> > > vars;
  vector< pair<string, vector<int> > > varsint;
  vector< pair<string, vector<string> > > varsstr;
};


//--------------------------------------------------------------------------

class PFCandidateProducer : public EDProducer 
{
public: 
   explicit PFCandidateProducer(const ParameterSet&);
   ~PFCandidateProducer();

private:
   virtual void beginJob();
   virtual void produce(edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;
   virtual void beginRun(edm::Run&, edm::EventSetup const&);
   virtual void endRun(edm::Run&, edm::EventSetup const&);
   virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
   virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
   virtual void writeraddList( );



   bool triggerFired(const string& triggerWildCard, const TriggerResults& triggerResults);
   unsigned int findTrigger(const string& triggerWildCard);


   bool file_exists(const std::string& name);

   
   InputTag srcCorrJets_;
   InputTag muonLabel_;


   HLTConfigProvider hltConfig_;
   InputTag hltInputTag_;
   InputTag rhoTag_;
   InputTag PFCandidateInputTag_;
   InputTag AK5PFInputTag_;
   
   InputTag lumiSummaryLabel_;
   
   int runNum;
   int eventNum;
   edm::LuminosityBlockNumber_t lumiBlockNumber_;
   
   long int startTime_;
   long int endTime_;

   int particleType;
   double px;
   double py;
   double pz;
   double energy;
   double mass;
   double area;
   
   long int eventSerialNumber_;
   
   FactorizedJetCorrector * AK5JetCorrector_;

   std::vector<std::string> filenames_;
   
   string mapFilename_;
   ifstream mapFile_;

   
   ifstream mapNumbersFile_;
   
   string outputDir_;
   ofstream fileOutput_;
   
   stringstream output_;
   
   string outputFilename_;
   string lastOutputFilename_;

   bool processFromTheBeginning_;
   
   string inputFile_;
   
   InputTag primaryVertices_;
   string dataVersion_;
   string dataSet_;



   string outputFilenameRoot_;
   string lastOutputFilenameRoot_;
   TFile* tFile_;
   TTree* tTree_;
   Writer writer;
   bool mIsMCarlo;
   bool skimFlag_;
};



PFCandidateProducer::PFCandidateProducer(const ParameterSet& iConfig)
: hltConfig_(),
  hltInputTag_("TriggerResults","","HLT"),
  rhoTag_(iConfig.getParameter<edm::InputTag>("rho")),
  PFCandidateInputTag_(iConfig.getParameter<InputTag>("PFCandidateInputTag")),
  AK5PFInputTag_(iConfig.getParameter<edm::InputTag>("AK5PFInputTag")),
  lumiSummaryLabel_(iConfig.getUntrackedParameter<edm::InputTag>("LumiSummaryLabel", InputTag("lumiProducer"))),
  primaryVertices_(iConfig.getParameter<InputTag>("primaryVertices")),
  dataVersion_(iConfig.getParameter<string>("dataVersion")),
  dataSet_(iConfig.getParameter<string>("dataSet")) 
{
  mapFilename_ = iConfig.getParameter<string>("mapFilename");
  mapFile_.open(mapFilename_.c_str()); 
  mIsMCarlo= iConfig.getUntrackedParameter<bool> ("isMCarlo",false);
  skimFlag_= iConfig.getUntrackedParameter<bool> ("skim",false);

  
  outputDir_ = iConfig.getParameter<string>("outputDir");
  outputFilename_ = "";
  lastOutputFilename_ = "";

  outputFilenameRoot_ = "";
  lastOutputFilenameRoot_ = "";

  cout << outputDir_ << endl;

  processFromTheBeginning_ = iConfig.getParameter<bool>("processFromTheBeginning");

  if ( ! processFromTheBeginning_)
    inputFile_ = iConfig.getParameter<string>("inputFile");
	  


   // muon
   //muonLabel_ = iConfig.getParameter<edm::InputTag>("muonLabel_");
}


PFCandidateProducer::~PFCandidateProducer() {

}

void PFCandidateProducer::produce(Event& iEvent, const EventSetup& iSetup) {

   string line, directory, filename;
   int fileEventNum, fileRunNum;

   getline(mapFile_, line);
   istringstream iss(line);
  	
   iss >> fileEventNum >> fileRunNum >> directory >> filename;

   runNum = iEvent.id().run();
   eventNum = iEvent.id().event();
   lumiBlockNumber_ = iEvent.luminosityBlock();
  // cout<< fileRunNum<< " "<< runNum<<endl;
   //cout<< fileEventNum<< " "<< eventNum<<endl;
   
   if ((fileRunNum == runNum) && (fileEventNum == eventNum)) {
   	
   	outputFilename_ = outputDir_ + "/" + filename.substr(0, filename.size() - strlen(".root")) + ".mod";
   	outputFilenameRoot_ = outputDir_ + "/" + filename;
      cout<< outputFilename_<< endl;
	
      if ((eventSerialNumber_ == 1) || (outputFilename_ != lastOutputFilename_)) {
         fileOutput_.close();
         fileOutput_.open(outputFilename_.c_str(), ios::out | ios::app );

         if ( lastOutputFilename_!= ""){
            cout<< "write TTree"<<endl;
            tFile_ -> Write( tTree_ -> GetName(), TObject::kOverwrite)  ;
            cout<< "close TFile"<< outputFilenameRoot_<<endl;
            tFile_ -> Close();
         }

         // should I use UPDATE or RECREAT
         cout<< "new TFile"<<endl;
         tFile_ = new TFile( outputFilenameRoot_.c_str(), "UPDATE")  ;
         tTree_ = new TTree ( "ntuple", "");
         writer.init(tTree_);
         

         lastOutputFilename_ = outputFilename_;
         lastOutputFilenameRoot_ = outputFilenameRoot_;
      }
   }
   
   writer.clear();
   output_.str("");
   output_.clear(); // Clear state flags.
   int dVersion;
   istringstream ( dataVersion_ ) >> dVersion;
   
   output_ << "BeginEvent Version " << dataVersion_ << " CMS_2010 Jet_Primary_Dataset" << endl;
   writer.varint( "version") = dVersion;
   


    Handle<vector<reco::Muon> > muons;
    iEvent.getByLabel(edm::InputTag("muons"), muons);



   
   
   // Primary Vertices.
   edm::Handle<VertexCollection> primaryVerticesHandle;
   iEvent.getByLabel( edm::InputTag("offlinePrimaryVertices"), primaryVerticesHandle);


        Vertex dummy;
        const Vertex *pv = &dummy;
        if (primaryVerticesHandle->size() != 0) {
           pv = &*primaryVerticesHandle->begin();
        } else { // create a dummy PV
          Vertex::Error e;
          e(0, 0) = 0.0015 * 0.0015;
          e(1, 1) = 0.0015 * 0.0015;
          e(2, 2) = 15. * 15.;
          Vertex::Point p(0, 0, 0);
          dummy = Vertex(p, e, 0, 0, 0);
        }

//   Handle<pat::MuonCollection> muonsHandle;
//   iEvent.getByLabel(edm::InputTag("selectedPatMuons"), muonsHandle);
//   const pat::MuonCollection & muons2 = *(muonsHandle.product());

 //   pat::MuonCollection::const_iterator amuon;
//    for (amuon = muons2.begin(); amuon != muons2.end(); amuon++) {
//       writer.varint("isGlobalMuon")  = amuon-> isGlobalMuon();
//       writer.varint("isPFMuon")  = amuon-> isPFMuon();
//       writer.varint("isTrackerMuon")  = amuon-> isPFMuon();
//       writer.var("Muon_px")  = amuon-> px();
//       writer.var("Muon_py")  = amuon-> py();
//       writer.var("Muon_pz")  = amuon-> pz();
//       writer.var("Muon_e")  = amuon-> energy();
//       cout<< amuon->px()<< " ,"<<  amuon->py()<< ", " ;
//       cout<< amuon->pz()<< " , " << amuon->energy()<<endl;
//       cout<< amuon->vertex().z()<< " , " << amuon-> charge()<<endl;
//       cout<< amuon->globalTrack()->hitPattern().numberOfValidMuonHits() <<endl;
//       cout<< amuon->globalTrack()->normalizedChi2()<<endl;
//       cout<< amuon->innerTrack()->hitPattern().numberOfValidPixelHits() <<endl;
//       cout<< amuon->innerTrack()->normalizedChi2()<<endl;
//       cout<< (*(amuon ->innerTrack())).dxy( pv-> position())<<endl;
//       cout<< pv->z() <<endl;
//       cout<<  amuon -> isolationR03().sumPt<<endl;

// }

   cout<< muons -> size()<<endl;
   writer.var("Pvertex_z")  = pv-> z();
    for (vector<reco::Muon>::const_iterator iMuon = muons->begin(); iMuon != muons->end(); ++iMuon) {
//        if (!iMuon->isGlobalMuon() || iMuon->pt() < 3.) continue;
         cout<< "Muon Pt "<<  iMuon -> pt()<<endl;
        //if (iMuon->pt() < 3.) continue;
        //if (!iMuon->isGlobalMuon()   && !iMuon->isPFMuon() && !iMuon->isTrackerMuon()  ) continue;
         writer.varint("isGlobalMuon")  = iMuon-> isGlobalMuon();
         
         writer.varint("isPFMuon")  = iMuon-> isPFMuon();
         writer.varint("isTrackerMuon")  = iMuon-> isTrackerMuon();
         writer.var("Muon_px")  = iMuon-> px();
         writer.var("Muon_py")  = iMuon-> py();
         writer.var("Muon_pz")  = iMuon-> pz();
         writer.var("Muon_e")  = iMuon-> energy();
         writer.varint("Muon_charge")  = iMuon-> charge();
         writer.var("Muon_vertex_z")  = iMuon-> vertex().z();
         //cout<<"Muon_vertex_z"<< iMuon-> vertex().z()<<endl;
         //cout<< (!iMuon->isGlobalMuon()   )<<endl; 
         //cout<< (!iMuon->isPFMuon()   )<<endl; 
         //cout<< (!iMuon->isTrackerMuon()   )<<endl; 
         //cout<< (!iMuon->isGlobalMuon()   && !iMuon->isPFMuon() && !iMuon->isTrackerMuon())<<endl; 
        if ( (!iMuon->isGlobalMuon())  )  {
         writer.var("Muon_global_chi2")  = -1 ;
         writer.varint("Muon_global_numhit")  = -1;
         }
         else {
         writer.var("Muon_global_chi2")  = iMuon->globalTrack()->normalizedChi2();
         writer.varint("Muon_global_numhit")  = iMuon->globalTrack()->hitPattern().numberOfValidMuonHits();

         }
         if (!iMuon->isGlobalMuon()   && !iMuon->isTrackerMuon()) {
            writer.var("Muon_inner_chi2")  = -1;
            writer.varint("Muon_inner_numhit")  =  -1 ;
            writer.var("Muon_ip")  = -1;
         }
         else {
            writer.var("Muon_inner_chi2")  = iMuon->innerTrack()->normalizedChi2();
            writer.varint("Muon_inner_numhit")  =  iMuon->innerTrack()->hitPattern().numberOfValidPixelHits();

            writer.var("Muon_ip")  = (*(iMuon ->innerTrack())).dxy( pv-> position());
         }
         //cout<< "end mu "<<endl;
         writer.var("Muon_isoR03_sumPt")  = iMuon -> isolationR03().sumPt;
         writer.var("Muon_isoR03_emEt")  = iMuon -> isolationR03().emEt;
         writer.var("Muon_isoR03_hadEt")  = iMuon -> isolationR03().hadEt;

   }



   // Muon
    //edm::Handle<vector<pat::Muon> > h_mu;
   //iEvent.getByLabel( edm::InputTag("offlinePrimaryVertices"), primaryVerticesHandle);
   
   
   // Luminosity Block Begins
   
   
   LuminosityBlock const& iLumi = iEvent.getLuminosityBlock();
   Handle<LumiSummary> lumi;
   iLumi.getByLabel(lumiSummaryLabel_, lumi);
      
   // Luminosity Block Ends
   if ( mIsMCarlo ) {
      output_ << "#   Cond          RunNum        EventNum     " << endl;
      output_ << "    Cond"
       << setw(16) << runNum
       << setw(16) << eventNum
	       << setw(16) << lumiBlockNumber_
 	       << endl;   
   // timeValue = timeHigh (unixTime); timeValue = timeValue << 32; timeValue += microsecondsOffset;
   writer.varint( "runNum") = runNum;
   writer.varint( "eventNum") = eventNum;
   writer.varint( "lumiBlock") = lumiBlockNumber_;
   }
   else 
   {
   output_ << "#   Cond          RunNum        EventNum       LumiBlock       validLumi     intgDelLumi     intgRecLumi     AvgInstLumi             NPV       timestamp        msOffset" << endl;
   output_ << "    Cond"
   	       << setw(16) << runNum
	       << setw(16) << eventNum
	       << setw(16) << lumiBlockNumber_
   	       << setw(16) << lumi->isValid()
   	       << setw(16) << lumi->intgDelLumi()
   	       << setw(16) << lumi->intgRecLumi()
   	       << setw(16) << lumi->avgInsDelLumi()
   	       << setw(16) << primaryVerticesHandle->size()
	       << setw(16) << iEvent.time().unixTime()
	       << setw(16) << iEvent.time().microsecondOffset()
 	       << endl;   
   // timeValue = timeHigh (unixTime); timeValue = timeValue << 32; timeValue += microsecondsOffset;
   writer.varint( "runNum") = runNum;
   writer.varint( "eventNum") = eventNum;
   writer.varint( "lumiBlock") = lumiBlockNumber_;
   writer.varint("validLumi") = lumi -> isValid();
   writer.var("intgDelLumi") = lumi -> intgDelLumi();
   writer.var("intgRecLumi") = lumi -> intgRecLumi();
   writer.var("avgInstLumi") = lumi -> avgInsDelLumi();
   writer.varint("NPV") = primaryVerticesHandle->size();
   writer.varint("timestamp") = iEvent.time().unixTime();
   writer.varint("msOffset") = iEvent.time().microsecondOffset();
   }

   



//566    writer.addint("msOffset");


   

   Handle<reco::PFCandidateCollection> PFCollection;
   iEvent.getByLabel(PFCandidateInputTag_, PFCollection);
   
   Handle<TriggerResults> trigResults; 
   iEvent.getByLabel(hltInputTag_, trigResults);
   
   edm::Handle<reco::PFJetCollection> AK5Collection;
   iEvent.getByLabel(AK5PFInputTag_, AK5Collection);
   
   if ( ! PFCollection.isValid()){
    cerr << "Invalid collection." << endl;
    return;
   }
  //if ( ! AK5Collection.isValid()){
   //std::cerr << "Invalid collection." << std::endl;
   //return;
   //}

   
   
   
   vector<PseudoJet> PFCForFastJet;
   
   double rapmin = std::numeric_limits<double>::min();
   double rapmax = std::numeric_limits<double>::max();
   for(reco::PFCandidateCollection::const_iterator it = PFCollection->begin(), end = PFCollection->end(); it != end; it++) {
      PFCForFastJet.push_back(PseudoJet(it->px(), it->py(), it->pz(), it->energy()));
      
      if (it->rapidity() < rapmin)
      	rapmin = it->rapidity();
      if (it->rapidity() > rapmax)
      	rapmax = it->rapidity();
   }
   
   
   // Record trigger information first.
   
   
   // Get all trigger names associated with the "Jet" dataset.
   const vector<string> triggerNames = hltConfig_.datasetContent(dataSet_.c_str());
   
   for (unsigned i = 0; i < triggerNames.size(); i++) {
      if (i == 0)
         output_ << "#   Trig                            Name      Prescale_1      Prescale_2          Fired?" << endl;
      
      string name = triggerNames[i];
      
      pair<int, int> prescale = hltConfig_.prescaleValues(iEvent, iSetup, name);

      bool fired = triggerFired(name, ( * trigResults));

      output_ << "    Trig"
       	          << setw(32) << name
	          << setw(16) << prescale.first
	          << setw(16) << prescale.second
                  << setw(16) << fired
                  << endl;
      writer.varstr( "trig") = name;
      writer.var( "pre1") = prescale.first;
      writer.var( "pre2") = prescale.second;
      writer.varint("fired") =  fired;
     // else  writer.addint("fired") =  0;
   
   
   }
   
  // Get AK5 Jets.
  
  // Setup background density for AK5 JEC.
  
  edm::Handle<double> rhoHandle;
  iEvent.getByLabel( edm::InputTag("kt6PFJetsForIsolation", "rho"), rhoHandle);
  double rho = * rhoHandle;
  
  for(reco::PFJetCollection::const_iterator it = AK5Collection->begin(), end = AK5Collection->end(); it != end; it++) {    
    if (it == AK5Collection->begin())
       output_ << "#    AK5" << "              px              py              pz          energy             jec            area     no_of_const     chrg_multip    neu_had_frac     neu_em_frac   chrg_had_frac    chrg_em_frac" << endl;
    
    px = it->px();
    py = it->py();
    pz = it->pz();
    energy = it->energy();
    area = it->jetArea();
    
    // JEC Factor.
    
    AK5JetCorrector_->setJetEta(it->eta());
    AK5JetCorrector_->setJetPt(it->pt());
    AK5JetCorrector_->setJetA(it->jetArea());
    AK5JetCorrector_->setRho(rho);
         
    double correction = AK5JetCorrector_->getCorrection();

    // Jet Quality Cut Parameters.

    double neutral_hadron_fraction = it->neutralHadronEnergy() / it->energy();
    double neutral_em_fraction = it->neutralEmEnergy() / it->energy();
    int number_of_constituents = it->nConstituents();
    double charged_hadron_fraction = it->chargedHadronEnergy() / it->energy();
    int charged_multiplicity = it->chargedMultiplicity();
    double charged_em_fraction = it->chargedEmEnergy() / it->energy();
 
    output_ << "     AK5"
        << setw(16) << fixed << setprecision(8) << px
        << setw(16) << fixed << setprecision(8) << py
        << setw(16) << fixed << setprecision(8) << pz
        << setw(16) << fixed << setprecision(8) << energy
        << setw(16) << fixed << setprecision(8) << correction
        << setw(16) << fixed << setprecision(8) << area   
        << setw(16) << fixed << setprecision(8) << number_of_constituents   
        << setw(16) << fixed << setprecision(8) << charged_multiplicity  
        << setw(16) << fixed << setprecision(8) << neutral_hadron_fraction   
        << setw(16) << fixed << setprecision(8) << neutral_em_fraction   
        << setw(16) << fixed << setprecision(8) << charged_hadron_fraction    
        << setw(16) << fixed << setprecision(8) << charged_em_fraction       
        << endl;

   

   writer.var("AK5_px") = px;
   writer.var("AK5_py") = py;
   writer.var("AK5_pz") = pz;
   writer.var("AK5_e")  = energy;
   writer.var("AK5_jec") = correction;
   writer.var("AK5_area") = area;
   writer.varint("AK5_num_const") = number_of_constituents;
   writer.varint("AK5_charg_multip") = charged_multiplicity;
   writer.var("AK5_neu_had_frac") = neutral_hadron_fraction;
   writer.var("AK5_neu_em_frac") = neutral_em_fraction;
   writer.var("AK5_chrg_had_frac") = charged_hadron_fraction;
   writer.var("AK5_chrg_em_frac") = charged_em_fraction;

  }
  
  
  // Get PFCandidates.
   if ( !skimFlag_) {
  for(reco::PFCandidateCollection::const_iterator it = PFCollection->begin(), end = PFCollection->end(); it != end; it++) {
    if (it == PFCollection->begin())
       output_ << "#    PFC" << "              px              py              pz          energy           pdgId" << endl;  
    
    px = it->px();
    py = it->py();
    pz = it->pz();
    energy = it->energy();
    int pdgId = it->pdgId();
    
    output_ << "     PFC"
        << setw(16) << fixed << setprecision(8) << px
        << setw(16) << fixed << setprecision(8) << py
        << setw(16) << fixed << setprecision(8) << pz
        << setw(16) << fixed << setprecision(8) << energy
        << setw(16) << noshowpos << pdgId
        << endl;
   writer.var("PFC_px") = px;
   writer.var("PFC_py") = py;
   writer.var("PFC_pz") = pz;
   writer.var("PFC_e") = energy;
   writer.varint("PFC_pid") = pdgId;
   }
   }

   if ( mIsMCarlo ) {

     edm::Handle<GenParticleCollection> genParticles;
     iEvent.getByLabel("genParticles", genParticles);
     for(unsigned int i = 0; i < genParticles->size(); ++ i) {
        const GenParticle & p = (*genParticles)[i];
        int id;
         id = p.pdgId();
        int st ;
        st = p.status();  
        double pt = p.pt();
         //double eta = p.eta(), phi = p.phi(), mass = p.mass();
        double vx = p.vx(), vy = p.vy(), vz = p.vz();
        int charge = p.charge();
        if ( abs( id ) == 13 && st == 1 ) {
            writer.var("GenMuon_px") = p.px();
            writer.var("GenMuon_py") = p.py();
            writer.var("GenMuon_pz") = p.pz();
            writer.var("GenMuon_x") = vx;
            writer.var("GenMuon_y") = vy;
            writer.var("GenMuon_z") = vz;
            writer.var("GenMuon_e") = p.energy();
            writer.varint("GenMuon_charge") = charge;
            writer.varint("GenMuon_pid") = id;
            cout<< "Gen Muon pt "<<  pt<<"  "<< st<< endl;
         }
      }
   }
   
   
   output_ << "EndEvent" << endl;
   
   fileOutput_ << output_.rdbuf();
   
   eventSerialNumber_++;
   tTree_  -> Fill();
}

void PFCandidateProducer::beginJob() {
   eventSerialNumber_ = 1;
   
   // Start timer.
   struct timeval tp;
   gettimeofday(&tp, NULL);
   startTime_ = tp.tv_sec * 1000 + tp.tv_usec / 1000;
   writeraddList( );
   
   // Figure out the JetCorrector objects for AK5 corrections.
   
   
   // AK5
   
   // Create the JetCorrectorParameter objects, the order does not matter.
   // YYYY is the first part of the txt files: usually the global tag from which they are retrieved
   JetCorrectorParameters *AK5ResJetPar = new JetCorrectorParameters("data/JEC/GR_R_42_V25_AK5PF_L2L3Residual.txt"); 
   JetCorrectorParameters *AK5L3JetPar  = new JetCorrectorParameters("data/JEC/GR_R_42_V25_AK5PF_L3Absolute.txt");
   JetCorrectorParameters *AK5L2JetPar  = new JetCorrectorParameters("data/JEC/GR_R_42_V25_AK5PF_L2Relative.txt");
   JetCorrectorParameters *AK5L1JetPar  = new JetCorrectorParameters("data/JEC/GR_R_42_V25_AK5PF_L1FastJet.txt");
   
   //  Load the JetCorrectorParameter objects into a vector, IMPORTANT: THE ORDER MATTERS HERE !!!! 
   vector<JetCorrectorParameters> vParAK5;
   vParAK5.push_back(*AK5L1JetPar);
   vParAK5.push_back(*AK5L2JetPar);
   vParAK5.push_back(*AK5L3JetPar);
   vParAK5.push_back(*AK5ResJetPar);
   
   AK5JetCorrector_ = new FactorizedJetCorrector(vParAK5);

   
   std::cout << "Processing PFCandidates." << std::endl;
   
   // Map thingy.

   if ( ! processFromTheBeginning_) {
   	
   	string line, directory;
   	int fileEventNum, fileRunNum;
   	int linesDown = 1;
   
	ifstream registryFile(mapFilename_.c_str());
   
   	string rootFilename = "";
   	
   	while((rootFilename != inputFile_)) {
   		
		getline(registryFile, line);
		istringstream iss(line);
   		iss >> fileEventNum >> fileRunNum >> directory >> rootFilename;
   		linesDown++;
	}
	
	cout << "Trying to find the correct line here!" << endl;
	cout << linesDown << endl;
	
	for(int i = 0; i < linesDown - 2; i++) {
		getline(mapFile_, line);
		istringstream iss(line);
		iss >> fileEventNum >> fileRunNum >> directory >> rootFilename;
	}
	
   }


}

void PFCandidateProducer::endJob() {
   struct timeval tp2;
   gettimeofday(&tp2, NULL);
   long int endTime_ = tp2.tv_sec * 1000 + tp2.tv_usec / 1000;   
   double elapsed_milliseconds = endTime_ - startTime_;


   cout<< "write TTree"<<endl;
   tFile_ -> Write( tTree_ -> GetName(), TObject::kOverwrite)  ;
   cout<< "close TFile"<< outputFilenameRoot_<<endl;
   tFile_ -> Close();
   
   cout << endl << endl << endl << "Finished processing " << (eventSerialNumber_ - 1) << " events in " << elapsed_milliseconds / (60*1000) << " minutes!" << endl;
}

void PFCandidateProducer::beginRun(edm::Run & iRun, edm::EventSetup const & iSetup){

   bool changed = true;
   if ( hltConfig_.init(iRun, iSetup, hltInputTag_.process(), changed) ) {
      // if init returns TRUE, initialisation has succeeded!
      edm::LogInfo("TopPairElectronPlusJetsSelectionFilter") << "HLT config with process name "
        << hltInputTag_.process() << " successfully extracted";
   }
   else {
      edm::LogError("TopPairElectronPlusJetsSelectionFilter_Error")
      << "Error! HLT config extraction with process name " << hltInputTag_.process() << " failed";
   }


  cout << "This is beginJob()" << endl;

}

void PFCandidateProducer::endRun(edm::Run&, edm::EventSetup const&) {

}

void PFCandidateProducer::beginLuminosityBlock(edm::LuminosityBlock& iLumi, edm::EventSetup const& iSetup) {
   
}

void PFCandidateProducer::endLuminosityBlock(edm::LuminosityBlock& iLumi, edm::EventSetup const& iSetup) {

}

bool PFCandidateProducer::triggerFired(const std::string& triggerWildCard, const edm::TriggerResults& triggerResults) {
   bool fired = false;
   unsigned int index = findTrigger(triggerWildCard);

   if (index < triggerResults.size()) {
      if (triggerResults.accept(index)) {
         fired = true;
      }
   }

   return fired;

}


void PFCandidateProducer::writeraddList( ) {
   writer.addint("version");
   writer.addint("eventNum");
   writer.addint("runNum");
   writer.addint("lumiBlock");
   writer.addint("validLumi");
   writer.add("intgDelLumi");
   writer.add("intgRecLumi");
   writer.add("avgInstLumi");
   writer.addint("NPV");
   writer.addint("timestamp");
   writer.addint("msOffset");

   // trigger
   writer.addstr("trig");
   writer.add("pre1");
   writer.add("pre2");
   writer.addint("fired");

   //jet 
   writer.add("AK5_px");
   writer.add("AK5_py");
   writer.add("AK5_pz");
   writer.add("AK5_e");
   writer.add("AK5_jec");
   writer.add("AK5_area");
   writer.addint("AK5_num_const");
   writer.addint("AK5_charg_multip");
   writer.add("AK5_neu_had_frac");
   writer.add("AK5_neu_em_frac");
   writer.add("AK5_chrg_had_frac");
   writer.add("AK5_chrg_em_frac");

   writer.add("PFC_px");
   writer.add("PFC_py");
   writer.add("PFC_pz");
   writer.add("PFC_e");
   writer.addint("PFC_pid");

   // Muon
   writer.addint("isGlobalMuon");
   writer.addint("isPFMuon");
   writer.addint("isTrackerMuon");
   writer.add("Muon_px");
   writer.add("Muon_py");
   writer.add("Muon_pz");
   writer.add("Muon_e");
   writer.addint("Muon_charge");

   writer.add("Muon_vertex_z");
   writer.add("Pvertex_z");

   writer.add("Muon_global_chi2");
   writer.addint("Muon_global_numhit");
   writer.addint("Muon_inner_numhit");
   writer.add("Muon_inner_chi2");

   writer.add("Muon_ip");

   writer.add("Muon_isoR03_sumPt");
   writer.add("Muon_isoR03_emEt");
   writer.add("Muon_isoR03_hadEt");

   if (mIsMCarlo){
      writer.add("GenMuon_px");
      writer.add("GenMuon_py");
      writer.add("GenMuon_pz");
      writer.add("GenMuon_x");
      writer.add("GenMuon_y");
      writer.add("GenMuon_z");
      writer.add("GenMuon_e");
      writer.addint("GenMuon_charge");
      writer.addint("GenMuon_pid");
   }


}

unsigned int PFCandidateProducer::findTrigger(const std::string& triggerWildCard) {
   const std::vector<std::string>& triggers = hltConfig_.triggerNames();
   unsigned int found = 9999;

   size_t length = triggerWildCard.size();
   for (unsigned int index = 0; index < triggers.size(); ++index) {
      if (length <= triggers[index].size() && triggerWildCard == triggers[index].substr(0, length)) {
         found = index;
         break;
      }
   }

   return found;
}

bool PFCandidateProducer::file_exists(const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}





DEFINE_FWK_MODULE(PFCandidateProducer);




//==========================================================================

// Reader class to access the variables from an input TTree with
// standard branches. Currently the following branch types are
// supported (see https://root.cern.ch/doc/master/classTTree.html for
// abbreviations):
//
//   I:  Int_t     
//   F:  Float_t 
//   D:  Double_t 
//   l:  ULong64_t
//   O:  Bool_t
//   VF: vector<Float_t>
//   VD: vector<Double_t>
//
// The init(TTree *tree) method sets the branch addresses for the tree
// to the internal variables.

// Writer class to produce a TTree ntuple of vector<Double_t>.

//--------------------------------------------------------------------------

// Initialize a tree.

void Writer::init(TTree *tree) {
  for (int i = 0; i < (int)vars.size(); i++)
    tree->Branch(vars[i].first.c_str(), &vars[i].second);
  for (int i = 0; i < (int)varsint.size(); i++)
    tree->Branch(varsint[i].first.c_str(), &varsint[i].second);
  for (int i = 0; i < (int)varsstr.size(); i++)
    tree->Branch(varsstr[i].first.c_str(), &varsstr[i].second);
}

//--------------------------------------------------------------------------

// Add a branch.

void Writer::add(string branch) {
  for (int i = 0; i < (int)vars.size(); i++) 
    if (vars[i].first == branch) return;
  vars.push_back(pair<string, vector<double> >(branch, vector<double>()));
}

void Writer::addint(string branch) {
  for (int i = 0; i < (int)varsint.size(); i++) 
    if (varsint[i].first == branch) return;
  varsint.push_back(pair<string, vector<int> >(branch, vector<int>()));
}

void Writer::addstr(string branch) {
  for (int i = 0; i < (int)varsstr.size(); i++) 
    if (varsstr[i].first == branch) return;
  varsstr.push_back(pair<string, vector<string> >(branch, vector<string>()));
}



//--------------------------------------------------------------------------

// Access a variable (-2 appends and returns, -1 is the last entry).

double& Writer::var(string branch, int idx) {
  for (unsigned int i = 0; i < vars.size(); i++) {
    if (vars[i].first != branch) continue;
    if (idx < 0) {
      if (idx < -1) vars[i].second.push_back(0);
      return vars[i].second.back();
    }
    else if (idx < (int)vars[i].second.size())
      return vars[i].second[idx];
  }
  return null;
}


int& Writer::varint(string branch, int idx) {
  for (unsigned int i = 0; i < varsint.size(); i++) {
    if (varsint[i].first != branch) continue;
    if (idx < 0) {
      if (idx < -1) varsint[i].second.push_back(0);
      return varsint[i].second.back();
    }
    else if (idx < (int)varsint[i].second.size())
      return varsint[i].second[idx];
  }
  return nullint;
}

string& Writer::varstr(string branch, int idx) {
  for (unsigned int i = 0; i < varsstr.size(); i++) {
    if (varsstr[i].first != branch) continue;
    if (idx < 0) {
      if (idx < -1) varsstr[i].second.push_back("");
      return varsstr[i].second.back();
    }
    else if (idx < (int)varsstr[i].second.size())
      return varsstr[i].second[idx];
  }
  return nullstr;
}
//--------------------------------------------------------------------------

// Determine the size of a branch vector.

int Writer::size(string branch) {
  for (unsigned int i = 0; i < vars.size(); i++)
    if (vars[i].first == branch) return vars[i].second.size();
  for (unsigned int i = 0; i < varsint.size(); i++)
    if (varsint[i].first == branch) return varsint[i].second.size();
  for (unsigned int i = 0; i < varsstr.size(); i++)
    if (varsstr[i].first == branch) return varsstr[i].second.size();
  return -1;
}

//--------------------------------------------------------------------------

// Clear the branches.
void Writer::clear() {
  for (int i = 0; i < (int)vars.size(); i++) vars[i].second.clear();
  for (int i = 0; i < (int)varsint.size(); i++) varsint[i].second.clear();
  for (int i = 0; i < (int)varsstr.size(); i++) varsstr[i].second.clear();
   //vars.clear();
  // varsint.clear();
}

//==========================================================================
