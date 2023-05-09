#define _USE_MATH_DEFINES
#include <cmath> 
#include <iostream>
using namespace std;

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include "Math/VectorUtil.h"
#include <stdlib.h>


#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

//Vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

//Electron ID stuff
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

//Photon ID stuff
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//Parton distribution and QCD scale variations stuff 
#include "FWCore/Framework/interface/Run.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" //LHE reader
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" //LHE reader

//#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" //JEC uncertainties
//#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h" //JEC uncertainties
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"

//JEC uncertainties
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
//#include "CondFormats/JetMETObjects/interface/JetResolution.h"
//#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//#include "CondFormats/DataRecord/interface/JetCorrectionsRecord.h"

//#include "CondFormats/JetMETObjects/interface/JetCorrectionsRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

typedef math::XYZTLorentzVector LorentzVector;
 
#include "ZMesonGamma.h"

 
// constructors and destructor
ZMesonGamma::ZMesonGamma(const edm::ParameterSet& iConfig) : 
runningOnData_(iConfig.getParameter<bool>("runningOnData")),
verboseIdFlag_(iConfig.getParameter<bool>("phoIdVerbose"))//,
//effectiveAreas_el_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_el")).fullPath() ),
//effectiveAreas_ph_( (iConfig.getParameter<edm::FileInPath>("effAreasConfigFile_ph")).fullPath() )
{
  prunedGenParticlesToken_            = consumes<std::vector<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  offlineSlimmedPrimaryVerticesToken_ = consumes<std::vector<reco::Vertex> > (edm::InputTag("offlineSlimmedPrimaryVertices"));
  pileupSummaryToken_                 = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  triggerBitsToken_                   = consumes<edm::TriggerResults> (edm::InputTag("TriggerResults","","HLT"));
  photonsMiniAODToken_                = consumes<std::vector<pat::Photon> > (edm::InputTag("slimmedPhotons"));
  slimmedMuonsToken_                  = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  electronsMiniAODToken_              = consumes<std::vector<pat::Electron> > (edm::InputTag("slimmedElectrons"));
  slimmedJetsToken_                   = consumes<std::vector<pat::Jet> >(edm::InputTag("slimmedJets"));
  GenInfoToken_                       = consumes<GenEventInfoProduct> (edm::InputTag("generator"));
  slimmedMETsToken_                   = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));
  slimmedMETsPuppiToken_              = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETsPuppi"));
  rhoToken_                           = consumes<double> (iConfig.getParameter <edm::InputTag>("rho"));



  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 8, 0., 8.);

  _Nevents_processed  = 0;
  _Nevents_triggered  = 0;
  _Nevents_isPhoton   = 0;
  _Nevents_isTwoKaons = 0;


  debug=false;  //DEBUG datamember 
  verbose=false; 

  h_pileup   = fs->make<TH1F>("pileup", "pileup", 75,0,75);

  create_trees();
}

ZMesonGamma::~ZMesonGamma()
{
}

// ------------ method called for each event  ------------
void ZMesonGamma::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<std::vector<reco::GenParticle>  > genParticles;
  iEvent.getByToken(prunedGenParticlesToken_, genParticles);

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBitsToken_, triggerBits);  

  edm::Handle<std::vector<pat::Photon> > slimmedPhotons;
  iEvent.getByToken(photonsMiniAODToken_,slimmedPhotons);

  edm::Handle<std::vector<pat::Muon>  > slimmedMuons;
  iEvent.getByToken(slimmedMuonsToken_, slimmedMuons);

  edm::Handle<std::vector<pat::Electron> > slimmedElectrons;
  iEvent.getByToken(electronsMiniAODToken_,slimmedElectrons);

  edm::Handle<std::vector<pat::Jet > > slimmedJets;
  iEvent.getByToken(slimmedJetsToken_, slimmedJets);

  edm::Handle<std::vector<reco::Vertex > > slimmedPV;
  iEvent.getByToken(offlineSlimmedPrimaryVerticesToken_, slimmedPV);

  edm::Handle<std::vector<pat::MET > > slimmedMETs;
  iEvent.getByToken(slimmedMETsToken_, slimmedMETs);

  edm::Handle<std::vector<pat::MET > > slimmedMETsPuppi;
  iEvent.getByToken(slimmedMETsPuppiToken_, slimmedMETsPuppi);


   _Nevents_processed++; //This will be saved in the output tree, giving the number of processed events
   
  //Retrieve the run number
  if(runningOnData_){
    run_number = iEvent.id().run();
    event_number = iEvent.id().event();
  }


  //*************************************************************//
  //                                                             //
  //-------------------------- Vertices -------------------------//
  //                                                             //
  //*************************************************************//

  //Count the number of vertices and return if there's no vertex
  nPV = 0;
  if(slimmedPV->size()<=0){
    if(verbose) cout<<"No primary vertex found, RETURN."<<endl;
    return;
  }

  for(reco::VertexCollection::const_iterator vtx=slimmedPV->begin();vtx!=slimmedPV->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      nPV++;
    }
  } 


  //*************************************************************//
  //                                                             //
  //--------------------------- Pile Up -------------------------//
  //                                                             //
  //*************************************************************//

  PU_Weight = -1.;
  float npT = -1.;

  if(!runningOnData_){
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI; 

    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      const int BX = PVI->getBunchCrossing();
      if(BX == 0) {
        npT  = PVI->getTrueNumInteractions();
      }
    }

    if(npT == -1) {
      std::cout << "!!!! npT = -1 !!!!" << std::endl;
      abort();
    }

  // Calculate weight using above code
      PU_Weight = Lumiweights_.weight(npT);

      // Fill histogram with PU distribution
      h_pileup->Fill(npT);
    }


  //*************************************************************//
  //                                                             //
  //-------------------------- MC Weight ------------------------//
  //                                                             //
  //*************************************************************//

  MC_Weight = -10000000.;

  if(!runningOnData_){
    edm::Handle<GenEventInfoProduct> GenInfo;
    iEvent.getByToken(GenInfoToken_, GenInfo);
    
    float _aMCatNLOweight = GenInfo->weight();
    MC_Weight = _aMCatNLOweight;

    if(MC_Weight == -10000000.) {
      std::cout << "!!!! MC_Weight = -10000000 !!!!" << std::endl;
      abort();
    }
  }

  //*************************************************************//
  //                                                             //
  //--------------------------- Trigger -------------------------//
  //                                                             //
  //*************************************************************//

  //Examine the trigger information, return if the trigger doesn't switch on and count the number of events where the trigger has switched on
  isTwoProngTrigger = false;

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  for(unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if(!triggerBits->accept(i)) continue;
    std::string tmp_triggername = names.triggerName(i);

    if( tmp_triggername.find("HLT_Photon35_TwoProngs35_v") != std::string::npos ){
      isTwoProngTrigger = true;
    }
  }

  if(!isTwoProngTrigger){
   if(verbose) cout<<"Event not triggered, RETURN."<<endl;
   return;
  }
  _Nevents_triggered++;


  //*************************************************************//
  //                                                             //
  //------------------ Variable initialization ------------------//
  //                                                             //
  //*************************************************************//

 
   nPhotons38WP80        = 0;
   nPhotons20WP90        = 0;
   nPhotonsChosen        = 0;
   
   
    //These variables will go in the tree
   ph_eT     = 0.;
   ph_eta    = 0.;
   ph_etaSC  = 0.;
   ph_phi    = 0.;
   LorentzVector ph_p4;

   eTphMax = -1000.;

  _Jet_Photon_invMass  = -1.;
  _MesonMass           = -1.;

  met_pT      = 0.;
  metpuppi_pT = 0.;

  _bestJet_pT                       =-1.;
  _bestJet_eta                      =-1.;
  _bestJet_phi                      =-1.;
  _bestJet_nDaughters               = 0.;
  _bestJet_pTMax                    =-1.;
  _bestJet_chargedEmEnergy          = 0.;
  _bestJet_neutralEmEnergy          = 0.;
  _bestJet_chargedHadEnergy         = 0.;
  _bestJet_neutralHadEnergy         = 0.;
  _bestJet_chargedEmEnergyFraction  = 0.;
  _bestJet_neutralEmEnergyFraction  = 0.;
  _bestJet_chargedHadEnergyFraction = 0.;
  _bestJet_neutralHadEnergyFraction = 0.;
  _bestJet_chargedHadMultiplicity   = 0.;
  _bestJet_invMass                  = 0.;
  _bestJet_Photon_invMass           = 0.;
  _bestJet_JECunc                   = 0.;
  _firstCandPt                      = 0.;
  _firstCandEta                     = 0.;
  _firstCandPhi                     = 0.;
  _firstCandCharge                  = 0.;
  _secondCandPt                     = 0.;
  _secondCandEta                    = 0.;
  _secondCandPhi                    = 0.;
  _secondCandCharge                 = 0.;
  _bestCouplePt                     = 0.;
  _bestCoupleEta                    = 0.;
  _bestCouplePhi                    = 0.;  


  //*************************************************************//
  //                                                             //
  //----------------------------- MET ---------------------------//
  //                                                             //
  //*************************************************************//
 for(auto met = slimmedMETs->begin(); met != slimmedMETs->end(); ++met){
  met_pT = met->pt();
}

for(auto metpuppi = slimmedMETsPuppi->begin(); metpuppi != slimmedMETsPuppi->end(); ++metpuppi){
  metpuppi_pT = metpuppi->pt();
}



  //*************************************************************//
  //                                                             //
  //---------------------------- Muons --------------------------//
  //                                                             //
  //*************************************************************//
 //Count muons for each event
for(auto mu = slimmedMuons->begin(); mu != slimmedMuons->end(); ++mu){
  if(mu->pt() < 10. || !mu->CutBasedIdMedium || fabs(mu->eta()) > 2.4 || fabs(mu->muonBestTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(mu->muonBestTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;
  if(!mu->PFIsoLoose) continue;
  nMuons10++;
  if(mu->pt() < 20.) continue;
  nMuons20++;
}



//*************************************************************//
  //                                                             //
  //-------------------------- Electrons ------------------------//
  //                                                             //
  //*************************************************************//
  //Count the number of electrons
  // Get rho value
edm::Handle< double > rhoH;
iEvent.getByToken(rhoToken_,rhoH);
rho_ = *rhoH;

float corr_pt = 0.;
//ph_en_sigmaUP = 0.;
//ph_en_sigmaDW = 0.;
//ph_en_scaleUP = 0.;
//ph_en_scaleDW = 0.;

for(auto el = slimmedElectrons->begin(); el != slimmedElectrons->end(); ++el){
    //Calculate electron p4, correct it with the Scale&Smearing correction and extract the pT
    LorentzVector el_p4 = el->p4() * el->userFloat("ecalTrkEnergyPostCorr")/el->energy();
    corr_pt = el_p4.pt();

    if(corr_pt < 10. || fabs(el->eta()) > 2.5 || fabs(el->gsfTrack()->dxy((&slimmedPV->at(0))->position())) >= 0.2 || fabs(el->gsfTrack()->dz((&slimmedPV->at(0))->position())) >= 0.5) continue;

     //-------------Conditions on loose/medium MVA electron ID-------------//
    if(el->electronID("mvaEleID-Fall17-iso-V2-wp80") == 0) continue;
    nElectrons10++;
    if (corr_pt < 20.) continue;
    nElectrons20++;
  }


  //*************************************************************//
  //                                                             //
  //--------------------------- Photons -------------------------//
  //                                                             //
  //*************************************************************//
  if(verbose) cout<< "PHOTONs"<<" --------------------------------"<<endl;
  
  bool cand_photon_found = false; //initialize this bool to false, return if it doesn't turn into true
  float corr_et = -1.;

  for(auto photon = slimmedPhotons->begin(); photon != slimmedPhotons->end(); ++photon){ //PHOTON FORLOOP START --------------------------------
    
    // Apply energy scale corrections, from 18Apr2023 the correction are embedded in the config file with the postReco tool
    corr_et   = photon->et(); //* photon->userFloat("ecalEnergyPostCorr") / photon->energy(); 

    if(corr_et < 20. || fabs(photon->eta()) > 2.5) continue; //loose selection to reject diphoton bkg 
    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp90") == 0) continue; //WP90
    if(!photon->passElectronVeto()) continue; 

    nPhotons20WP90++;

    if(corr_et < 35.) continue;
    if(photon->photonID("mvaPhoID-RunIIFall17-v2-wp80") == 0) continue; //WP80

    //float abseta = fabs(photon->superCluster()->eta());
    //float eA = effectiveAreas_ph_.getEffectiveArea(abseta);

    //photon_iso = (pfIso.sumChargedHadronPt + std::max( 0.0f, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - eA*rho_))/photon->et();

    //if(photon->chargedHadronIso()/corr_et > 0.3 || photon->photonIso() > 4.) continue; //|| photon->trackIso() > 6

    nPhotons38WP80++;

    // Apply energy scale corrections
    ph_p4 = photon->p4();// * photon->userFloat("ecalEnergyPostCorr") / photon->energy();

    if(corr_et < eTphMax) continue; //choose as best photon the one with highest eT
    
    eTphMax = corr_et;
    //ph_iso_ChargedHadron = photon->chargedHadronIso();
    //ph_iso_NeutralHadron = photon->neutralHadronIso();
    ph_iso_Photon        = photon->photonIso();
    //ph_iso_eArho         = eA * rho_;
    ph_eT     = corr_et;
    //ph_energy = photon->energy();
    ph_eta    = photon->eta();
    ph_etaSC  = photon->superCluster()->eta();
    ph_phi    = photon->phi();

    photonRegressionError = photon->getCorrectedEnergyError(reco::Photon::P4type::regression2);
    if(debug) cout << "Regression2 Energy Error: " << photonRegressionError << endl;
    
    cand_photon_found = true;
    nPhotonsChosen++;

  }//PHOTON FORLOOP END -------------------------------------------------------------------------------------------------------------------------

  //Return if there are no photons chosen
  if(!cand_photon_found) {
    cout<<"No best photon found, RETURN."<<endl;
    return;
  } 

  _Nevents_isPhoton++;
  
  mytree->Fill();
}

//*************************************************************//
//                                                             //
//---------------------- Create the tree ----------------------//
//                                                             //
//*************************************************************//

void ZMesonGamma::create_trees()
{
  mytree = fs->make<TTree>("mytree", "Tree containing gen&reco");
  
  mytree->Branch("isTwoProngTrigger",&isTwoProngTrigger);

  mytree->Branch("nPhotons38WP80",&nPhotons38WP80);
  mytree->Branch("nPhotons20WP90",&nPhotons20WP90);
  mytree->Branch("nPhotonsChosen",&nPhotonsChosen);

  mytree->Branch("photon_eT",&ph_eT);
  mytree->Branch("photon_eta",&ph_eta);
  mytree->Branch("photon_etaSC",&ph_etaSC);
  mytree->Branch("photon_phi",&ph_phi);

  
    
}

void ZMesonGamma::beginJob()
{
  //Flag for PileUp reweighting
  if (!runningOnData_){ // PU reweighting for 2017
    Lumiweights_ = edm::LumiReWeighting("MCpileUp_2018_25ns_UltraLegacy_PoissonOOTPU.root", "MyDataPileupHistogram.root", "pileup", "pileup");
  }
}

void ZMesonGamma::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_triggered);
  h_Events->Fill(2.5,_Nevents_isPhoton);
  h_Events->GetXaxis()->SetBinLabel(1,"processed");
  h_Events->GetXaxis()->SetBinLabel(2,"triggered");
  h_Events->GetXaxis()->SetBinLabel(3,"best photon");


}

//define this as a plug-in
DEFINE_FWK_MODULE(ZMesonGamma);
