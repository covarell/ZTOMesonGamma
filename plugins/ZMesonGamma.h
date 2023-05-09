//---------- class declaration----------

class ZMesonGamma : public edm::EDAnalyzer {

 public:
  explicit ZMesonGamma(const edm::ParameterSet&);
  ~ZMesonGamma();

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() override;

  bool runningOnData_;
  const edm::InputTag packedPFCandidates_;
  const edm::InputTag prunedGenParticles_;
  const edm::InputTag slimmedPhotons_;
  const edm::InputTag slimmedElectrons_;
  const edm::InputTag slimmedMuons_; 
  const edm::InputTag slimmedJets_;
  const edm::InputTag GenInfo_;


  edm::LumiReWeighting Lumiweights_;



  edm::Service<TFileService> fs;

  void create_trees();

  // ---------- member data ----------- //
  TH1F* h_Events;
  TH1F* h_pileup;

  //Counters
  int nPV;

  int nMuons10;
  int nMuons20;
  int nElectrons10;
  int nElectrons20;
  int nPhotonsChosen;
  int nPhotons20WP90;
  int nPhotons38WP80;
  int nJets30;
  int nJets25;

  int _Nevents_triggered;
  int _Nevents_processed;
  int _Nevents_isPhoton;
  int _Nevents_isTwoKaons;


  //bools
  bool isTwoProngTrigger;

  //debug
  bool debug;
  bool verbose;

  //TTree and TTree variables
  TTree *mytree;

  int run_number;
  int event_number;


  float ph_eT;
  float ph_energy;
  float ph_eta;
  float ph_etaSC;
  float ph_phi;
  float eTphMax;
  float photonRegressionError;
  float ph_iso_ChargedHadron;
  float ph_iso_NeutralHadron;
  float ph_iso_Photon;
  float ph_iso_eArho;

  float firstCandPx;
  float firstCandPy;
  float firstCandPz;
  float secondCandPx;
  float secondCandPy;
  float secondCandPz;
  float firstCandEnergy;
  float secondCandEnergy;
  float firstCandEnergy_K;
  float secondCandEnergy_K;
  float firstCandEnergy_Pi;
  float secondCandEnergy_Pi;
  float _firstCandPt;
  float _firstCandEta;
  float _firstCandPhi;
  float _firstCandCharge;
  float _secondCandPt;
  float _secondCandEta;
  float _secondCandPhi;
  float _secondCandCharge;
  float _bestCouplePt;
  float _bestCoupleEta;
  float _bestCouplePhi;

  float _Jet_Photon_invMass;
  float _MesonMass;

  float met_pT;
  float metpuppi_pT;

  //Jet datamember  
  float _jet_invMass;
  float _bestJet_pT;
  float _bestJet_eta;
  float _bestJet_phi;
  int _bestJet_nDaughters;
  float _bestJet_pTMax;
  float _bestJet_chargedEmEnergy;
  float _bestJet_neutralEmEnergy;
  float _bestJet_chargedHadEnergy;
  float _bestJet_neutralHadEnergy;
  float _bestJet_chargedEmEnergyFraction;
  float _bestJet_neutralEmEnergyFraction;
  float _bestJet_chargedHadEnergyFraction;
  float _bestJet_neutralHadEnergyFraction;
  int _bestJet_chargedHadMultiplicity;
  float _bestJet_invMass;
  float _bestJet_Photon_invMass;
  float _bestJet_JECunc;

  //MC truth
  float PU_Weight;
  float MC_Weight;

  float rho_;

  //Tokens
  edm::EDGetTokenT<std::vector<reco::GenParticle> > prunedGenParticlesToken_; 
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > packedPFCandidatesToken_; 
  edm::EDGetTokenT<std::vector<reco::Vertex> > offlineSlimmedPrimaryVerticesToken_;  
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
  edm::EDGetTokenT<GenEventInfoProduct> GenInfoToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
  edm::EDGetTokenT<std::vector<pat::Muon> > slimmedMuonsToken_;   
  edm::EDGetTokenT<std::vector<pat::Jet> > slimmedJetsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > slimmedMETsPuppiToken_;
  //Ele ID decisions objects
  edm::EDGetToken electronsMiniAODToken_;

  //Photon ID decisions
  edm::EDGetToken photonsMiniAODToken_;

  bool verboseIdFlag_;


};