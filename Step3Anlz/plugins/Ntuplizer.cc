#include "TestBeamReconstruction/Step3Anlz/plugins/Ntuplizer.h"

Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
  hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
  hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
  hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH")))
{

  recHitTools.reset(new hgcal::RecHitTools());

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree", "pidtree");
  void initialize_branches();
}

void Ntuplizer::initialize_branches()
{
  tree->Branch("run", &run, "run/I");
  tree->Branch("event", &event, "event/I");
  tree->Branch("lumi", &lumi, "lumi/I");
  tree->Branch("weight", &weight, "weight/F");
}

void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<HGCRecHitCollection> recHitHandleEE;
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  // init vars
  edm::ESHandle<CaloGeometry> geom; 
  iSetup.get<CaloGeometryRecord>().get(geom); 
  recHitTools->setGeometry(*geom);

  //tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void Ntuplizer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
