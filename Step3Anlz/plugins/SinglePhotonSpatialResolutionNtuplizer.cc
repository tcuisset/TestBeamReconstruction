#include "UserCode/Step3Anlz/plugins/SinglePhotonSpatialResolutionNtuplizer.h"

SinglePhotonSpatialResolutionNtuplizer::SinglePhotonSpatialResolutionNtuplizer(const edm::ParameterSet& iConfig) :
  hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
  hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
  hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
  simVertexesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertexes")))
{
  recHitTools.reset(new hgcal::RecHitTools());

  usesResource("TFileService");
  edm::Service<TFileService> file;
  tree_ = file->make<TTree>("relevant_branches", "relevant_banches");
  initialize_branches();
}

void SinglePhotonSpatialResolutionNtuplizer::initialize_branches()
{
  tree_->Branch("NRechits", &nhits_);
  tree_->Branch("beamEnergy", &beam_energy_);

  //std::vector
  tree_->Branch("ce_clean_energy_MeV", &rechit_energy_);
  tree_->Branch("ce_clean_detid", &rechit_detid_);
  tree_->Branch("ce_clean_layer", &rechit_layer_);
  tree_->Branch("ce_clean_x", &rechit_x_);
  tree_->Branch("ce_clean_y", &rechit_y_);
  tree_->Branch("ce_clean_z", &rechit_z_);
  tree_->Branch("impactX_shifted", &impactX_); //"shifted" refers to the shift applied when considering testbeam data; here (CMSSW simulation) it makes no sense
  tree_->Branch("impactY_shifted", &impactY_);
}

void SinglePhotonSpatialResolutionNtuplizer::clear()
{
  nhits_ = 0;
  
  //std::vector
  rechit_energy_.clear();
  rechit_detid_.clear();
  rechit_layer_.clear();
  rechit_x_.clear();
  rechit_y_.clear();
  rechit_z_.clear();
  impactX_.clear();
  impactY_.clear();

  impactX_.reserve(nlayers);
  impactX_.reserve(nlayers);
}

void SinglePhotonSpatialResolutionNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  clear();

  edm::Handle<HGCRecHitCollection> recHitHandleEE;
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);
  eeRecHits = recHitHandleEE.product();

  //geometry
  edm::ESHandle<HGCalGeometry> hgceeGeomHandle;
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive", hgceeGeomHandle);
  ddd_ = &( hgceeGeomHandle->topology().dddConstants() );

  //rechits
  edm::ESHandle<CaloGeometry> geom; 
  iSetup.get<CaloGeometryRecord>().get(geom); 
  //recHitTools->setGeometry(*geom); UNCOMMENT!!!!!!!!!!!!
    
  //vertexes
  edm::Handle<std::vector<SimVertex>> simVertexesHandle;
  iEvent.getByToken(simVertexesToken_, simVertexesHandle);
  simVertexes = simVertexesHandle.product();

  for(int i = 0; i<nlayers; ++i) 
    {
      impactX_.push_back(simVertexes->at(0).position().X());
      impactY_.push_back(simVertexes->at(0).position().Y()); 
    }
  
  for (auto it = eeRecHits->begin(); it != eeRecHits->end(); ++it)
    {
      HGCSiliconDetId detid(it->detid());
      int layer = detid.layer();
      if(layer<0) //use only one endcap
	continue;

      ++nhits_;

      rechit_energy_.push_back( it->energy()*1e3 ); //convert to MeV
      rechit_detid_.push_back( it->detid() );
      rechit_layer_.push_back( static_cast<unsigned int>(layer) );
      
      std::pair<float, float> xy = std::make_pair(0.f,0.f); //UNCOMMENT!!!!! ddd_->locateCell(detid, false);
      rechit_x_.push_back(xy.first);
      rechit_y_.push_back(xy.second);
      rechit_z_.push_back( abs(ddd_->waferZ(layer, true)) );
    }
  tree_->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void SinglePhotonSpatialResolutionNtuplizer::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SinglePhotonSpatialResolutionNtuplizer::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SinglePhotonSpatialResolutionNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SinglePhotonSpatialResolutionNtuplizer);
