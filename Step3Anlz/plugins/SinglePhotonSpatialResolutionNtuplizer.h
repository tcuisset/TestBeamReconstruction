// system include files
#include <memory>
#include <algorithm>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/HGCRecHit/interface/HGCRecHit.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"
#include "RecoLocalCalo/HGCalRecProducers/interface/HGCalClusteringAlgoBase.h"

// from HGC Validator code
#include "Validation/HGCalValidation/interface/HGVHistoProducerAlgo.h"

//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1.h>

#include "UserCode/Step3Anlz/interface/CommonDataFormats.h"

class SinglePhotonSpatialResolutionNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SinglePhotonSpatialResolutionNtuplizer(const edm::ParameterSet&);
  ~SinglePhotonSpatialResolutionNtuplizer() override {};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  void initialize_branches();
  void clear();
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

  // ----------member data ---------------------------
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertexesToken_;

  TTree* tree_;

  edm::RunNumber_t irun;
  edm::EventNumber_t ievent;
  edm::LuminosityBlockNumber_t ilumiblock;
  edm::Timestamp itime;

  float beam_energy_ = 50.f;
  unsigned nhits_;
  int nlayers = 28;

  std::vector<float> rechit_energy_;
  std::vector<unsigned int> rechit_detid_;
  std::vector<unsigned int> rechit_layer_;
  std::vector<float> rechit_x_;
  std::vector<float> rechit_y_;
  std::vector<float> rechit_z_;
  std::vector<float> impactX_;
  std::vector<float> impactY_;

  const HGCalDDDConstants* ddd_;
  const HGCeeRecHitCollection* eeRecHits = nullptr;
  const std::vector<SimVertex>* simVertexes = nullptr;
};
