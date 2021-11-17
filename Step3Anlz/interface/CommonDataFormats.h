#ifndef HGCAL_RECO_ANALYSIS_HGCMLANALYZER_H
#define HGCAL_RECO_ANALYSIS_HGCMLANALYZER_H

#include <vector>
#include "DataFormats/DetId/interface/DetId.h"

struct caloparticle {
  int pdgid_;
  float energy_;
  float pt_;
  float eta_;
  float phi_;
  std::vector<DetId> rechitdetid_;
  std::vector<float> rechitenergy_;
};

struct layercluster {
  float energy_;
  float eta_;
  float phi_;
  float x_;
  float y_;
  float z_;
  int nrechits_;
  int layer_;
  int idx2Trackster_;
};

struct trackster {
  int idx_;
  int type_;  // pdgid
  float energy_;
  float eta_;
  float phi_;
  float x_;
  float y_;
  float z_;
  float pcaeigval0_;
  float pcasig0_;
  float pcaeigval1_;
  float pcasig1_;
  float pcaeigval2_;
  float pcasig2_;
  float cpenergy_;
  float cpeta_;
  int cppdgid_;
};

#endif
