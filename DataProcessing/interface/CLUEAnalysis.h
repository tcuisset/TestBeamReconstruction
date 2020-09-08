#ifndef CLUEAnalysis_h
#define CLUEAnalysis_h

// C/C++ headers
#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <chrono>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

//loop over two vectors
#include "UserCode/DataProcessing/interface/range.h"

namespace detectorConstants {
  constexpr std::array<float, 2> energyDepositedByMIP = {{0.0850, 0.0567}}; //value given by Thorben [MeV]
  constexpr float sigmaNoiseSiSensor = energyDepositedByMIP[0] / 6.f; //value given by Thorben [MeV]
  constexpr unsigned int nlayers_emshowers = 28;
  constexpr unsigned int totalnlayers = 40;
  constexpr unsigned int layerBoundary = 26;
  constexpr std::array<float, nlayers_emshowers> dEdX = {{11.289,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,11.360,11.360,11.360,11.360,10.995,10.995,11.153,7.470}}; //values in DN-19-019 [MeV/MIP]  
}


namespace dataformats {
  using data = std::tuple<float,float,float,float>; //posx, posy, poz/layer and energy
  using position = std::tuple<float,float,float>; //posx, posy and poz/layer
  using layerfracs = std::array< std::tuple<float, float>, detectorConstants::nlayers_emshowers>;
  using layervars = std::array< std::tuple<unsigned int, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<bool>, std::vector<float>, std::vector<float>, std::vector<unsigned int>>, detectorConstants::nlayers_emshowers>;
  using layerhitvars = std::array< std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<bool>, std::vector<float>, std::vector<float>, std::vector<unsigned int>>, detectorConstants::nlayers_emshowers>;
  using clustervars = std::array< std::tuple< std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float> >, detectorConstants::nlayers_emshowers>;
}

class CLUEAnalysis {
  /*Outliers are all identified to the 'cluster' of index = 0*/
private:
  const float W0_ = 2.9f;
  std::vector<dataformats::position> pos_;
  std::vector< float > en_;
  dataformats::layervars layerdep_vars_; //clusterized nhits and clusterized energy per event and per layer, densities, distances, isSeed boolean flag, x position and y position
  std::array< std::tuple< std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float> >, detectorConstants::nlayers_emshowers> clusterdep_vars_; //clusterized nhits and clusterized energy per event, per layer and per cluster
  std::vector<float> frac_clust_hits_;
  
public:
  CLUEAnalysis(){};
  void calculatePositionsAndEnergy(const std::vector<float>&, const std::vector<float>&, 
				   const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculatePositions(const std::vector<float>&, const std::vector<float>&, 
			    const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculateEnergy(const std::vector<float>&, const std::vector<int>&);
  void verboseResults(std::string&);
  void calculateLayerDepVars(const std::vector<float>&, const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const std::vector<int>&, const std::vector<float>&, const std::vector<float>&, const std::vector<bool>&, const std::vector<unsigned int>&);
  void calculateClusterDepVars(const std::vector<float>&, const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  std::vector<dataformats::data> getTotalPositionsAndEnergyOutput(std::string& outputFileName, bool verbose=0);
  float getTotalEnergyOutput(const std::string& outputFileName, bool verbose=0);
  dataformats::layervars getTotalLayerDepOutput();
  std::array< std::tuple< std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float> >, detectorConstants::nlayers_emshowers> getTotalClusterDepOutput();
};

#endif //CLUEAnalysis_h
