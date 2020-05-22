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

namespace dataformats {
  using data = std::tuple<float,float,float,float>; //posx, posy, poz/layer and energy
  using position = std::tuple<float,float,float>; //posx, posy and poz/layer
}

namespace detectorConstants {
    constexpr float sigmaNoise = 1e-3 * 86.f / 6.f; //value given by Thorben [MeV]
    constexpr unsigned int nlayers = 28; //needs to be static to be used as array dimension
}

class CLUEAnalysis {
  /*Outliers are all identified to the 'cluster' of index = 0*/
private:
  const float W0_ = 2.9f;
  std::vector<dataformats::position> pos_;
  std::vector< float > en_;
  std::array< std::tuple<unsigned int, float, std::vector<float>, std::vector<float>>, detectorConstants::nlayers> layerdep_vars_; //clusterized nhits and clusterized energy per event and per layer
  std::array< std::tuple< std::vector<unsigned int>, std::vector<float> >, detectorConstants::nlayers> clusterdep_vars_; //clusterized nhits and clusterized energy per event, per layer and per cluster
  
public:
  CLUEAnalysis(){};
  void calculatePositionsAndEnergy(const std::vector<float>&, const std::vector<float>&, 
				   const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculatePositions(const std::vector<float>&, const std::vector<float>&, 
			    const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculateEnergy(const std::vector<float>&, const std::vector<int>&);
  void verboseResults(std::string&);
  void calculateLayerDepVars(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&, 
			     const std::vector<float>&, const std::vector<float>&);
  void calculateClusterDepVars(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  std::vector<dataformats::data> getIndividualClusterOutput(std::string& outputFileName, bool verbose=0);
  float getTotalEnergyOutput(const std::string& outputFileName, bool verbose=0);
  std::array< std::tuple<unsigned int, float, std::vector<float>, std::vector<float>>, detectorConstants::nlayers> getTotalLayerDepOutput();
  std::array< std::tuple< std::vector<unsigned int>, std::vector<float> >, detectorConstants::nlayers> getTotalClusterDepOutput();
};

#endif //CLUEAnalysis_h
