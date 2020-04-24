#ifndef CLUEAnalysis_h
#define CLUEAnalysis_h

// C/C++ headers
#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <vector>
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

class CLUEAnalysis {
  /*Outliers are all identified to the 'cluster' of index = 0*/
private:
  static const unsigned int nlayers_ = 28;
  std::vector<dataformats::position> pos_;
  std::vector< float > en_;
  std::array< std::tuple<unsigned int, float>, nlayers_> layerdep_vars_; //clusterized nhits and << ", " << lusterized energy per event
  
public:
  CLUEAnalysis(){};
  void calculatePositionsAndEnergy(const std::vector<float>&, const std::vector<float>&, 
				   const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculatePositions(const std::vector<float>&, const std::vector<float>&, 
			    const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculateEnergy(const std::vector<float>&, const std::vector<int>&);
  void verboseResults(std::string&);
  void calculateLayerDepVars(const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  std::vector<dataformats::data> getIndividualClusterOutput(std::string& outputFileName, bool verbose=0);
  float getTotalClusterEnergyOutput(const std::string& outputFileName, bool verbose=0);
  std::array< std::tuple<unsigned int, float>, nlayers_> getTotalClusterLayerDepOutput();
};

#endif //CLUEAnalysis_h
