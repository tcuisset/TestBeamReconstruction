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
 public:
  CLUEAnalysis(){};
  void calculatePositionsAndEnergy(const std::vector<float>&, const std::vector<float>&, 
				   const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculatePositions(const std::vector<float>&, const std::vector<float>&, 
			    const std::vector<float>&, const std::vector<int>&, const std::vector<int>&);
  void calculateEnergy(const std::vector<float>&, const std::vector<int>&);
  void verboseResults(std::string&);
  std::vector<dataformats::data> getIndividualClusterOutput(std::string& outputFileName, bool verbose=0);
  float getTotalClusterOutput(const std::string& outputFileName, bool verbose=0);

 private:
  /*
  constexpr nlayers_ = 50;
  TTree *t_;
  std::array< float, nlayers_ > outx_;
  std::array< float, nlayers_ > outy_;
  std::array< TBranch*, nlayers_ > branchesx_;
  std::array< TBranch*, nlayers_ > branchesy_;
  */
  std::vector<dataformats::position> pos_;
  std::vector< float > en_;
};

#endif //CLUEAnalysis_h
