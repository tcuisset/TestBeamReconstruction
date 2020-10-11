#ifndef CLUEAnalysis_h
#define CLUEAnalysis_h

// C/C++ headers
#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <cmath>
#include <chrono>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

//loop over two vectors
#include "UserCode/DataProcessing/interface/range.h"

enum SHOWERTYPE { EM, HAD };
enum DATATYPE { DATA, MC };

namespace detectorConstants {
  constexpr std::array<float, 2> energyDepositedByMIP = {{0.0850, 0.0567}}; //value given by Thorben [MeV]
  constexpr float sigmaNoiseSiSensor = energyDepositedByMIP[0] / 6.f; //value given by Thorben [MeV]
  constexpr unsigned int nlayers_emshowers = 28;
  constexpr unsigned int totalnlayers = 40;
  constexpr unsigned int layerBoundary = 26;
  constexpr std::array<float, nlayers_emshowers> dEdX = {{11.289,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,11.360,11.360,11.360,11.360,10.995,10.995,11.153,7.470}}; //values in DN-19-019 [MeV/MIP]

  //Shubham: https://indico.cern.ch/event/923102/contributions/3878372/attachments/2069945/3474694/PionAnalysis_EE_FH_AHCAL_energyScale_short_07072020.pdf (slide #5)
  constexpr float globalWeightCEE = 10.6; // [MeV/MIP]
  constexpr float globalWeightCEH = 78.9; // [MeV/MIP]
  constexpr float globalWeightRelative = 0.4; // [dimensionless]

  //Shubham: https://indico.cern.ch/event/923097/contributions/3878337/attachments/2048763/3433394/pion_analysis_shower_start_finder_algorithm_optimization_2June2020.pdf (slide #16)
  //constexpr std::unordered_map<int, int> Ethresh = {{20,12}, {50,20}, {80,25}, {100,30}, {120,30}, {200,40}, {250,40}, {300,40}};
}

namespace dataformats {
  using data = std::tuple<float,float,float,float>; //posx, posy, poz/layer and energy
  using position = std::tuple<float,float,float>; //posx, posy and poz/layer
  using layerfracs = std::vector< std::tuple<float, float> >;
  using layervars = std::vector< std::tuple<unsigned int, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<bool>, std::vector<float>, std::vector<float>, std::vector<unsigned int>> >;
  using layerhitvars = std::vector< std::tuple<std::vector<float>, std::vector<float>, std::vector<float>, std::vector<bool>, std::vector<float>, std::vector<float>, std::vector<unsigned int>> >;
  using clustervars = std::vector< std::tuple< std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>> >;
}

class CLUEAnalysis {
  /*Outliers are all identified to the 'cluster' of index = 0*/
private:
  SHOWERTYPE showertype;
  unsigned lmax;
  const float W0_ = 2.3f, dpos_ = 1.3f; //tunable parameters for cluster position measurement
  std::vector<dataformats::position> pos_;
  std::vector< float > en_;
  dataformats::layervars layerdep_vars_; //clusterized nhits and clusterized energy per event and per layer, densities, distances, isSeed boolean flag, x position and y position
  dataformats::clustervars clusterdep_vars_; //clusterized nhits and clusterized energy per event, per layer and per cluster
  std::vector<float> frac_clust_hits_;

  float hit_distance(const float&, const float&, const float&, const float&);
  
public:
  CLUEAnalysis(const SHOWERTYPE&);
  unsigned getLayerMax() {return lmax;}
  void calculateEnergy(const std::vector<float>&, const std::vector<int>&);
  void verboseResults(std::string&);
  void calculateLayerDepVars(const std::vector<float>&, const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const std::vector<int>&, const std::vector<float>&, const std::vector<float>&, const std::vector<bool>&, const std::vector<unsigned int>&);
  void calculateClusterDepVars(const std::vector<float>&, const std::vector<float>&, const std::vector<float>&, const std::vector<int>&, const std::vector<int>&, const std::vector<float>&, const std::vector<float>&);
  std::vector<dataformats::data> getTotalPositionsAndEnergyOutput(std::string& outputFileName, bool verbose=0);
  float getTotalEnergyOutput(const std::string& outputFileName, bool verbose=0);
  dataformats::layervars getTotalLayerDepOutput();
  dataformats::clustervars getTotalClusterDepOutput();
};

#endif //CLUEAnalysis_h
