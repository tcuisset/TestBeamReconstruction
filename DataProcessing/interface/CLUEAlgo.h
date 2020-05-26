#ifndef CLUEAlgo_h
#define CLUEAlgo_h

// C/C++ headers
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <functional>
#include <chrono>

#include "CLUEAnalysis.h"
#include "LayerTiles.h"
#include "Points.h"

class CLUEAlgo{

  public:
    // constructor
  CLUEAlgo(float dc, float kappa, float ecut, bool verbose=false ){ 
      dc_ = dc; 
      ecut_ = ecut;
      kappa_ = kappa;
      //dm_ =  std::max(deltao_, deltac_);
      outlierDeltaFactor_ = 2.;
      verbose_ = verbose;
    
    }
    // destructor
    ~CLUEAlgo(){} 
    
    // public variables
    float dc_, ecut_, kappa_, outlierDeltaFactor_;
    bool verbose_;
    
    Points points_;

    std::vector<float> getHitsClusterX();
    std::vector<float> getHitsClusterY();
    std::vector<float> getHitsWeight();
    std::vector<int> getHitsClusterId();
    std::vector<int> getHitsLayerId();
    std::vector<float> getHitsRho();
    std::vector<float> getHitsDistanceToHighest();

    // public methods
    //Note: The layer input and output (see getHitsLayerId()) start counting at 1, but the calculations inside use a 0-based index
    void setPoints(int n, float* x, float* y, unsigned int* layer, float* weight) {
      points_.clear();

      // input variables
      for(int i=0; i<n; ++i)
	{
	  if(layer[i] > detectorConstants::nlayers)
	    continue;
	  float endeposited_mip = layer[i] < 27 ? detectorConstants::energyDepositedByMIP[0] : detectorConstants::energyDepositedByMIP[1];
	  if( weight[i] < ecut_ * detectorConstants::sigmaNoiseSiSensor / endeposited_mip * detectorConstants::dEdX.at(layer[i]-1) )
	      continue;
	  points_.x.push_back(x[i]);
	  points_.y.push_back(y[i]);
	  points_.layer.push_back(layer[i]-1);
	  points_.weight.push_back(weight[i]);
	}

      points_.n = points_.x.size();

      // result variables
      points_.rho.resize(points_.n,0);
      points_.delta.resize(points_.n,std::numeric_limits<float>::max());
      points_.nearestHigher.resize(points_.n,-1);
      points_.isSeed.resize(points_.n,0);
      points_.followers.resize(points_.n);
      points_.clusterIndex.resize(points_.n,-1);
    }

    void clearPoints(){ points_.clear(); }

    void makeClusters();

    void verboseResults( std::string outputFileName = "cout", int nVerbose = -1) { 
      
      if (verbose_) {
      
        if (nVerbose ==-1) nVerbose=points_.n;

        // verbose to screen
        if (outputFileName.compare("cout") == 0 )  {
          std::cout << "index,x,y,layer,weight,rho,delta,nh,isSeed,clusterId"<< std::endl;
          for(int i = 0; i < nVerbose; i++) {
            std::cout << i << ","<<points_.x[i]<< ","<<points_.y[i]<< ","<<points_.layer[i]+1 << ","<<points_.weight[i];
            std::cout << "," << points_.rho[i];
            if (points_.delta[i] <= 999) 
              std::cout << ","<<points_.delta[i];
            else
              std::cout << ",999"; // convert +inf to 999 in verbose
            std::cout << ","<<points_.nearestHigher[i];
            std::cout << "," << points_.isSeed[i];
            std::cout << ","<<points_.clusterIndex[i];
            std::cout << std::endl;
          }
        }

        // verbose to file
        else{
          std::ofstream oFile(outputFileName);
          oFile << "index,x,y,layer,weight,rho,delta,nh,isSeed,clusterId\n";
          for(int i = 0; i < nVerbose; i++) {
            oFile << i << ","<<points_.x[i]<< ","<<points_.y[i]<< ","<<points_.layer[i]+1 << ","<<points_.weight[i];
            oFile << "," << points_.rho[i];
            if (points_.delta[i] <= 999) 
              oFile << ","<<points_.delta[i];
            else
              oFile << ",999"; // convert +inf to 999 in verbose
            oFile << ","<<points_.nearestHigher[i];
            oFile << "," << points_.isSeed[i];
            oFile << ","<<points_.clusterIndex[i];
            oFile << "\n";
          }
          oFile.close();
        }
      }// end of if verbose_
    }

    void verboseResults(  std::vector<unsigned int> rechits_id, std::string outputFileName = "cout", int nVerbose = -1) { 
      assert(points_.x.size() == rechits_id.size());

      if (verbose_) {
      
	if (nVerbose ==-1) nVerbose=points_.n;

	// verbose to screens
	if (outputFileName.compare("cout") == 0 )  {
	  std::cout << "index,rechit_id,x,y,layer,weight,rho,delta,nh,isSeed,clusterId"<< std::endl;
	  for(int i = 0; i < nVerbose; i++) {
	    std::cout << i << ","<<rechits_id[i]<<","<<points_.x[i]<< ","<<points_.y[i]<< ","<<points_.layer[i]+1 << ","<<points_.weight[i];
	    std::cout << "," << points_.rho[i];
	    if (points_.delta[i] <= 999) 
	      std::cout << ","<<points_.delta[i];
	    else
	      std::cout << ",999"; // convert +inf to 999 in verbose
	    std::cout << ","<<points_.nearestHigher[i];
	    std::cout << "," << points_.isSeed[i];
	    std::cout << ","<<points_.clusterIndex[i];
	    std::cout << std::endl;
	  }
	}

	// verbose to file
	else{
	  std::ofstream oFile(outputFileName);
	  oFile << "index,rechit_id,x,y,layer,weight,rho,delta,nh,isSeed,clusterId\n";
	  for(int i = 0; i < nVerbose; i++) {
	    oFile << i << ","<<rechits_id[i]<< ","<<points_.x[i]<< ","<<points_.y[i]<< ","<<points_.layer[i]+1 << ","<<points_.weight[i];
	    oFile << "," << points_.rho[i];
	    if (points_.delta[i] <= 999) 
	      oFile << ","<<points_.delta[i];
	    else
	      oFile << ",999"; // convert +inf to 999 in verbose
	    oFile << ","<<points_.nearestHigher[i];
	    oFile << "," << points_.isSeed[i];
	    oFile << ","<<points_.clusterIndex[i];
	    oFile << "\n";
	  }
	  oFile.close();
	}
      }// end of if verbose_
        
    }
        
  private:
    // private member methods
    void prepareDataStructures(std::array<LayerTiles, detectorConstants::nlayers> & );
    void calculateLocalDensity(std::array<LayerTiles, detectorConstants::nlayers> & );
    void calculateDistanceToHigher(std::array<LayerTiles, detectorConstants::nlayers> & );
    void findAndAssignClusters();
    inline float distance(int , int) const ;
};

#endif
