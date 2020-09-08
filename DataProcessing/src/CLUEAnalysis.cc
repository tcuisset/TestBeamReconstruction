#include "UserCode/DataProcessing/interface/CLUEAnalysis.h"

//Simplified with respect to the CMSSW version, since
//   1) Only considers the EE section
//   2) Has W0=2.9 hardcoded
//Returns x and y position of the cluster
void CLUEAnalysis::calculatePositionsAndEnergy(const std::vector<float>& xpos, const std::vector<float>& ypos, const std::vector<float>& weights, const std::vector<int>& clusterid, const std::vector<int>& layerid) {
  assert(!xpos.empty() && !ypos.empty() && !weights.empty() && !clusterid.empty()&& !layerid.empty());
  auto start = std::chrono::high_resolution_clock::now();
  const int nclusters = ( *( std::max_element(clusterid.begin(), clusterid.end()) ) 
			  + 1 /*cluster index starts at zero*/  + 1 /*outliers*/ );
  std::vector<float> total_weight(nclusters, 0.);
  std::vector<float> total_weight_log(nclusters, 0.) ;
  std::vector<float> x(nclusters, 0.);
  std::vector<float> y(nclusters, 0.);
  std::vector<float> layers(nclusters, 0.);
  for(auto i: util::lang::indices(weights))
    {
      unsigned int weight_index = clusterid.at(i) + 1; //outliers will correspond to total_weight[0]
      total_weight.at(weight_index) += weights.at(i);
    }
  en_ = total_weight; //copy

  for (auto i: util::lang::indices(weights))
    {
      unsigned int weight_index = clusterid.at(i) + 1; //outliers will correspond to total_weight[0]
      float Wi = std::max(W0_ + std::log(weights.at(i) / total_weight.at(weight_index)), 0.f);
      x.at(weight_index) += xpos.at(i) * Wi;
      y.at(weight_index) += ypos.at(i) * Wi;
      total_weight_log.at(weight_index) += Wi;
      if(layers[weight_index] == 0 and weight_index != 0) 
	layers[weight_index] = layerid.at(i); //all hits in a cluster will belong to the same layer
    }
  for(auto i : util::lang::indices(total_weight_log))
    {
      if (total_weight_log.at(i) != 0.) {
	float inv_tot_weight_log = 1.f / total_weight_log.at(i);
	pos_.push_back( std::make_tuple( x.at(i) * inv_tot_weight_log, y.at(i) * inv_tot_weight_log, layers.at(i) ) );
      } 
      else
	pos_.push_back( std::make_tuple(0.f, 0.f, 0.f) );
    }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  //std::cout << "--- calculatePositions      " << elapsed.count() *1000 << " ms\n\n";
}


//Simplified with respect to the CMSSW version, since
//   1) Only considers the EE section
//   2) Has W0=2.9 hardcoded
//Returns x and y position of the cluster
void CLUEAnalysis::calculatePositions(const std::vector<float>& xpos, const std::vector<float>& ypos, const std::vector<float>& weights, const std::vector<int>& clusterid, const std::vector<int>& layerid) {
  auto start = std::chrono::high_resolution_clock::now();

  const int nclusters = *( std::max_element(clusterid.begin(), clusterid.end()) ) + 1; //number of clusters plus outliers
  std::vector<float> total_weight(nclusters, 0.);
  std::vector<float> total_weight_log(nclusters, 0.) ;
  std::vector<float> x(nclusters, 0.);
  std::vector<float> y(nclusters, 0.);
  std::vector<float> layers(nclusters, 0.);

  for(auto i: util::lang::indices(weights))
    {
      unsigned int weight_index = clusterid[i] + 1; //outliers will correspond to total_weight[0]
      total_weight[weight_index] += weights[i];
    }

  for (auto i: util::lang::indices(weights))
    {
      unsigned int weight_index = clusterid[i] + 1; //outliers will correspond to total_weight[0]
      float Wi = std::max(W0_ + std::log(weights[i] / total_weight[weight_index]), 0.f);
      x[weight_index] += xpos[i] * Wi;
      y[weight_index] += ypos[i] * Wi;
      total_weight_log[weight_index] += Wi;
      if(layers[weight_index] == 0) 
	layers[weight_index] = layerid[i]; //all hits in a cluster will belong to the same layer
    }

  for(auto i : util::lang::indices(total_weight_log))
    {
      if (total_weight_log[i] != 0.) {
	float inv_tot_weight_log = 1.f / total_weight_log[i];
	pos_.push_back( std::make_tuple( x[i] * inv_tot_weight_log, y[i] * inv_tot_weight_log, layers[i] ) );
      } 
      else
	pos_.push_back( std::make_tuple(0.f, 0.f, 0.f) );
    }

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "--- calculatePositions      " << elapsed.count() *1000 << " ms\n\n";
}

void CLUEAnalysis::calculateEnergy( const std::vector<float>& weights, const std::vector<int>& clusterid ) {
  auto start = std::chrono::high_resolution_clock::now();

  const int nclusters = *( std::max_element(clusterid.begin(), clusterid.end()) ) + 1; //number of clusters plus outliers
  std::vector<float> total_weight(nclusters, 0.);

  for(auto i: util::lang::indices(weights))
    {
      unsigned int weight_index = clusterid[i] + 1; //outliers will correspond to total_weight[0]
      total_weight[weight_index] += weights[i];
    }

  en_ = total_weight;

  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "--- calculateEnergy      " << elapsed.count() *1000 << " ms\n\n";
}

//calculate the number of clusterized hits and clusterized energy per layer
void CLUEAnalysis::calculateLayerDepVars(const std::vector<float>& xpos, const std::vector<float>& ypos, const std::vector<float>& weights, const std::vector<int>& clusterid, const std::vector<int>& layerid, const std::vector<float>& rhos, const std::vector<float>& deltas, const std::vector<bool>& seeds, const std::vector<unsigned int>& nhitsincluster) {
  assert(!weights.empty() && !clusterid.empty() && !layerid.empty() && !rhos.empty() && !deltas.empty());

  //calculate the number of rechits and clusterized energy per layer
  std::array<unsigned int, detectorConstants::nlayers_emshowers> hits_per_layer = {{0}};
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> en_per_layer;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> rhos_per_layer;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> deltas_per_layer;
  std::array< std::vector<bool>,  detectorConstants::nlayers_emshowers> seeds_per_layer;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> xpos_per_layer;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> ypos_per_layer;
  std::array< std::vector<unsigned int>, detectorConstants::nlayers_emshowers> cluster_size_per_layer; //number of hits in the cluster for each hit and layer

  //unsigned int this_cluster_size = 1;
  //int old_i = -1;
  for(auto i: util::lang::indices(weights))
    {
      if(clusterid[i] != -1)  //outliers are not considered
	{
	  int layeridx = layerid[i]-1; //layers start at 1

	  //////////////////////////////////////////////////////////////////////////////////////
	  /*Algorithm for finding the number of hits in the cluster to which each hit is associated*/

	  /*
	  if( i >= old_i + this_cluster_size ) //'==' does not work when clusterid == -1
	    {
	      this_cluster_size = 1;
	      unsigned int next_idx = i+1;
	      if(next_idx < weights.size()) //if it overflows the cluster size is 1 (last hit in the vector)
		{
		  while( clusterid[i] == clusterid[next_idx])
		    {
		      this_cluster_size += 1;
		      ++next_idx;
		      if( next_idx >= weights.size()) //exit if it overflows
			break;
		    }
		}
	      for(unsigned j=0; j<this_cluster_size; ++j)
		cluster_size_per_layer.at(layeridx).push_back( this_cluster_size );
	      old_i = i;
	    }
	  */
	  //////////////////////////////////////////////////////////////////////////////////////
	  
	  hits_per_layer.at(layeridx) += 1;
	  en_per_layer.at(layeridx).push_back( weights[i] );
	  rhos_per_layer.at(layeridx).push_back( rhos[i] );
	  deltas_per_layer.at(layeridx).push_back( deltas[i] );
	  seeds_per_layer.at(layeridx).push_back( seeds[i] );
	  xpos_per_layer.at(layeridx).push_back( xpos[i] );
	  ypos_per_layer.at(layeridx).push_back( ypos[i] );
	  cluster_size_per_layer.at(layeridx).push_back( nhitsincluster[i] );
	}
      //Note: We should get an out-of-bounds error for trying to access info at layers > 28. 
      //      It does not happen since all hits not in the CEE are marked as outliers by CLUE (clusterid == -1).
    }
  //fill std::array with fractions
  for(unsigned int ilayer=0; ilayer<detectorConstants::nlayers_emshowers; ++ilayer)
    layerdep_vars_.at(ilayer) = std::make_tuple(hits_per_layer.at(ilayer), en_per_layer.at(ilayer), 
						rhos_per_layer.at(ilayer), deltas_per_layer.at(ilayer), seeds_per_layer.at(ilayer),
						xpos_per_layer.at(ilayer), ypos_per_layer.at(ilayer),
						cluster_size_per_layer.at(ilayer)); 
}

//calculate the number of clusterized hits and clusterized energy per layer and per cluster
void CLUEAnalysis::calculateClusterDepVars(const std::vector<float>& xpos, const std::vector<float>& ypos, const std::vector<float>& weights, const std::vector<int>& clusterid, const std::vector<int>& layerid) {
  assert(!weights.empty() && !clusterid.empty() && !layerid.empty());

  std::array< unsigned int, detectorConstants::nlayers_emshowers> nclusters_per_layer = {{0}}; //number of clusters per layer for resizing the vectors
  std::vector<unsigned int> clusterIds;
  std::unordered_map<unsigned int, unsigned int> clusterIndexMap;

  //calculate number of clusters per layer
  for(auto i: util::lang::indices(weights))
    {
      if(clusterid[i]!=-1 /*outliers are not considered*/ and clusterIndexMap.find(clusterid[i]) == clusterIndexMap.end())
	{
	  int layeridx = layerid[i]-1; //layers start at 1
	  clusterIndexMap.emplace(clusterid[i], nclusters_per_layer[layeridx]);
	  nclusters_per_layer[layeridx] += 1;
	}
    }
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> en_per_cluster;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers> en_per_cluster_log; //helper for calculating the positions
  std::array< std::vector<unsigned int>, detectorConstants::nlayers_emshowers > hits_per_cluster;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers > x_per_cluster;
  std::array< std::vector<float>, detectorConstants::nlayers_emshowers > y_per_cluster;
  for(unsigned int i=0; i<detectorConstants::nlayers_emshowers; ++i) 
    {
      en_per_cluster[i].resize( nclusters_per_layer[i], 0.f ); //resizes and default-initializes to zero
      en_per_cluster_log[i].resize( nclusters_per_layer[i], 0.f ); //resizes and default-initializes to zero
      hits_per_cluster[i].resize( nclusters_per_layer[i], 0 ); //resizes and default-initializes to zero
      x_per_cluster[i].resize( nclusters_per_layer[i], 0.f ); //resizes and default-initializes to zero
      y_per_cluster[i].resize( nclusters_per_layer[i], 0.f ); //resizes and default-initializes to zero
    }

  for(auto i: util::lang::indices(weights))
    {
      if(clusterid[i] != -1)  //outliers are not considered
	{
	  unsigned int layeridx = static_cast<unsigned int>(layerid[i]) - 1; //layers start at 1
	  unsigned int vectoridx = clusterIndexMap[clusterid[i]];
	  en_per_cluster.at(layeridx).at(vectoridx) += weights[i];
	  hits_per_cluster.at(layeridx).at(vectoridx) += 1;
	}
      //Note: We should get an out-of-bounds error for trying to access info at layers > 28. 
      //      It does not happen since all hits not in the CEE are marked as outliers by CLUE (clusterid == -1).
    }
  
   for (auto i: util::lang::indices(weights))
    {
      if(clusterid[i] != -1)  //outliers are not considered
	{
	  unsigned int layeridx = static_cast<unsigned int>(layerid[i]) - 1; //layers start at 1
	  unsigned int vectoridx = clusterIndexMap[clusterid[i]];
	  float Wi = std::max(W0_ + std::log(weights[i] / en_per_cluster.at(layeridx).at(vectoridx)), 0.f);
	  x_per_cluster.at(layeridx).at(vectoridx) += xpos[i] * Wi;
	  y_per_cluster.at(layeridx).at(vectoridx) += ypos[i] * Wi;
	  en_per_cluster_log.at(layeridx).at(vectoridx) += Wi;
	}
    }

   for(unsigned int s1=0; s1<en_per_cluster_log.size(); ++s1) {
     for(unsigned int s2=0; s2<en_per_cluster_log[s1].size(); ++s2)
       {
	 if (en_per_cluster_log.at(s1).at(s2) != 0.)
	   {
	     float inv_log = 1.f / en_per_cluster_log.at(s1).at(s2);
	     x_per_cluster.at(s1).at(s2) *= inv_log;
	     y_per_cluster.at(s1).at(s2) *= inv_log;
	   }
	 else
	   {
	     x_per_cluster.at(s1).at(s2) = -99.f;
	     y_per_cluster.at(s1).at(s2) = -99.f;
	   }
       }
   }

   
  //fill std::array with clusterized cluster number of hits and energy
  //both vectors might well be empty, in case there was no cluster in a particular layer
  for(unsigned int ilayer=0; ilayer<detectorConstants::nlayers_emshowers; ++ilayer)
    clusterdep_vars_.at(ilayer) = std::make_tuple( hits_per_cluster.at(ilayer), en_per_cluster.at(ilayer), x_per_cluster.at(ilayer), y_per_cluster.at(ilayer));
}

//Returns quantities of interest of individual clusters
std::vector<dataformats::data> CLUEAnalysis::getTotalPositionsAndEnergyOutput(std::string& outputFileName, bool verbose) {
  bool pos_set = !pos_.empty(), en_set = !en_.empty();

  if(verbose) 
    {
      std::ofstream oFile(outputFileName, std::ios::out);
      if(pos_set && en_set)
	assert(en_.size() == pos_.size());
      oFile << "###Cluster Analysis: individual clusters (posx, posy, posz, energy)\n";
      for(unsigned int i=0; i<pos_.size(); ++i)
	{
	  if(pos_set)
	    oFile << std::get<0>(pos_[i]) << "," << std::get<1>(pos_[i]) << "," << std::get<2>(pos_[i]);
	  if(en_set)
	    oFile << en_[i];
	  oFile << std::endl;
	}
      oFile.close();
    }
  
  std::vector<dataformats::data> d;
  for(auto i : util::lang::indices(pos_))
    {
      float posx=0., posy=0., posz=0., en=0.;
      if(pos_set)
	{
	  posx = std::get<0>(pos_[i]);
	  posy = std::get<1>(pos_[i]);
	  posz = std::get<2>(pos_[i]);
	}
      if(en_set)
	en = en_[i];
      d.push_back( std::make_tuple(posx,posy,posz,en) );
    }
  return d;
}

//Returns the sum of all quantities of interest of individual clusters
float CLUEAnalysis::getTotalEnergyOutput(const std::string& outputFileName, bool verbose) {
  assert( !en_.empty() );
  float toten = 0.;
  for(unsigned int i=1; i<en_.size(); ++i) //outliers (i=0) are skipped!
    toten += en_[i];
  if(verbose) 
    {
      std::ofstream oFile(outputFileName, std::ios::out);
      oFile << "###Cluster Analysis: total energy\n";
      oFile << toten << std::endl;
      oFile.close();
    }
  return toten;
}

//Returns the number of clusterized hits and clusterized energy per layer
dataformats::layervars CLUEAnalysis::getTotalLayerDepOutput() {
  return this->layerdep_vars_;
}

//Returns the number of clusterized hits and clusterized energy per layer and per cluster
std::array< std::tuple< std::vector<unsigned int>, std::vector<float>, std::vector<float>, std::vector<float> >, detectorConstants::nlayers_emshowers> CLUEAnalysis::getTotalClusterDepOutput() {
  return this->clusterdep_vars_;
}
