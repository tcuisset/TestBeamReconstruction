#include "UserCode/DataProcessing/interface/analyzer.h"

Analyzer::Analyzer(const std::vector< std::string >& in_file_path, const std::string& in_tree_name, const float& dc, const float& kappa, const float& ecut, const SHOWERTYPE& st): dc_(dc), kappa_(kappa), ecut_(ecut), st_(st)
{
  nfiles_ = in_file_path.size();
  std::cout << "Number of files being processed: " << nfiles_ << std::endl;
  std::cout << "Files: " << std::endl;
  beam_energies_.resize(nfiles_, 0.f);
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      sanity_checks(in_file_path[i]);
      std::cout << "#" << std::to_string(i+1) << " " << in_file_path[i] << std::endl;
      names_.push_back( std::make_pair(in_file_path[i], in_tree_name) );
      en_total_.push_back( std::vector< std::tuple<float, float> >() );
      layer_fracs_.push_back( std::vector< dataformats::layerfracs >() );
      layer_hitvars_.push_back( std::vector< dataformats::layerhitvars >() );
      clusterdep_.push_back( std::vector< dataformats::clustervars >() );
    }
}

//Overloaded constructor for job submission. Each job processes one file only.
Analyzer::Analyzer(const std::string& in_file_path, const std::string& in_tree_name, const float& dc, const float& kappa, const float& ecut, const SHOWERTYPE& st): dc_(dc), kappa_(kappa), ecut_(ecut), st_(st)
{
  nfiles_ = 1;
  std::cout << "Number of files being processed: " << nfiles_ << std::endl;
  std::cout << "File: " << in_file_path << std::endl;
  sanity_checks(in_file_path);
  names_.push_back( std::make_pair(in_file_path, in_tree_name) );
  beam_energies_.resize(1, 0.f);
  en_total_.push_back( std::vector< std::tuple<float, float> >() );
  layer_fracs_.push_back( std::vector< dataformats::layerfracs >() );
  layer_hitvars_.push_back( std::vector< dataformats::layerhitvars >() );
  clusterdep_.push_back( std::vector< dataformats::clustervars >() );
}

Analyzer::~Analyzer()
{
}

void Analyzer::resize_vectors() 
{
  if(this->lmax == 0)
    throw std::runtime_error("Please define the maximum layer before calling resize_vectors().");
  this->layer_fracs_.resize(this->lmax);
  this->layer_hitvars_.resize(this->lmax);
  this->clusterdep_.resize(this->lmax);
}

void Analyzer::runCLUE() {
  float tot_en;
  dataformats::layervars layerdep_vars;
  dataformats::clustervars clusterdep_vars;

  std::vector< std::vector<float> > x_;
  std::vector< std::vector<float> > y_;
  std::vector< std::vector<unsigned int> > layer_;
  std::vector< std::vector<float> > weight_;
  std::vector< std::vector<unsigned int> > rechits_id_; //not required by CLUE
  std::vector< std::vector<float> > impactX_;
  std::vector< std::vector<float> > impactY_;
  
  std::pair<unsigned int, float> out_pair;
  unsigned int nevents = 0;
  float beam_energy = -1;
  CLUEAlgo clueAlgo(dc_, kappa_, ecut_); //non-verbose
  CLUEAnalysis clueAna(this->st_);
  this->lmax = clueAna.getLayerMax();
  resize_vectors();

  for(unsigned int i=0; i<nfiles_; ++i) 
    {
      std::cout << "Processing file number " << i+1 << std::endl;
      //clear data vectors reading a new file
      x_.clear();
      y_.clear();
      layer_.clear();
      weight_.clear();
      rechits_id_.clear();

      out_pair = this->_readTree( this->names_[i].first, x_, y_, layer_, weight_, rechits_id_, impactX_, impactY_ );
      nevents = out_pair.first;
      beam_energy = out_pair.second;
      beam_energies_[i] = beam_energy;

      for(unsigned iEvent=0; iEvent<nevents; ++iEvent)
	{	  
	  //std::cout << "Inside this tree there are " << nevents << " events: ";
	  //std::cout << iEvent/static_cast<float>(nevents)*100 << "% \r";

	  if( x_[iEvent].size() == 0) //empty event
	    continue;

	  //calculate quantities including outliers
	  std::vector<unsigned> tot_hits_per_layer(this->lmax, 0);
	  std::vector<float> tot_en_per_layer(this->lmax, 0.f);
	  for(unsigned int j=0; j<layer_[iEvent].size(); ++j)
	    {
	      unsigned int layeridx = layer_.at(iEvent).at(j) - 1;
	      if(layeridx > this->lmax - 1)
		continue;
	      if( ! ecut_selection(weight_.at(iEvent).at(j), layeridx) )
		continue;
	      tot_hits_per_layer.at(layeridx) += 1;
	      tot_en_per_layer.at(layeridx) += weight_.at(iEvent).at(j);
	    }

	  //run the algorithm per event
	  if ( clueAlgo.setPoints(x_[iEvent].size(), &x_[iEvent][0], &y_[iEvent][0], &layer_[iEvent][0], &weight_[iEvent][0]) )
	    continue; //no event passed the initial energy cut
	  clueAlgo.makeClusters();
	  clueAlgo.infoSeeds();
	  clueAlgo.infoHits();

	  //calculate the total energy that was clusterized (excluding outliers)
	  clueAna.calculateEnergy( clueAlgo.getHitsWeight(), clueAlgo.getHitsClusterId() );
	  tot_en = clueAna.getTotalEnergyOutput("", false); //non-verbose
	  this->en_total_[i].push_back( std::make_tuple( tot_en, beam_energy) ); 
	  //calculate per layer fraction of clusterized number of hits and energy
	  clueAna.calculateLayerDepVars( clueAlgo.getHitsPosX(), clueAlgo.getHitsPosY(),
					 clueAlgo.getHitsWeight(), clueAlgo.getHitsClusterId(), clueAlgo.getHitsLayerId(),
					 clueAlgo.getHitsRho(), clueAlgo.getHitsDistanceToHighest(), clueAlgo.getHitsSeeds(), clueAlgo.getNHitsInCluster());
	  layerdep_vars = clueAna.getTotalLayerDepOutput();

	  //fill fractions (the denominators include outliers!)
	  dataformats::layerfracs fracs_tmp(lmax);
	  dataformats::layerhitvars hitvars_tmp(lmax);
	  for(unsigned int j=0; j<this->lmax; ++j)
	    {
	      if (tot_hits_per_layer[j] != 0 and tot_en_per_layer[j] != 0)
		{
		  std::vector<float> energy_in_this_layer = std::get<1>(layerdep_vars[j]);
		  float energy_sum = std::accumulate( energy_in_this_layer.begin(), energy_in_this_layer.end(), decltype(energy_in_this_layer)::value_type(0) );
		  fracs_tmp[j] = std::make_tuple( static_cast<float>( std::get<0>(layerdep_vars[j]) ) / tot_hits_per_layer[j], energy_sum / tot_en_per_layer[j]);
		  hitvars_tmp[j] = std::make_tuple(energy_in_this_layer, std::get<2>(layerdep_vars[j]), std::get<3>(layerdep_vars[j]), std::get<4>(layerdep_vars[j]), std::get<5>(layerdep_vars[j]), std::get<6>(layerdep_vars[j]), std::get<7>(layerdep_vars[j]) );
		}
	      else
		{
		  std::vector<float> tmpv0(0,-.1f);
		  std::vector<float> tmpv1(0,-.1f);
		  std::vector<float> tmpv2(0,-.1f);
		  std::vector<bool>  tmpv3(0,false);
		  std::vector<float> tmpv4(0,-.1f);
		  std::vector<float> tmpv5(0,-.1f);
		  std::vector<unsigned int> tmpv6(0,0);
		  fracs_tmp[j] = std::make_tuple( -.1f, -.1f );
		  hitvars_tmp[j] = std::make_tuple( std::move(tmpv0), std::move(tmpv1), std::move(tmpv2), std::move(tmpv3), std::move(tmpv4), std::move(tmpv5), std::move(tmpv6) );
		}
	    }
	  this->layer_fracs_.at(i).push_back( fracs_tmp );
	  this->layer_hitvars_.at(i).push_back( hitvars_tmp );

	  //calculate per cluster and per layer clusterized number of hits and energy
	  clueAna.calculateClusterDepVars( clueAlgo.getHitsPosX(), clueAlgo.getHitsPosY(), clueAlgo.getHitsWeight(), clueAlgo.getHitsClusterId(), clueAlgo.getHitsLayerId(), impactX_[iEvent], impactY_[iEvent] );
	  clusterdep_vars = clueAna.getTotalClusterDepOutput();
	  this->clusterdep_.at(i).push_back( clusterdep_vars );
	}
    }
}

std::pair<unsigned int, float> Analyzer::_readTree( const std::string& infile,
			std::vector< std::vector<float> >& x, std::vector< std::vector<float> >& y, 
			std::vector< std::vector<unsigned int> >& layer, std::vector< std::vector<float> >& weight, 
	                std::vector< std::vector<unsigned int> >& rechits_id,
		        std::vector< std::vector<float> >& impactX, std::vector< std::vector<float> >& impactY) {
  //enable parallel execution
  ROOT::EnableImplicitMT( ncpus_ );
  //creates RDataFrame object
  std::string intree = "relevant_branches";
  ROOT::RDataFrame d( intree.c_str(), infile.c_str() );
  //declare data vectors per event to be separately filled by independent cpu threads. dimension: (ncpus, nentries)
  std::array< std::vector< std::vector<float> >, ncpus_ > x_split;
  std::array< std::vector< std::vector<float> >, ncpus_ > y_split;
  std::array< std::vector< std::vector<unsigned int> >, ncpus_ > layer_split;
  std::array< std::vector< std::vector<float> >, ncpus_ > weight_split;
  std::array< std::vector< std::vector<unsigned int> >, ncpus_ > rechits_id_split;
  std::array< std::vector< std::vector<float> >, ncpus_ > impactX_split;
  std::array< std::vector< std::vector<float> >, ncpus_ > impactY_split;

  float beam_energy = 0;
  //lambda function passed to RDataFrame.ForeachSlot(); the first parameters gives the thread number contained in [0;ncpus[
  auto fill = [&x_split, &y_split, &layer_split, &weight_split, &rechits_id_split, &beam_energy, 
	       &impactX_split, &impactY_split](unsigned int slot, std::vector<float>& x_, std::vector<float>& y_, 
					       std::vector<unsigned int>& layer_, std::vector<float>& weight_,
					       std::vector<unsigned int>& rechits_id_, float beamen, 
					       std::vector<float> impactX_, std::vector<float> impactY_) {

    x_split[slot].push_back(x_);
    y_split[slot].push_back(y_);
    layer_split[slot].push_back(layer_);
    weight_split[slot].push_back(weight_);
    rechits_id_split[slot].push_back(rechits_id_);
    impactX_split[slot].push_back(impactX_);
    impactY_split[slot].push_back(impactY_);
    if(beam_energy == 0)
      beam_energy = beamen; //only changes the first time to avoid extra operations
  };

  //loop over the TTree pointed by the RDataFrame
  d.ForeachSlot(fill, {"ce_clean_x", "ce_clean_y", "ce_clean_layer", "ce_clean_energy_MeV", "ce_clean_detid", "beamEnergy",
	"impactX_shifted", "impactY_shifted"});

  //calculate number of events taking into account that ncpus is just a hint to EnableImplicitMT
  std::vector<unsigned int> nevents_v;
  for(unsigned int iThread=0; iThread<ncpus_; ++iThread) {
    nevents_v.push_back( x_split[iThread].size() ); //assumes all the sizes are the same
  }
  unsigned int nevents = std::accumulate(nevents_v.begin(), nevents_v.end(), 0);
  std::vector<unsigned int> nevents_cumsum( nevents_v.size() );
  std::partial_sum(nevents_v.begin(), nevents_v.end(), nevents_cumsum.begin());

  //merge arrays filled by independent cpu threads into the ones passed by reference
  unsigned int iThread_new=0, iThread_old=0, padding=0;
  for(unsigned int iEvent=0; iEvent<nevents; ++iEvent) {

    //get cpu thread corresponding to the event being processed
    unsigned int a = 0;
    for(auto i: nevents_cumsum) {
      if(iEvent<i) iThread_new = a;
      else ++a;
    }

    if(iThread_new != iThread_old) {
      padding += nevents_v[iThread_old];
      assert(padding == iEvent);
    }
    iThread_old = iThread_new;
    x.push_back(x_split[iThread_new][iEvent - padding]);
    y.push_back(y_split[iThread_new][iEvent - padding]);
    layer.push_back(layer_split[iThread_new][iEvent - padding]);
    weight.push_back(weight_split[iThread_new][iEvent - padding]);
    rechits_id.push_back(rechits_id_split[iThread_new][iEvent - padding]);
    impactX.push_back(impactX_split[iThread_new][iEvent - padding]);
    impactY.push_back(impactY_split[iThread_new][iEvent - padding]);
  }
  assert( x.size() == y.size() );
  assert( x.size() == layer.size() );
  assert( x.size() == weight.size() );
  assert( x.size() == rechits_id.size() );
  assert( x.size() == impactX.size() );
  assert( x.size() == impactY.size() );
  return std::make_pair(nevents, beam_energy);
}

bool Analyzer::ecut_selection(const float& energy, const unsigned int& layer)
{
  float endeposited_mip = layer < detectorConstants::layerBoundary ? detectorConstants::energyDepositedByMIP[0] : detectorConstants::energyDepositedByMIP[1];

  float weight_tmp;
  if(layer >= detectorConstants::nlayers_emshowers)
    weight_tmp = detectorConstants::globalWeightCEH;
  else
    weight_tmp = detectorConstants::dEdX.at(layer);
  
  return energy > ecut_ * detectorConstants::sigmaNoiseSiSensor / endeposited_mip * weight_tmp;
}

void Analyzer::sum_energy(const bool& with_ecut)
{
  std::mutex mut; //anonymous function to pass to RDataFrame.ForEach()

  ROOT::EnableImplicitMT( ncpus_ ); //enable parallelism

  for(unsigned int i=0; i<nfiles_; ++i)
    {
      auto sum_ce = [&](const std::vector<float>& ce_en, const std::vector<unsigned int>& ce_layer, float beamen)
		    { 
		      float entot = 0.f;
		      for(unsigned int ien=0; ien<ce_en.size(); ++ien)
			{
			  if(ce_layer[ien] > this->lmax)
			    continue;
			  if(with_ecut and ! ecut_selection(ce_en[ien], ce_layer[ien]-1))
			    continue;

			  entot += ce_en[ien];
			}
		      { //mutex lock scope
			std::lock_guard lock(mut);
			this->en_total_[i].push_back(std::make_tuple(entot, beamen));
		      }
		    };

      auto sum_ahc = [&](const std::vector<float>& ce_en, const std::vector<unsigned int>& ce_layer, const std::vector<float>& ahc_en, float beamen)
		     { 
		       float entot = 0.f;
		       for(unsigned int ien=0; ien<ce_en.size(); ++ien)
			 {
			   if(ce_layer[ien] > this->lmax)
			     continue;
			   if(with_ecut and ! ecut_selection(ce_en[ien], ce_layer[ien]-1))
			     continue;

			   entot += ce_en[ien];
			 }
		       for(unsigned int ien=0; ien<ahc_en.size(); ++ien)
			 entot += ahc_en[ien];

		       { //mutex lock scope
			 std::lock_guard lock(mut);
			 this->en_total_[i].push_back(std::make_tuple(entot, beamen));
		       }
		     };

      //define dataframe that owns the TTree
      ROOT::RDataFrame d(this->names_[i].second.c_str(), this->names_[i].first.c_str());
      //store the contents of the TTree according to the specified columns
      en_total_[i].clear();
      if(this->st_ == SHOWERTYPE::EM)
	d.Foreach(sum_ce, {"ce_clean_energy_MeV", "ce_clean_layer", "beamEnergy"});
      else if(this->st_ == SHOWERTYPE::HAD)
	d.Foreach(sum_ahc, {"ce_clean_energy_MeV", "ce_clean_layer", "ahc_clean_energy_MeV", "beamEnergy"});
    }
};

int Analyzer::sanity_checks(const std::string& fname) 
{
  if(fname.substr(fname.length()-5) != ".root")
    throw std::invalid_argument("The file has to be i the ROOT format.");
  return 1;
}

void Analyzer::save_to_file(const std::string& filename) {
  std::ofstream oFile(filename);
  std::cout << "SAVE: " << filename << std::endl;
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      std::string curr_name = std::get<0>(names_[i]);
      curr_name = curr_name.substr(curr_name.length()-10, 5);

      oFile << "ensum" << curr_name << ",";
      oFile << "beamen" << curr_name;
      if(i<nfiles_-1)
	oFile << ",";
    }
  oFile << std::endl;

  //get size of larger energy vector
  std::vector<unsigned int> en_sizes(nfiles_, 0);
  for(unsigned int i=0; i<nfiles_; ++i)
    en_sizes[i] = en_total_[i].size();
  const unsigned int max = *( std::max_element(en_sizes.begin(), en_sizes.end()) );
  
  for(unsigned int k=0; k<max; ++k)
    {
      for(unsigned int i=0; i<nfiles_; ++i)
	{
	  if(k<en_total_[i].size())
	    {
	      oFile << std::to_string( std::get<0>(en_total_[i][k]) ) << ",";
	      oFile << std::to_string( std::get<1>(en_total_[i][k]) );
	    }
	  else
	    oFile << "-99., -99.";
	  if(i<nfiles_-1)
	    oFile << ",";
	}
      oFile << std::endl;
    }
}

void Analyzer::save_to_file_layer_dependent(const std::string& filename) {
  std::cout << "SAVE: " << filename << std::endl;
  std::cout << "NFILES: " << nfiles_ << std::endl;
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      //create TTree
      const std::string str_start = filename.substr(0,filename.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
      const std::string str_end = filename.substr(filename.find('.', 20), 5);
      const std::string filename_with_energy = str_start + "_beamen" + std::to_string( static_cast<int>(beam_energies_.at(i)) ) + str_end;

      TFile file( filename_with_energy.c_str(), "RECREATE");
      std::string name = "tree" + std::to_string(i);
      TTree tmptree( name.c_str(), name.c_str());

      //define branches
      tmptree.Branch("BeamEnergy", &beam_energies_.at(i));
      std::vector< float > fracs_hits(this->lmax);
      std::vector< float > fracs_en(this->lmax);
      std::vector< std::vector<float> > energies(this->lmax);
      std::vector< std::vector<float> > posx(this->lmax);
      std::vector< std::vector<float> > posy(this->lmax);
      std::vector< std::vector<float> > rhos(this->lmax);
      std::vector< std::vector<float> > deltas(this->lmax);
      std::vector< std::vector<bool> > seeds(this->lmax);
      std::vector< std::vector<unsigned> > clustersizes(this->lmax);
      for(unsigned int ilayer=0; ilayer<this->lmax; ++ilayer) 
	{
	  std::string bname_hits     = "NhitsFrac_layer"  + std::to_string(ilayer + 1);
	  std::string bname_energy   = "EnergyFrac_layer" + std::to_string(ilayer + 1);
	  std::string bname_energies = "Energies_layer"   + std::to_string(ilayer + 1);
	  std::string bname_posx     = "PosX_layer"       + std::to_string(ilayer + 1);
	  std::string bname_posy     = "PosY_layer"       + std::to_string(ilayer + 1);
	  std::string bname_rho      = "Densities_layer"  + std::to_string(ilayer + 1);
	  std::string bname_delta    = "Distances_layer"  + std::to_string(ilayer + 1);
	  std::string bname_seeds    = "isSeed_layer"     + std::to_string(ilayer + 1);
	  std::string bname_clsize   = "ClustSize_layer"  + std::to_string(ilayer + 1);
	  tmptree.Branch(bname_hits.c_str(),     &fracs_hits[ilayer]);
	  tmptree.Branch(bname_energy.c_str(),   &fracs_en[ilayer]);
	  tmptree.Branch(bname_energies.c_str(), &energies[ilayer]);
	  tmptree.Branch(bname_posx.c_str(),     &posx[ilayer]);
	  tmptree.Branch(bname_posy.c_str(),     &posy[ilayer]);
	  tmptree.Branch(bname_rho.c_str(),      &rhos[ilayer]);
	  tmptree.Branch(bname_delta.c_str(),    &deltas[ilayer]);
	  tmptree.Branch(bname_seeds.c_str(),    &seeds[ilayer]);
	  tmptree.Branch(bname_clsize.c_str(),   &clustersizes[ilayer]);
	}

      //loop over TTree and fill branches
      assert( this->layer_fracs_.at(i).size() == this->layer_hitvars_.at(i).size() );
      unsigned int nentries = this->layer_fracs_.at(i).size(); // read the number of entries in the t3
      for (unsigned int ientry = 0; ientry<nentries; ++ientry) 
	{
	  for(unsigned int ilayer=0; ilayer<this->lmax; ++ilayer) 
	    {
	      fracs_hits[ilayer]   = std::get<0>( this->layer_fracs_.at(i).at(ientry).at(ilayer) );
	      fracs_en[ilayer]     = std::get<1>( this->layer_fracs_.at(i).at(ientry).at(ilayer) );
	      energies[ilayer]     = std::get<0>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      rhos[ilayer]         = std::get<1>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      deltas[ilayer]       = std::get<2>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      seeds[ilayer]        = std::get<3>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      posx[ilayer]         = std::get<4>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      posy[ilayer]         = std::get<5>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	      clustersizes[ilayer] = std::get<6>( this->layer_hitvars_.at(i).at(ientry).at(ilayer) );
	    }
	  tmptree.Fill();
	}
      file.Write();
      file.Close();
    }
}

void Analyzer::save_to_file_cluster_dependent(const std::string& filename) {
  std::cout << std::endl;
  std::cout << "SAVE: " << filename << std::endl;
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      //create TTree
      std::string name = "tree" + std::to_string(i);
      TFile file( filename.c_str(), "RECREATE");
      TTree tmptree( name.c_str(), name.c_str());

      //define branches
      tmptree.Branch("BeamEnergy", &beam_energies_.at(i));
      std::vector< std::vector<unsigned> > arr_hits(this->lmax);
      std::vector< std::vector<float> > arr_en(this->lmax);
      std::vector< std::vector<float> > arr_x(this->lmax);
      std::vector< std::vector<float> > arr_y(this->lmax);
      std::vector< std::vector<float> > arr_dx(this->lmax);
      std::vector< std::vector<float> > arr_dy(this->lmax);
      for(unsigned int ilayer=0; ilayer<this->lmax; ++ilayer) 
	{
	  std::string bname_hits   = "Nhits_layer"  + std::to_string(ilayer + 1);
	  std::string bname_energy = "Energy_layer" + std::to_string(ilayer + 1);
	  std::string bname_x      = "X_layer"      + std::to_string(ilayer + 1);
	  std::string bname_y      = "Y_layer"      + std::to_string(ilayer + 1);
	  std::string bname_dx     = "dX_layer"     + std::to_string(ilayer + 1);
	  std::string bname_dy     = "dY_layer"     + std::to_string(ilayer + 1);
	  tmptree.Branch(bname_hits.c_str(),   &arr_hits[ilayer]);
	  tmptree.Branch(bname_energy.c_str(), &arr_en[ilayer]);
	  tmptree.Branch(bname_x.c_str(),      &arr_x[ilayer]);
	  tmptree.Branch(bname_y.c_str(),      &arr_y[ilayer]);	 
	  tmptree.Branch(bname_dx.c_str(),     &arr_dx[ilayer]);
	  tmptree.Branch(bname_dy.c_str(),     &arr_dy[ilayer]);	 
	}

      //loop over TTree and fill branches
      unsigned int nentries = this->clusterdep_.at(i).size(); // read the number of entries in the t3
      for (unsigned int ientry = 0; ientry<nentries; ++ientry) 
	{
	  for(unsigned int ilayer=0; ilayer<this->lmax; ++ilayer) 
	    {
	      arr_hits[ilayer] = std::get<0>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	      arr_en[ilayer]   = std::get<1>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	      arr_x[ilayer]    = std::get<2>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	      arr_y[ilayer]    = std::get<3>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	      arr_dx[ilayer]   = std::get<4>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	      arr_dy[ilayer]   = std::get<5>( this->clusterdep_.at(i).at(ientry).at(ilayer) );
	    }
	  tmptree.Fill();
	}
      //tmptree.Write("", TObject::kOverwrite); //save only the new version of the tree
      file.Write();
      file.Close();
    }

}
