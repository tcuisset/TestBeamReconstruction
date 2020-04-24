#include "UserCode/DataProcessing/interface/analyzer.h"

Analyzer::Analyzer(const std::vector< std::string >& in_file_path, const std::string& in_tree_name)
{
  nfiles_ = in_file_path.size();
  std::cout << "Number of files being processed: " << nfiles_ << std::endl;
  std::cout << "Files: " << std::endl;
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      sanity_checks(in_file_path[i]);
      std::cout << "#" << std::to_string(i+1) << " " << in_file_path[i] << std::endl;
      names_.push_back( std::make_pair(in_file_path[i], in_tree_name) );
      en_total_.push_back( std::vector< std::tuple<float, float> >() );
      fracs_.push_back( std::vector< std::array< std::tuple<float, float>, 28 > >() );
    }
}

//Overloaded constructor for job submission. Each job processes one file only.
Analyzer::Analyzer(const std::string& in_file_path, const std::string& in_tree_name)
{
  nfiles_ = 1;
  std::cout << "Number of files being processed: " << nfiles_ << std::endl;
  std::cout << "File: " << in_file_path << std::endl;
  sanity_checks(in_file_path);
  names_.push_back( std::make_pair(in_file_path, in_tree_name) );
  en_total_.push_back( std::vector< std::tuple<float, float> >() );
  fracs_.push_back( std::vector< std::array< std::tuple<float, float>, nlayers_ > >() );
}

Analyzer::~Analyzer()
{
}

void Analyzer::runCLUE(float dc, float rhoc_300, float rhoc_200) {
  float tot_en = 0;
  std::array< std::tuple<unsigned int, float>, nlayers_> layerdep_vars;

  std::vector< std::vector<float> > x_;
  std::vector< std::vector<float> > y_;
  std::vector< std::vector<unsigned int> > layer_;
  std::vector< std::vector<float> > weight_;
  std::vector< std::vector<unsigned int> > rechits_id_; //not required by CLUE
  
  std::pair<unsigned int, float> out_pair;
  unsigned int nevents = 0;
  float beam_energy = -1;
  CLUEAlgo clueAlgo(dc, rhoc_300, rhoc_200); //non-verbose
  CLUEAnalysis clueAna;
  for(unsigned int i=0; i<nfiles_; ++i) 
    {
      std::cout << "Progress: " << i/static_cast<float>(nfiles_)*100 << "%" << std::endl;
      //clear data vectors reading a new file
      x_.clear();
      y_.clear();
      layer_.clear();
      weight_.clear();
      rechits_id_.clear();

      out_pair = this->_readTree( this->names_[i].first, x_, y_, layer_, weight_, rechits_id_ );
      nevents = out_pair.first;
      beam_energy = out_pair.second;

      for(unsigned int iEvent=0; iEvent<nevents; ++iEvent) //otherwise error for unused 'nevents' variable
	{
	  std::cout << "Inside this tree there are " << nevents << " events: ";
	  std::cout << iEvent/static_cast<float>(nevents)*100 << "% \r";

	  //run the algorithm per event
	  clueAlgo.setPoints(x_[iEvent].size(), &x_[iEvent][0], &y_[iEvent][0], &layer_[iEvent][0], &weight_[iEvent][0]);
	  clueAlgo.makeClusters();

	  //calculate quantities including outliers
	  std::array<unsigned int, nlayers_> tot_hits_per_layer = {{0}};
	  std::array<float, nlayers_> tot_en_per_layer = {{0}};
	  for(unsigned int j=0; j<layer_[iEvent].size(); ++j)
	    {
	      unsigned int layeridx = layer_[iEvent][j] - 1;
	      tot_hits_per_layer[layeridx] += 1;
	      tot_en_per_layer[layeridx] += weight_[iEvent][j];
	    }
	  
	  //calculate the total energy that was clusterized (excluding outliers)
	  clueAna.calculatePositionsAndEnergy( clueAlgo.getHitsClusterX(), clueAlgo.getHitsClusterY(), clueAlgo.getHitsWeight(), clueAlgo.getHitsClusterId(), clueAlgo.getHitsLayerId() );
	  tot_en = clueAna.getTotalClusterEnergyOutput("", false); //non-verbose
	  this->en_total_[i].push_back( std::make_tuple(tot_en, beam_energy) ); 
	  clueAna.calculateLayerDepVars( clueAlgo.getHitsWeight(), clueAlgo.getHitsClusterId(), clueAlgo.getHitsLayerId() );
	  layerdep_vars = clueAna.getTotalClusterLayerDepOutput(); //non-verbose
	  //fill fractions (the denominators include outliers!)
	  std::array< std::tuple<float, float>, nlayers_> eventarray_tmp;
	  for(unsigned int j=0; j<nlayers_; ++j)
	    {
	      std::cout << std::endl;
	      std::cout << std::get<0>(layerdep_vars[j]) << ", " << tot_hits_per_layer[j] << std::endl;
	      std::cout << std::get<1>(layerdep_vars[j]) << ", " << tot_en_per_layer[j] << std::endl;
	      if (tot_hits_per_layer[j] != 0 and tot_en_per_layer[j] != 0)
		eventarray_tmp[j] = std::make_tuple(static_cast<float>(std::get<0>(layerdep_vars[j]))/tot_hits_per_layer[j], std::get<1>(layerdep_vars[j])/tot_en_per_layer[j] );
	      else
		eventarray_tmp[j] = std::make_tuple(0., 0.);
	    }
	  this->fracs_.at(i).push_back( eventarray_tmp );
	}
    }
}

std::pair<unsigned int, float> Analyzer::_readTree( const std::string& infile,
			std::vector< std::vector<float> >& x, std::vector< std::vector<float> >& y, 
			std::vector< std::vector<unsigned int> >& layer, std::vector< std::vector<float> >& weight, 
			std::vector< std::vector<unsigned int> >& rechits_id ) {
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

  float beam_energy = 0;
  //lambda function passed to RDataFrame.ForeachSlot(); the first parameters gives the thread number contained in [0;ncpus[
  auto fill = [&x_split, &y_split, &layer_split, &weight_split, &rechits_id_split, &beam_energy](unsigned int slot, 
			  std::vector<float>& x_, std::vector<float>& y_, std::vector<unsigned int>& layer_, std::vector<float>& weight_,
                          std::vector<unsigned int>& rechits_id_, float beamen) {
    x_split[slot].push_back(x_);
    y_split[slot].push_back(y_);
    layer_split[slot].push_back(layer_);
    weight_split[slot].push_back(weight_);
    rechits_id_split[slot].push_back(rechits_id_);
    if(beam_energy == 0)
      beam_energy = beamen; //only changes the first time to avoid extra operations
  };

  //loop over the TTree pointed by the RDataFrame
  d.ForeachSlot(fill, {"rechit_x", "rechit_y", "rechit_layer", "rechit_weighted_energy_MeV", "rechit_detid", "beamEnergy"});

  //calculate number of events taking into account that ncpus is just a hint to EnableImplicitMT
  std::vector<unsigned int> nevents_v;
  for(unsigned int iThread=0; iThread<ncpus_; ++iThread) {
    nevents_v.push_back( x_split[iThread].size() ); //assumes x, y, layer and weight sizes are equal
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
      
  }
  return std::make_pair(nevents, beam_energy);
}

void Analyzer::sum_energy()
{
  //anonymous function to pass to RDataFrame.ForEach()
  std::mutex mut;

  //enable parallelism
  ROOT::EnableImplicitMT( ncpus_ );

  for(unsigned int i=0; i<nfiles_; ++i)
    {

      auto sum = [&](const std::vector<float>& en, float beamen)
	{
	  float entot = std::accumulate(en.begin(), en.end(), 0.);
	  { //mutex lock scope
	    std::lock_guard lock(mut);
	    this->en_total_[i].push_back(std::make_tuple(entot, beamen));
	  }
	};

      //define dataframe that owns the TTree
      ROOT::RDataFrame d(this->names_[i].second.c_str(), this->names_[i].first.c_str());
      //store the contents of the TTree according to the specified columns
      en_total_[i].clear();
      d.Foreach(sum, {"rechit_weighted_energy_MeV", "beamEnergy"});
    }
};

void Analyzer::histogram_checks()
{
  //enable parallelism
  ROOT::EnableImplicitMT( ncpus_ );

  auto sum_energies = [&](std::vector<float> en, std::vector<unsigned int> l) {
    std::vector<float> ensum(en.size(), 0.f);
    unsigned int layer = 0;
    //float thick_corr = 0.;
    float weight = 1.;
    for(unsigned int i=0; i<en.size(); ++i)
      {
	layer = l[i];
	if(layer>0 && layer<27)
	  {
	    //thick_corr = this->thickness_correction_[0];
	    weight = this->energy_weights_[layer];
	  }
	else if(layer >= 27 && layer<29)
	  {
	    //thick_corr = this->thickness_correction_[1];
	    weight = this->energy_weights_[layer];
	  }
	else if(layer >= 29 && layer<=50)
	  {
	    //thick_corr = 0; //ignore hits in the hadronic section
	    weight = 0;
	  }
	else
	  throw std::out_of_range("Unphysical layer number: "+std::to_string(layer));
	ensum.at(i) =  en[i] * /*thick_corr **/ weight;
      }
    return std::accumulate(ensum.begin(), ensum.end(), 0.);
  };
  
  //define dataframe that owns the TTree
  assert(this->names_.size() == 1);
  ROOT::RDataFrame d(this->names_[0].second.c_str(), this->names_[0].first.c_str());
  auto h1 = d.Histo2D<std::vector<unsigned int>, std::vector<float>>({"h1", "RecHit Energy vs. RecHit Layer", 41u, 0, 42, 200u, 0., 14.}, "rechit_layer", "rechit_energy");
  auto d_def = d.Define("energy_sum", sum_energies, {"rechit_energy_MeV", "rechit_layer"}); 
  auto h_sum = d_def.Histo1D("energy_sum");

  TCanvas *c = new TCanvas("c", "c", 1600, 600);
  c->Divide(2,1);
  c->cd(1);
  h_sum->SetLineColor(4);
  h_sum->SetTitle("RecHit Energy Sum per Event");
  h_sum->SetXTitle("Energy Sum [MeV]");
  h_sum->SetYTitle("Counts");
  h_sum->Draw();
  c->cd(2);
  h1->SetLineColor(4);
  h1->SetYTitle("Rechit energy [MIP]");
  h1->SetXTitle("Layer");
  h1->Draw("colz");
  c->SaveAs("histo.png");
  delete c;
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
      curr_name = curr_name.substr(curr_name.length()-8, 3);
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
  std::ofstream oFile(filename);
  std::cout << "SAVE: " << filename << std::endl;
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      std::string curr_name = std::get<0>(names_[i]);
      curr_name = curr_name.substr(curr_name.length()-8, 3);
      for(unsigned int ilayer=0; ilayer<28; ++ilayer)
	{
	  oFile << "nhitsfrac_layer" << ilayer << "_" << curr_name << ",";
	  oFile << "enfrac_" << ilayer << "_" << curr_name << ",";
	}
      if(i<nfiles_-1)
	oFile << ",";
    }
  oFile << std::endl;

  //get size of larger energy vector
  std::vector<unsigned int> fracs_sizes(nfiles_, 0);
  for(unsigned int i=0; i<nfiles_; ++i)
    {
      fracs_sizes[i] = fracs_[i].size();
      assert(fracs_sizes[i] == this->en_total_[i].size());
    }
  const unsigned int max = *( std::max_element(fracs_sizes.begin(), fracs_sizes.end()) );
  
  for(unsigned int k=0; k<max; ++k)
    {
      for(unsigned int i=0; i<nfiles_; ++i)
	{
	  if(k<fracs_sizes[i])
	    {
	      for(unsigned int ilayer=0; ilayer<28; ++ilayer)
		{
		  oFile << std::to_string( std::get<0>(fracs_[i][k][ilayer]) ) << ",";
		  oFile << std::to_string( std::get<1>(fracs_[i][k][ilayer]) );
		}
	    }
	  else
	    oFile << "-99., -99.";
	  if(i<nfiles_-1)
	    oFile << ",";
	}
      oFile << std::endl;
    }
}
