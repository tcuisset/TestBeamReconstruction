#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector(const std::string& in_file_path, const std::string& out_file_path,
		   const std::string& datatype_, const std::string& showertype_, const int& beam_energy_,
		   std::optional<std::string> in_tree_name, std::optional<std::string> out_tree_name)
{
  sanity_checks(in_file_path);
  this->indata_.file_path = in_file_path;
  this->outdata_.file_path = out_file_path;
  this->beam_energy = beam_energy_;

  if(datatype_.find("sim") != std::string::npos)
    this->datatype = DATATYPE::DATA;
  else if(datatype_.find("data") != std::string::npos)
    this->datatype = DATATYPE::MC;
  else
    throw std::invalid_argument("Wrong data type: " + datatype_);
    
  if(showertype_ == "em")
    this->showertype = SHOWERTYPE::EM;
  else if(showertype_ == "had")
    this->showertype = SHOWERTYPE::HAD;
  
  if(in_tree_name != std::nullopt)
    this->indata_.tree_name = in_tree_name.value();
  if(out_tree_name != std::nullopt)
    this->indata_.tree_name = out_tree_name.value();
}

Selector::~Selector()
{
}

bool Selector::common_selection(const unsigned& layer, const float& energy, const unsigned& chip, const unsigned& channel)
{
  if(layer<1 || layer>detectorConstants::totalnlayers)
    throw std::out_of_range("Unphysical layer number: " + std::to_string(layer));
  else if(this->showertype == SHOWERTYPE::EM and layer>detectorConstants::nlayers_emshowers && layer<=detectorConstants::totalnlayers)
    return false;
  else if(this->showertype == SHOWERTYPE::HAD and layer>36 && layer<=37)
    return false;

  if(energy<=0.5) //MIP noise cut
    return false;
  if( (chip==0 and layer==1) or (chip==3 and channel==22) ) //mask noisy chip
    return false;
  return true;
}

void Selector::select_relevant_branches()
{
  ROOT::EnableImplicitMT( ncpus_ ); //enable parallelism
  
  auto clean_float_arrays = [&](std::vector<float> var, std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel) {
	unsigned layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<float> var_clean;
	var_clean.reserve(ensize);
	
	for(unsigned i=0; i<ensize; ++i) {
	  energy = en[i];
	  layer = l[i];
	  if( common_selection(layer, energy, chip[i], channel[i]) )
	    var_clean.push_back( var[i] );
	}
	return var_clean;  
      };

  auto clean_detids = [&](std::vector<unsigned> detid, std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel) {
	unsigned layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<unsigned> detid_clean;
	detid_clean.reserve(ensize);
	
	for(unsigned i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];
	    if( common_selection(layer, energy, chip[i], channel[i]) )
	      detid_clean.push_back( detid[i] );
	  }
	return detid_clean;  
      };

  auto clean_layers = [&](std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel) {
	unsigned layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<unsigned> layer_clean;
	layer_clean.reserve(ensize);

	for(unsigned i=0; i<ensize; ++i) {
	  energy = en[i];
	  layer = l[i];
	  if( common_selection(layer, energy, chip[i], channel[i]) )
	    layer_clean.push_back( layer );
	}
	return layer_clean;  
      };
      
  auto weight_and_clean_energy = [&](std::vector<float> en, std::vector<unsigned> l, std::vector<unsigned> chip, std::vector<unsigned> channel) {
	unsigned layer = 0;
	float energy = 0.f;
	float weight = 1.f;
	size_t ensize = en.size();
	std::vector<float> en_weighted;
	en_weighted.reserve(ensize); //maximum possible size

	for(unsigned i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];

	    if( (this->showertype==SHOWERTYPE::EM and layer>0 and layer<=detectorConstants::nlayers_emshowers)
		or (this->showertype==SHOWERTYPE::HAD and layer>0 and layer<=detectorConstants::totalnlayers) )
	      {
		if(this->showertype==SHOWERTYPE::HAD and layer>detectorConstants::nlayers_emshowers) //remove as soon as the extra dEdX are known
		  weight = 1.;
		else
		  {
		    //remove this bit once the FH weights are established
		    float XXXXweight;
		    if(layer-1 > detectorConstants::nlayers_emshowers)
		      XXXXweight = 1.f;
		    else
		      XXXXweight = detectorConstants::dEdX.at(layer-1);
		    ///////////////////////////////////////////////////

		    weight = XXXXweight;
		  }
	      }
	    else
	      continue;

	    if( common_selection(layer, energy, chip[i], channel[i]) )
	      en_weighted.push_back( energy * weight );
	  }
	return en_weighted;  
      };

  //filters events with pion contamination in the FH region
  auto limit_max_hits_fh = [&](const std::vector<unsigned> l) {
    unsigned counter=0;
    for(unsigned i=0; i<l.size(); ++i)
      {
	if( (l[i]>detectorConstants::nlayers_emshowers and l[i]<36) or (l[i]>37 and l[i]<=detectorConstants::totalnlayers) )
	  ++counter;
      }

    if(counter<50)
      return true;
    return false;
  };

  auto guarantee_95_containment = [&](const std::vector<float> en) {
    if( std::accumulate( en.begin(), en.end(), 0.) < 0.95*this->beam_energy*1000 )
      return false;
    return true;
  };
  
  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  auto partial_process = d.Filter(limit_max_hits_fh, {"rechit_layer"})
    .Define(this->newcol_clean_detid_,     clean_detids,              {"rechit_detid", "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel"})
    .Define(this->newcol_clean_x_,         clean_float_arrays,        {"rechit_x",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel" })
    .Define(this->newcol_clean_y_,         clean_float_arrays,        {"rechit_y",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel"})
    .Define(this->newcol_clean_z_,         clean_float_arrays,        {"rechit_z",     "rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel"})
    .Define(this->newcol_clean_layer_,     clean_layers,              {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel"})
    .Define(this->newcol_clean_energy_,    weight_and_clean_energy,   {"rechit_energy", "rechit_layer", "rechit_chip", "rechit_channel"})
    .Filter(guarantee_95_containment, {"rechit_energy"});

    if(this->datatype == DATATYPE::MC)
      partial_process.Filter("ahc_energySum==0").Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), this->cols_);
    else
      partial_process.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), this->cols_);
};

void Selector::print_relevant_branches(const int& nrows=5, std::optional<std::string> filename)
{
  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  //prints TTree content according to the specified columns; a multi-threaded version is not supported
  auto display = d.Display(cols_, nrows);
  //write contents to file
  std::ofstream f;
  if(filename != std::nullopt)
    f.open(filename.value());
  else
    f.open("data.txt");
  f << display->AsString();
  f.close();
};

int Selector::sanity_checks(const std::string& fname) 
{
  if(fname.substr(fname.length()-5) != ".root")
    throw std::invalid_argument("The file has to be i the ROOT format.");
  return 1;
}
