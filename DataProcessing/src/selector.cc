#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector(const std::string& in_file_path, std::optional<std::string> out_file_path,
		   std::optional<std::string> in_tree_name, std::optional<std::string> out_tree_name)
{
  sanity_checks(in_file_path);
  this->indata_.file_path = in_file_path;
  if(out_file_path != std::nullopt)
    {
      this->outdata_.file_path = out_file_path.value();
      sanity_checks(this->outdata_.file_path);
    }
  if(in_tree_name != std::nullopt)
    this->indata_.tree_name = in_tree_name.value();
  if(out_tree_name != std::nullopt)
    this->indata_.tree_name = out_tree_name.value();
}

Selector::~Selector()
{
}

void Selector::select_relevant_branches()
{
  if(this->outdata_.file_path == "/eos/user/b/bfontana/TestBeamReconstruction/default_output.txt")
    throw std::invalid_argument("In order to perform TTree skimming you must specify an output file path.");
  //enable parallelism
  ROOT::EnableImplicitMT( ncpus_ );

  auto clean_float_arrays = [&](std::vector<float> var, std::vector<float> en, std::vector<unsigned int> l, std::vector<unsigned int> chip) {
	unsigned int layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<float> var_clean;
	var_clean.reserve(ensize);
	
	for(unsigned int i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];
	    if(layer>detectorConstants::nlayers_emshowers && layer<=detectorConstants::totalnlayers)
	      continue;
	    else if(layer<1 || layer>detectorConstants::totalnlayers)
	      throw std::out_of_range("Unphysical layer number: "+std::to_string(layer));
	    if(energy<=0.5) //MIP noise cut
	      continue;
	    if(chip[i]==0 and layer==1) //mask noisy chip
	      continue;
	    var_clean.push_back( var[i] );
	  }
	return var_clean;  
      };

  auto clean_detids = [&](std::vector<unsigned int> detid, std::vector<float> en, std::vector<unsigned int> l, std::vector<unsigned int> chip) {
	unsigned int layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<unsigned int> detid_clean;
	detid_clean.reserve(ensize);
	
	for(unsigned int i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];
	    if(layer>detectorConstants::nlayers_emshowers && layer<=detectorConstants::totalnlayers)
	      continue;
	    else if(layer<1 || layer>detectorConstants::totalnlayers)
	      throw std::out_of_range("Unphysical layer number: "+std::to_string(layer));
	    if(energy<=0.5) //MIP noise cut
	      continue;
	    if(chip[i]==0 and layer==1) //mask noisy chip
	      continue;
	    detid_clean.push_back( detid[i] );
	  }
	return detid_clean;  
      };

  auto clean_layers = [&](std::vector<float> en, std::vector<unsigned int> l, std::vector<unsigned int> chip) {
	unsigned int layer = 0;
	float energy = 0.f;
	size_t ensize = en.size();
	std::vector<unsigned int> layer_clean;
	layer_clean.reserve(ensize);

	for(unsigned int i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];
	    if(layer>detectorConstants::nlayers_emshowers and layer<=detectorConstants::totalnlayers)
	      continue;
	    else if(layer<1 or layer>detectorConstants::totalnlayers)
	      {
		throw std::out_of_range("Unphysical layer number: "+std::to_string(layer));
		continue;
	      }
	    if(energy<=0.5) //MIP noise cut
	      continue;
	    if(chip[i]==0 and layer==1) //mask noisy chip
	      continue;

	    layer_clean.push_back( layer );
	  }
	return layer_clean;  
      };
      
  auto weight_and_clean_energy = [&](std::vector<float> en, std::vector<unsigned int> l, std::vector<unsigned int> chip) {
	unsigned int layer = 0;
	float energy = 0.f;
	float weight = 1.f;
	size_t ensize = en.size();
	std::vector<float> en_weighted;
	en_weighted.reserve(ensize); //maximum possible size

	for(unsigned int i=0; i<ensize; ++i)
	  {
	    energy = en[i];
	    layer = l[i];
	    if(energy<=0.5)  //MIP noise cut
	      continue;
	    if(chip[i]==0 and layer==1) //mask noisy chip
	      continue;
	    
	    if(layer>0 and layer<=detectorConstants::nlayers_emshowers)
	      weight = detectorConstants::dEdX.at(layer-1);
	    else
	      continue;
	    en_weighted.push_back( energy * weight );
	  }
	return en_weighted;  
      };
 

  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  //convert from MIPs to MeV
  d.Define(this->newcol_clean_detid_,      clean_detids,              {"rechit_detid", "rechit_energy", "rechit_layer", "rechit_chip"})
    .Define(this->newcol_clean_x_,         clean_float_arrays,        {"rechit_x",     "rechit_energy", "rechit_layer", "rechit_chip"})
    .Define(this->newcol_clean_y_,         clean_float_arrays,        {"rechit_y",     "rechit_energy", "rechit_layer", "rechit_chip"})
    .Define(this->newcol_clean_z_,         clean_float_arrays,        {"rechit_z",     "rechit_energy", "rechit_layer", "rechit_chip"})
    .Define(this->newcol_clean_layer_,     clean_layers,              {"rechit_energy", "rechit_layer", "rechit_chip"})
    .Define(this->newcol_clean_energy_,    weight_and_clean_energy,   {"rechit_energy", "rechit_layer", "rechit_chip"})
    .Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), this->cols_);
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
