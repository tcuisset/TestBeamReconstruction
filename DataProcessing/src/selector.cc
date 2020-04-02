#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector(const std::string& in_file_path, std::optional<std::string> out_file_path,
		   std::optional<std::string> in_tree_name, std::optional<std::string> out_tree_name)
{
  sanity_checks(in_file_path);
  this->indata_.file_path = in_file_path;
  if(in_tree_name != std::nullopt)
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

  auto convert_energy = [&](std::vector<float> en, std::vector<unsigned int> l) {
    std::vector<float> en_conv(en.size(), 0.);
    float thick_corr = 0;
    for(unsigned int i=0; i<en.size(); ++i)
      {
	if(l[i]>0 && l[i]<27)
	  thick_corr = this->thickness_correction_[0];
	else if(l[i]<29)
	  thick_corr = this->thickness_correction_[1];
	else
	  throw std::out_of_range("The layer number is too large.");
	en_conv.at(i) = en[i] * thick_corr;
      }
    return en_conv;
  };
  
  auto weight_energy = [&](std::vector<float> en, std::vector<unsigned int> l) {
	unsigned int layer = 0;
	std::vector<float> en_conv(en.size(), 0.);
	for(unsigned int i=0; i<en.size(); ++i)
	  {
	    layer = l[i];
	    if(layer > 28 || layer < 1)
	      throw std::out_of_range("The layer number is too large.");
	    en_conv.at(i) = en[i] * this->energy_weights_[layer];
	  }
	return en_conv;  
      };
      
  //define dataframe that owns the TTree
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  //convert from MIPs to MeV
  auto d_def = d.Define(newcol1_, convert_energy, {"rechit_energy", "rechit_layer"})
      .Define(newcol2_, weight_energy, {"rechit_energy_MeV", "rechit_layer"});
  //store the contents of the TTree according to the specified columns
  d.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), cols_);
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
