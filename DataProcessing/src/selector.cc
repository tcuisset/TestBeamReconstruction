#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector(const std::string& in_file_path, const std::string& out_file_path,
		   std::optional<std::string> in_tree_name, std::optional<std::string> out_tree_name)
{
  sanity_checks(in_file_path);
  sanity_checks(out_file_path);
  this->indata_.file_path = in_file_path;
  this->outdata_.file_path = out_file_path;
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
  const int ncpus = 4;
  ROOT::EnableImplicitMT( ncpus );
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  std::cout << this->indata_.tree_name.c_str() << ", " << this->indata_.file_path.c_str() << std::endl;
  std::cout << this->outdata_.tree_name.c_str() << ", " << this->outdata_.file_path.c_str() << std::endl;
  d.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), {"rechit_x", "rechit_y", "rechit_z", "rechit_layer", "rechit_iu", "rechit_iv", "rechit_iU", "rechit_iV", "rechit_type", "rechit_energy", "beamEnergy", "pdgID"});
};

int Selector::sanity_checks(const std::string& fname) 
{
  if(fname.substr(fname.length()-5) != ".root")
    throw std::invalid_argument("The file has to be i the ROOT format.");
  return 1;
}
