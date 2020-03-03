#include "UserCode/DataProcessing/interface/selector.h"

Selector::Selector()
{
}

Selector::~Selector()
{
}

void Selector::select_relevant_branches()
{
  const int ncpus = 4;
  ROOT::EnableImplicitMT( ncpus );
  ROOT::RDataFrame d(this->indata_.tree_name.c_str(), this->indata_.file_path.c_str());
  d.Snapshot(this->outdata_.tree_name.c_str(), this->outdata_.file_path.c_str(), {"rechit_x", "rechit_y", "rechit_z", "rechit_layer", "rechit_iu", "rechit_iv", "rechit_iU", "rechit_iV", "rechit_type", "rechit_energy", "beamEnergy", "pdgID"});
};
