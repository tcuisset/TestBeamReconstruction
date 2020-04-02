#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"

class Selector {
 public:
  Selector(const std::string&, std::optional<std::string> out_file_path = std::nullopt, std::optional<std::string> in_tree_name = std::nullopt, std::optional<std::string> out_tree_name = std::nullopt);
  ~Selector();
  void select_relevant_branches();
  void print_relevant_branches(const int&, std::optional<std::string> filename = std::nullopt);

 private:
  int sanity_checks(const std::string&);
  
  static const int ncpus_ = 4;
  //weights and thickness corrections taken from the third column of Table 3 of CMS DN-19-019
  std::array<float, 28> energy_weights_ = {{11.289,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,11.360,11.360,11.360,11.360,10.995,10.995,11.153,7.470}};
  std::array<float, 2> thickness_correction_ = {{0.0850, 0.0567}};

  std::string newcol1_ = "rechit_energy_MeV";
  std::string newcol2_ = "rechit_weighted_energy_MeV";
  const ROOT::Detail::RDF::ColumnNames_t cols_ = {"event", "run", "NRechits", "rechit_detid", "rechit_x", "rechit_y", "rechit_z", "rechit_layer", "rechit_iu", "rechit_iv", "rechit_iU", "rechit_iV", "rechit_type", "rechit_energy", newcol1_, newcol2_, "beamEnergy", "trueBeamEnergy", "pdgID"};
  struct indata {
    std::string file_path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_1000.root";
    std::string tree_name = "rechitntupler/hits";
  }indata_;
  struct outdata {
    std::string file_path = "/eos/user/b/bfontana/TestBeamReconstruction/default_output.txt";
    std::string tree_name = "relevant_branches";
  }outdata_;
};
