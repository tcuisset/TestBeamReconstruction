#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "UserCode/DataProcessing/interface/range.h"
#include "UserCode/DataProcessing/interface/CLUEAnalysis.h"

class Selector {
 public:
  Selector(const std::string&, std::optional<std::string> out_file_path = std::nullopt, std::optional<std::string> in_tree_name = std::nullopt, std::optional<std::string> out_tree_name = std::nullopt);
  ~Selector();
  void select_relevant_branches();
  void print_relevant_branches(const int&, std::optional<std::string> filename = std::nullopt);

 private:
  int sanity_checks(const std::string&);
  
  static const int ncpus_ = 4;

  std::string newcol_clean_detid_     = "rechit_clean_detid";
  std::string newcol_clean_x_         = "rechit_clean_x";
  std::string newcol_clean_y_         = "rechit_clean_y";
  std::string newcol_clean_z_         = "rechit_clean_z";
  std::string newcol_clean_layer_     = "rechit_clean_layer";
  std::string newcol_clean_energy_    = "rechit_clean_energy_MeV";
  //const ROOT::Detail::RDF::ColumnNames_t cols_ = {"event", "run", "NRechits", "rechit_detid", "rechit_x", "rechit_y", "rechit_z", "rechit_layer", "rechit_iu", "rechit_iv", "rechit_iU", "rechit_iV", "rechit_type", "rechit_energy", newcol1_, "beamEnergy", "pdgID"};
  const ROOT::Detail::RDF::ColumnNames_t cols_ = {"event", "run", "NRechits", newcol_clean_detid_, newcol_clean_x_, newcol_clean_y_, newcol_clean_z_, newcol_clean_layer_, newcol_clean_energy_, "beamEnergy"};

  struct indata {
    std::string file_path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_1000.root";
    std::string tree_name = "rechitntupler/hits";
  }indata_;

  struct outdata {
    std::string file_path = "/eos/user/b/bfontana/TestBeamReconstruction/default_output.txt";
    std::string tree_name = "relevant_branches";
  }outdata_;
};
