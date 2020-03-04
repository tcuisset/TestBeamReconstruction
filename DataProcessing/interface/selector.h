#include <iostream>
#include <string>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

class Selector {
 public:
  Selector(const std::string&, const std::string&, std::optional<std::string> in_tree_name = std::nullopt, std::optional<std::string> out_tree_name = std::nullopt);
  ~Selector();
  void select_relevant_branches();

 private:
  int sanity_checks(const std::string&);
  struct indata {
    std::string file_path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_1000.root";
    std::string tree_name = "rechitntupler/hits";
  }indata_;
  struct outdata {
    std::string file_path = "/eos/user/b/bfontana/TestBeamReconstruction/output.txt";
    std::string tree_name = "relevant_branches";
  }outdata_;
};
