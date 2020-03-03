#include <iostream>
#include <string>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

class Selector {
 public:
  Selector();
  ~Selector();
  void select_relevant_branches();

 private:
  struct indata {
    std::string file_path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/";
    std::string tree_name = "rechitntupler/hits";
  }indata_;
  struct outdata {
    std::string file_path = "output.root";
    std::string tree_name = "relevant_branches";
  }outdata_;
};
