#include <iostream>
#include <fstream>
#include <string>
#include <optional>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "TestBeamReconstruction/DataProcessing/interface/range.h"
#include "TestBeamReconstruction/DataProcessing/interface/CLUEAnalysis.h"

using mapT = std::map< std::pair<unsigned,unsigned>, float >;

class Selector {
 public:
  Selector(const std::string&, const std::string&, const std::string&, const std::string&, const int&, std::optional<std::string> in_tree_name = std::nullopt, std::optional<std::string> in_tree_name_friend = std::nullopt, std::optional<std::string> out_tree_name = std::nullopt);
  ~Selector();
  void select_relevant_branches();
  void print_relevant_branches(const int&, std::optional<std::string> filename = std::nullopt);

 private:
  int sanity_checks(const std::string&);
  static bool common_selection(const unsigned& layer, const float& energy, const unsigned& chip, const unsigned& channel, const unsigned& module, const float& amplitude, const bool& noise_flag, const mapT& map, const bool& showertype);
  template<typename T> static std::vector<T> clean_ce(const std::vector<T>&, const std::vector<float>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<unsigned>&, const std::vector<float>&, const std::vector<bool>&, const mapT&, const bool&);
  template<typename T> static std::vector<T> clean_ahc(const std::vector<T>&, const std::vector<float>&, const std::vector<int>&, const bool&);
  static std::vector<float> weight_energy_ce(const std::vector<float>&, const std::vector<unsigned>&, const bool&);
  static std::vector<float> weight_energy_ahc(const std::vector<float>&, const bool&);
  void load_noise_values();
  void load_shift_values();
  static bool reject_noise(const mapT& map, const unsigned& mod, const unsigned& chip, const unsigned& l, const float& amp, const bool& st);
  static float ahc_energy_sum(const std::vector<float>&);
  static bool remove_missing_dwc(const std::vector<float>&);
  
  SHOWERTYPE showertype;
  DATATYPE datatype;
  int beam_energy;
  mapT noise_map_;
  std::vector< std::pair<float,float> > shifts_map_;
  
  std::string new_detid_    = "ce_clean_detid";
  std::string new_x_        = "ce_clean_x";
  std::string new_y_        = "ce_clean_y";
  std::string new_z_        = "ce_clean_z";
  std::string new_layer_    = "ce_clean_layer";
  std::string new_en_       = "ce_clean_energy";
  std::string new_en_MeV_   = "ce_clean_energy_MeV";
  std::string new_impX_     = "impactX_shifted";
  std::string new_impY_     = "impactY_shifted";
  //had showers only
  std::string new_ahc_en_     = "ahc_clean_energy";
  std::string new_ahc_en_MeV_ = "ahc_clean_energy_MeV";
  //columns to save
  ROOT::Detail::RDF::ColumnNames_t savedcols_;
  ROOT::Detail::RDF::ColumnNames_t impactcols_, impactXcols_, impactYcols_;

  struct indata {
    std::string file_path = "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/ntuple_1000.root";
    std::string tree_name = "rechitntupler/hits";
    std::string tree_name_friend = "trackimpactntupler/impactPoints";
  } indata_;

  struct outdata {
    std::string file_path = "/eos/user/b/bfontana/TestBeamReconstruction/default_output.txt";
    std::string tree_name = "relevant_branches";
  } outdata_;
};
