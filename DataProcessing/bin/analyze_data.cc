#include "UserCode/DataProcessing/interface/analyzer.h"

void analysis_CLUE(const std::string& in_fname, const std::string& out_fname, const std::string& out_fname2, const std::string& out_fname3, const std::string& in_tname) {
  const float dc = 1.3f /*centimeters*/;
  const float kappa = 9.f;
  const float ecut = 3.f;
  /*////////////////////////
    Run custom analyzer
  *////////////////////////
  std::cout << "check" << std::endl;
  Analyzer ana(in_fname, in_tname, dc, kappa, ecut);
  ana.runCLUE();
  //ana.save_to_file(out_fname);
  ana.save_to_file_layer_dependent(out_fname2);
  //ana.save_to_file_cluster_dependent(out_fname3);

  //sum rechit energy directly without clustering
  /*
  bool sum_with_ecut = true;
  ana.sum_energy(sum_with_ecut);
  std::string first_half = out_fname.substr(0,out_fname.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  std::string second_half = out_fname.substr(out_fname.find('.', 20), 4);
  std::string first_half2 = out_fname2.substr(0,out_fname2.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  std::string second_half2 = out_fname2.substr(out_fname2.find('.', 20), 4);
  ana.save_to_file(first_half + "_noclusters" + second_half);
  */
}

void analysis_histos(std::string& in_fname, std::string& out_fname, std::string& in_tname) {
  Analyzer ana(in_fname, in_tname, 0.f, 0.f, 0.f);
  ana.histogram_checks();
}

//run example: analyze_data_exe /eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_437.root out_TEST.csv
int main(int argc, char **argv) {
  const std::string in_tname = "relevant_branches";
  const std::string in_fname = std::string(argv[1]);
  const std::string out_fname = std::string(argv[2]);
  const std::string out_fname2 = std::string(argv[3]);
  const std::string out_fname3 = std::string(argv[4]);  

  const std::string str2 = out_fname2.substr(0,out_fname2.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  const std::string str3 = out_fname3.substr(0,out_fname3.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  const std::string root_end = out_fname3.substr(out_fname3.find('.', 20), 5);
  const std::string out_fname_layer_dependent = str2 + "_layerdep" + root_end;
  const std::string out_fname_cluster_dependent = str3 + "_clusterdep" + root_end; //must end with '.root'!

  analysis_CLUE(in_fname, out_fname, out_fname_layer_dependent, out_fname_cluster_dependent, in_tname);
  //analysis_histos(in_fname, out_fname, in_tname);
  return 0;
}
