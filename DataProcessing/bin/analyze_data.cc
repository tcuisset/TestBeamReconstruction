#include "UserCode/DataProcessing/interface/analyzer.h"

void analysis_CLUE(std::string& in_fname, std::string& out_fname, std::string& out_fname2, std::string& in_tname) {
  float dc = 1.3 /*centimeters*/;
  float kappa = 9;
  std::array<float, 2> snratio = {{7.2 /*300 micron sensors*/, 4.8 /*200 micron sensors*/}}; //values given by Thorben

  /*////////////////////////
    Run custom analyzer
  *////////////////////////
  Analyzer ana(in_fname, in_tname);

  ana.runCLUE(dc, kappa * snratio[0], kappa * snratio[1]);
  ana.save_to_file(out_fname);
  ana.save_to_file_layer_dependent(out_fname2);
  //sum rechit energy directly without clustering
  ana.sum_energy();
  std::string first_half = out_fname.substr(0,out_fname.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  std::string second_half = out_fname.substr(out_fname.find('.', 20), 4);
  ana.save_to_file(first_half+"_noclusters"+second_half);
  std::string first_half2 = out_fname.substr(0,out_fname2.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  std::string second_half2 = out_fname.substr(out_fname2.find('.', 20), 4);
  ana.save_to_file(first_half2+"_noclusters"+second_half2);
}

void analysis_histos(std::string& in_fname, std::string& out_fname, std::string& in_tname) {
  Analyzer ana(in_fname, in_tname);
  ana.histogram_checks();
}

//run example: analyze_data_exe /eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_437.root out_TEST.csv
int main(int argc, char **argv) {
  std::string in_fname = std::string(argv[1]);
  std::string out_fname = std::string(argv[2]);
  std::string out_fname_layer_dependent = "layer_dep_" + std::string(argv[3]);
  std::string in_tname = "relevant_branches";

  analysis_CLUE(in_fname, out_fname, out_fname_layer_dependent, in_tname);
  //analysis_histos(in_fname, out_fname, in_tname);
  return 0;
}
