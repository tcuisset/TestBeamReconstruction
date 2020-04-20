#include "UserCode/DataProcessing/interface/analyzer.h"

void analysis_CLUE(std::string& in_fname, std::string& out_fname, std::string& in_tname) {
  float dc = 1.3 /*centimeters*/;
  float kappa = 9;
  std::array<float, 2> snratio = {{7.2 /*300 micron sensors*/, 4.8 /*200 micron sensors*/}}; //values given by Thorben

  /*////////////////////////
    Run custom analyzer
  *////////////////////////
  Analyzer ana(in_fname, out_fname, in_tname);

  ana.runCLUE(dc, kappa * snratio[0], kappa * snratio[1]);
  ana.save_to_file(out_fname);
  //sum rechit energy directly without clustering
  ana.sum_energy();
  std::string first_half = out_fname.substr(0,out_fname.find('.', 20)); //the 20 avoids the '.' in 'cern.ch'
  std::string second_half = out_fname.substr(out_fname.find('.', 20), 4);
  ana.save_to_file(first_half+"_noclusters"+second_half);
}

void analysis_histos(std::string& in_fname, std::string& out_fname, std::string& in_tname) {
  Analyzer ana(in_fname, out_fname, in_tname);
  ana.histogram_checks();
}

//run example: analyze_data_exe /eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_437.root out_TEST.csv
int main(int argc, char **argv) {
  std::string in_fname = std::string(argv[1]);
  std::string out_fname = std::string(argv[2]);
  std::string in_tname = "relevant_branches";

  analysis_CLUE(in_fname, out_fname, in_tname);
  //analysis_histos(in_fname, out_fname, in_tname);
  return 0;
}
