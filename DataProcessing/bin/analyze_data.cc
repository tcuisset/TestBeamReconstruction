#include "UserCode/DataProcessing/interface/analyzer.h"

int main(int argc, char **argv) {
  std::string in_fname = std::string(argv[1]);
  std::string out_fname = std::string(argv[2]);
  std::string in_tname = "relevant_branches";

  /*////////////////////////
    Run custom analyzer
  *////////////////////////
  Analyzer ana(in_fname, out_fname, in_tname);
  ana.runCLUE(20, 20, 20, 20); //ana.sum_energy();
  ana.save_to_file(out_fname);
  return 0;
}
