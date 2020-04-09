#include "UserCode/DataProcessing/interface/selector.h"

void selection_phase(const std::string& in_fname, const std::string& out_fname) {
  Selector selector(in_fname, out_fname);
  selector.select_relevant_branches();
}

int main(int argc, char **argv) {
  std::string input_file = std::string(argv[1]);
  std::string output_file = std::string(argv[2]);
  selection_phase(input_file, output_file);
  return 0;
}
