#include "TestBeamReconstruction/DataProcessing/interface/selector.h"

int main(int argc, char **argv) {
  std::string input_file = std::string(argv[1]);
  std::string output_file = std::string(argv[2]);
  std::string datatype = std::string(argv[3]);
  std::string showertype = std::string(argv[4]);
  int beam_energy = std::stoi( std::string(argv[5]) );

  Selector selector(input_file, output_file, datatype, showertype, beam_energy);
  selector.select_relevant_branches();
  
  return 0;
}
