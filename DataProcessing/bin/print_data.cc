#include "UserCode/DataProcessing/interface/selector.h"

void print_tree_contents(const std::string& in_fname) {
  Selector selector(in_fname, "", "", "", 0);
  selector.print_relevant_branches(1);
}

int main(int argc, char **argv) {
  print_tree_contents(std::string(argv[1]));
  return 0;
}
