#include <cstdlib>
#include "UserCode/DataProcessing/interface/analyzer.h"

int main(int argc, char **argv) {
  /*////////////////////////
    Select only the rechits which originated from an electron or positron
    Thorben recommended only using files with configuration 22
  *////////////////////////
  std::string cmssw_base = std::getenv("CMSSW_BASE");
  std::ifstream infile(cmssw_base + "/src/UserCode/CondorJobs/ntuple_ids.txt");
  std::vector<int> a;
  int _a;
  while (infile >> _a)
      a.push_back(_a);
  //some jobs failed, therefore some input files are lacking
  std::vector<int> avoid = {453,508,601,607,639};

  std::vector<int> file_id;
  for(unsigned int i=435; i<=509; ++i)
    {
      if( std::find(a.begin(), a.end(), i) != a.end() and std::find(avoid.begin(), avoid.end(), i) == avoid.end() )
	file_id.push_back(i);
    }
  for(unsigned int i=594; i<=676; ++i)
    {
      if( std::find(a.begin(), a.end(), i) != a.end() and std::find(avoid.begin(), avoid.end(), i) == avoid.end() )
	file_id.push_back(i);
    }
  std::string init_str = "/eos/user/b/bfontana/TestBeamReconstruction/ntuple_selection_";
  std::string end_str = ".root";
  std::vector<std::string> in_fname(file_id.size());
  std::ofstream f_write(cmssw_base + "/src/UserCode/CondorJobs/file_names.txt");
  for(auto i: util::lang::indices(file_id))
    {
      in_fname[i] = init_str + std::to_string(file_id[i]) + end_str;
      std::ifstream f(in_fname[i].c_str());
      M_Assert( f.good(), in_fname[i].c_str());
      f_write << in_fname[i] << std::endl;
    }
}
