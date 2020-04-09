#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "UserCode/DataProcessing/interface/range.h"
#include "UserCode/DataProcessing/interface/CLUEAlgo.h"
#include "UserCode/DataProcessing/interface/CLUEAnalysis.h"

#ifndef NDEBUG
#   define M_Assert(Expr, Msg) \
    __M_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#   define M_Assert(Expr, Msg) ;
#endif

void __M_Assert(const char* expr_str, bool expr, const char* file, int line, const char* msg)
{
    if (!expr)
    {
        std::cerr << "Assert failed:\t" << msg << "\n"
            << "Expected:\t" << expr_str << "\n"
            << "Source:\t\t" << file << ", line " << line << "\n";
        abort();
    }
}

class Analyzer {
 public:
  Analyzer(const std::vector< std::string >&, const std::string& out_file_path, std::string in_tree_name);
  Analyzer(const std::string&, const std::string& out_file_path, std::string in_tree_name);
  ~Analyzer();
  void runCLUE(float dc, float rhoc_300, float rhoc_200);
  void sum_energy();
  void save_to_file(const std::string&);

 private:
  //methods 
  std::pair<unsigned int, float> _readTree( const std::string&, std::vector< std::vector<float> >& x, std::vector< std::vector<float> >&, std::vector< std::vector<unsigned int> >&, std::vector< std::vector<float> >&, std::vector< std::vector<unsigned int> >&);
  int sanity_checks(const std::string&);
  
  //data
  size_t nfiles_;
  static const int ncpus_ = 4;
  std::vector< std::pair<std::string, std::string> > names_; //file and tree names
  std::vector< std::vector< std::tuple<float, float> > > en_total_; //total energy per event (vector of RecHits) per file (run)
};
