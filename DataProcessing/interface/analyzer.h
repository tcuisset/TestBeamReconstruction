#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <numeric>
#include "TFile.h"
#include "TCanvas.h"
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
  Analyzer(const std::vector< std::string >&, const std::string&, const float&, const float&, const float&);
  Analyzer(const std::string&, const std::string&, const float&, const float&, const float&);
  ~Analyzer();
  void runCLUE();
  void histogram_checks();
  void sum_energy(const bool&);
  void save_to_file(const std::string&);
  void save_to_file_layer_dependent(const std::string&);
  void save_to_file_cluster_dependent(const std::string&);
  
 private:
  //methods 
  std::pair<unsigned int, float> _readTree( const std::string&, std::vector< std::vector<float> >& x, std::vector< std::vector<float> >&, std::vector< std::vector<unsigned int> >&, std::vector< std::vector<float> >&, std::vector< std::vector<unsigned int> >&);
  int sanity_checks(const std::string&);
  bool ecut_selection(const unsigned int&, const float&);
  
  //data
  size_t nfiles_;
  static const int ncpus_ = 4;
  static const unsigned int nlayers_ = 28;
  float dc_, kappa_, ecut_;
  //weights and thickness corrections taken from the third column of Table 3 of CMS DN-19-019
  std::array<float, nlayers_> energy_weights_ = {{11.289,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,9.851,11.360,11.360,11.360,11.360,10.995,10.995,11.153,7.470}};
  std::array<float, 2> thickness_correction_ = {{0.0850, 0.0567}};
  std::vector< std::pair<std::string, std::string> > names_; //file and tree names
  std::vector<float> beam_energies_;
  std::vector< std::vector< std::tuple<float, float> > > en_total_; //total energy per event (vector of RecHits) per file (run) and corresponding beam energy
  std::vector< std::vector< std::array< std::tuple<float, float>, nlayers_> > > fracs_; //fraction of clusterized nhits and clusterized energy per event
  using cluster_dependent_type = std::array< std::tuple< std::vector<unsigned int>, std::vector<float> >, nlayers_>;
  std::vector< std::vector< cluster_dependent_type > > clusterdep_;
};
