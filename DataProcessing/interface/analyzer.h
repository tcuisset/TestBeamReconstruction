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
#include "TestBeamReconstruction/DataProcessing/interface/range.h"
#include "TestBeamReconstruction/DataProcessing/interface/CLUEAlgo.h"
#include "TestBeamReconstruction/DataProcessing/interface/CLUEAnalysis.h"

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
  Analyzer(const std::vector< std::string >&, const std::string&, const float&, const float&, const float&, const SHOWERTYPE&, const float&, const float&);
  Analyzer(const std::string&, const std::string&, const float&, const float&, const float&, const SHOWERTYPE&, const float&, const float&);
  ~Analyzer();
  void runCLUE();
  void sum_energy(const bool&);
  void save_to_file(const std::string&);
  void save_to_file_layer_dependent(const std::string&);
  void save_to_file_cluster_dependent(const std::string&);
  
 private:
  //methods 
  /**
   * Fills the vectors in parameters from the input tree of file infile
  */
  std::pair<unsigned int, float> _readTree( const std::string &infile,
                                                   std::vector<std::vector<float>> &x, std::vector<std::vector<float>> &y,
                                                   std::vector<std::vector<unsigned int>> &layer, std::vector<std::vector<float>> &weight,
                                                   std::vector<std::vector<unsigned int>> &rechits_id,
                                                   std::vector<std::vector<float>> &impactX, std::vector<std::vector<float>> &impactY);
  int sanity_checks(const std::string&);
  bool ecut_selection(const float&, const unsigned int&);
  void resize_vectors();
  
  //data
  size_t nfiles_;
  unsigned lmax=0;
  float dc_, kappa_, ecut_;
  SHOWERTYPE st_;
  float W0_, dpos_;
  //weights and thickness corrections taken from the third column of Table 3 of CMS DN-19-019
  std::vector< std::pair<std::string, std::string> > names_; //file and tree names
  std::vector<float> beam_energies_;
  std::vector< std::vector< std::tuple<float, float> > > en_total_; //total energy per event (vector of RecHits) per file (run) and corresponding beam energy
  std::vector< std::vector< dataformats::layerfracs > > layer_fracs_; //fraction of clusterized nhits and clusterized energy per event
  std::vector< std::vector< dataformats::layerhitvars > > layer_hitvars_; //hit-dependent variables that will be plotted in the layer-level analysis: energy, density, distance and isSeed boolena flag
  std::vector< std::vector< dataformats::clustervars > > clusterdep_;
};
