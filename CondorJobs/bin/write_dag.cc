#include <cstdlib>
#include "UserCode/DataProcessing/interface/analyzer.h"

int main(int argc, char **argv) {
  /*////////////////////////
    Select only the rechits which originated from an electron or positron
    Thorben recommended only using files with configuration 22
  *////////////////////////
  std::string cmssw_base = std::getenv("CMSSW_BASE");
  std::string condorjobs = "/src/UserCode/CondorJobs/";
  std::string submission_folder = "submission/";
  std::string condorjobs_base = cmssw_base + condorjobs;
  std::ifstream infile(condorjobs_base + "ntuple_ids.txt");
  std::string end_str = ".sub";
  std::vector<int> a;
  int _a;
  while (infile >> _a)
      a.push_back(_a);
  //std::vector<int> avoid = {453,508,601,607,639};   //some jobs failed, therefore some input files are lacking

  std::vector<int> file_id;
  for(unsigned int i=435; i<=509; ++i)
    {
      if( std::find(a.begin(), a.end(), i) != a.end() /*and std::find(avoid.begin(), avoid.end(), i) == avoid.end()*/ )
	file_id.push_back(i);
    }
  for(unsigned int i=594; i<=676; ++i)
    {
      if( std::find(a.begin(), a.end(), i) != a.end() /*and std::find(avoid.begin(), avoid.end(), i) == avoid.end()*/ )
	file_id.push_back(i);
    }

  //write Direct Acyclic Graph (DAG) submission file
  std::string filename = "clue";
  std::ofstream f_write(condorjobs_base + filename + ".dag");
  for(auto i: util::lang::indices(file_id))
    {
      std::string n = std::to_string(file_id[i]);
      std::string in = "JOB selection" + n + "\t" + condorjobs_base + submission_folder + "selection/selection" + n + end_str;
      f_write << in << std::endl;
      in = "JOB analysis" + n + "\t\t" + condorjobs_base + submission_folder + "analysis/analysis" + n + end_str;
      f_write << in << std::endl;
    }
  f_write << std::endl;
  for(auto i: util::lang::indices(file_id))
    {
      std::string n = std::to_string(file_id[i]);
      std::string in2 = "PARENT selection" + n + " CHILD analysis" + n;
      f_write << in2 << std::endl;
    }
  f_write << "DOT " + filename + ".dot" << std::endl; //for visualization: dot -Tps clue.dot -o clue.ps

  //write test Direct Acyclic Graph (DAG) submission file: one ntuple only
  std::string filename_test = "clue_TEST";
  std::string n = std::to_string(437); //test number
  std::ofstream f_write_test(condorjobs_base + filename_test + ".dag");
  std::string in = "JOB selection" + n + "\t" + condorjobs_base + submission_folder + "selection/selection" + n + end_str;
  f_write_test << in << std::endl;
  in = "JOB analysis" + n + "\t\t" + condorjobs_base + submission_folder + "analysis/analysis" + n + end_str;
  f_write_test << in << std::endl;
  f_write_test << std::endl;
  std::string in2 = "PARENT selection" + n + " CHILD analysis" + n;
  f_write_test << in2 << std::endl;
  f_write_test << "DOT " + filename_test + ".dot" << std::endl; //for visualization: dot -Tps name.dot -o name.ps

  //write all individual submission jobs: selection stage
  for(auto i: util::lang::indices(file_id))
    {
      std::string n = std::to_string(file_id[i]);
      std::string subname = condorjobs_base + submission_folder + "selection/selection" + n + ".sub";
      std::ofstream f_write_sel(subname);
      f_write_sel << "script_name=selector" << std::endl;
      f_write_sel << std::endl;
      f_write_sel << "executable = " + condorjobs_base + "selector.sh" << std::endl;
      f_write_sel << "arguments = " + n << std::endl;
      f_write_sel << "universe = vanilla" << std::endl;
      f_write_sel << "requirements = (OpSysAndVer =?= \"CentOS7\")" << std::endl;
      f_write_sel << "output = " + condorjobs_base + "out/$(script_name)." + n + ".out" << std::endl;
      f_write_sel << "error = " + condorjobs_base + "out/$(script_name)." + n + ".err" << std::endl;
      f_write_sel << "log = " + condorjobs_base + "log/$(script_name)." + n + ".log" << std::endl;
      f_write_sel << "RequestMemory = 1GB" << std::endl;
      f_write_sel << "+JobFlavour = \"workday\"" << std::endl;
      f_write_sel << "queue" << std::endl;
    }

  //write all individual submission jobs: analysis stage
  for(auto i: util::lang::indices(file_id))
    {
      std::string n = std::to_string(file_id[i]);
      std::string subname = condorjobs_base + submission_folder + "analysis/analysis" + n + ".sub";
      std::ofstream f_write_ana(subname);
      f_write_ana << "script_name=analyzer" << std::endl;
      f_write_ana << std::endl;
      f_write_ana << "executable = " + condorjobs_base + "analyzer.sh" << std::endl;
      f_write_ana << "arguments = " + n << std::endl;
      f_write_ana << "universe = vanilla" << std::endl;
      f_write_ana << "requirements = (OpSysAndVer =?= \"CentOS7\")" << std::endl;
      f_write_ana << "output = " + condorjobs_base + "out/$(script_name)." + n + ".out" << std::endl;
      f_write_ana << "error = " + condorjobs_base + "out/$(script_name)." + n + ".err" << std::endl;
      f_write_ana << "log = " + condorjobs_base + "log/$(script_name)." + n + ".log" << std::endl;
      f_write_ana << "RequestMemory = 80MB" << std::endl;
      f_write_ana << "+JobFlavour = \"workday\"" << std::endl;
      f_write_ana << "queue" << std::endl;
    }

}
