#include <cstdlib>
#include "UserCode/DataProcessing/interface/analyzer.h"
#include "UserCode/CondorJobs/interface/run_en_map.h"

//convenience function which prints all the elements in a vector of strings to std::cout
void print_vector_elements(const std::vector<std::string>& v)
{
  for(auto& elem : v) {
    if(elem == v.back())
      std::cout << elem << std::endl;
    else
      std::cout << elem << " / ";
  }
}

struct DataParameters {
  std::string datatype;
  std::string showertype;
  bool last_step_only = false;
  std::string tag;
  std::string w0;
  std::string dpos;
};

//write all individual submission jobs: selection stage
void write_submission_file(const int& id, const std::string& jobpath, const std::string& base, std::string mode,
			   const unsigned int& energy, const DataParameters& p)
{
  std::string cpp_exec;
  if(mode == "selection")
    cpp_exec = "process_data_exe";
  else if(mode == "analysis")
    cpp_exec = "analyze_data_exe";
  else
    throw std::invalid_argument("Mode " + mode + " does not exist.");
  std::string n = std::to_string(id);
  std::ofstream fw(jobpath, std::ios_base::ate);

  std::string memory, flavour;
  if(mode == "selection") {
    memory = "1.3GB";
    flavour = "\"workday\"";
  }
  else {
    memory = "400MB";
    flavour = "\"longlunch\"";
  }
  fw << "indir = " + base << std::endl;
  fw << "executable = launcher.sh" << std::endl;

  fw << "should_transfer_files = YES" << std::endl;
  fw << "when_to_transfer_output = ON_EXIT" << std::endl;
  fw << "public_input_files = ../../../../TestBeamAnalysis/bin/" + std::string(getenv("SCRAM_ARCH")) + "/" + cpp_exec << std::endl;
  
  fw << "arguments = --ntupleid " + n + " --datatype " + p.datatype + " --showertype " + p.showertype + " --tag " + p.tag;
  fw << " --w0 " + p.w0 + " --dpos " + p.dpos;
  fw << " --energy " + std::to_string(energy);
  fw << " --step " + mode;
  fw << std::endl;

  fw << "universe = vanilla" << std::endl;
  fw << "requirements = (OpSysAndVer =?= \"CentOS7\")" << std::endl;

  std::string outname = mode + "_" + p.datatype + "_" + p.showertype + "_" + p.tag + "." + n;
  fw << "output = out/" + outname + ".out" << std::endl;
  fw << "error =  out/" + outname + ".err" << std::endl;
  fw << "log =    log/" + outname + ".log" << std::endl;

  fw << "getenv = True" << std::endl;
  
  fw << "RequestMemory = " + memory << std::endl;
  fw << "+JobFlavour = " + flavour << std::endl;
  fw << "queue" << std::endl;
}

//write Direct Acyclic Graph (DAG) submission file: jobs definition
void write_dag_jobs(const std::string& filepath, const std::vector<std::string>& jobnames, const std::vector<std::string>& jobpaths, const std::ios_base::openmode& mode)
{
  assert( jobnames.size() == jobpaths.size() );
  
  std::ofstream f_write(filepath, mode);
  
  for(auto i: util::lang::indices(jobnames))
    {
      std::string in = "JOB " + jobnames[i] + "\t" + jobpaths[i];
      f_write << in << std::endl;
    }
  f_write << std::endl;
}

//write Direct Acyclic Graph (DAG) submission file: jobs parent/child relationships
void write_dag_hierarchy(const std::string& filepath, const std::vector<std::string>& jobdads, const std::vector<std::string>& jobchilds)
{
  assert( jobchilds.size() == jobdads.size() );

  std::ofstream f_write(filepath, std::ios_base::app);
  
  for(auto i: util::lang::indices(jobchilds))
    {
      std::string in2 = "PARENT " + jobdads[i] + " CHILD " + jobchilds[i];
      f_write << in2 << std::endl;
    }
  f_write << std::endl;
}

void write_dag_repetitions(const std::string& filepath, const std::vector<std::string>& jobnames, const int& nrepetitions)
{
  std::ofstream f_write(filepath, std::ios_base::app);
  
  for(auto i: util::lang::indices(jobnames))
    {
      std::string in = "RETRY " + jobnames[i] + " " + std::to_string(nrepetitions);
      f_write << in << std::endl;
    }
  f_write << std::endl;
}

void write_data(const std::string& submission_folder, const std::string& base, const DataParameters& p)
{
  std::ifstream infile(base + "ntuple_ids.txt");
  std::vector<int> a;
  int _a;
  while (infile >> _a)
      a.push_back(_a);
  //std::vector<int> avoid = {620,621,622};

  std::vector<int> file_id;
  unsigned keymin=run_en_map.cbegin()->first, keymax=run_en_map.crbegin()->first;
  for(unsigned i=keymin; i<=keymax; ++i) //min and max run numbers for data configuration #22
    {
      if( ( p.showertype == "em" and ( (i>=435 and i<=509) or (i>=594 and i<=676) ) )
	  or
	  ( p.showertype == "had" and ( i==321 or (i>=512 and i<=523) or (i>=525 and i<= 593) or (i>=679 and i<=697) ) ) )
	{
	  if( std::find(a.begin(), a.end(), i) != a.end() /*and std::find(avoid.begin(), avoid.end(), i) == avoid.end()*/ )
	    file_id.push_back(i);
	}
    }

  std::string filepath = base + "clue_data_" + p.showertype + "_" + p.tag;
  std::vector<std::string> steps;
  if(p.last_step_only) {
    steps = {"analysis"};
    filepath += "_" + steps[0] + "_only.dag";
  }
  else {
    steps = {"selection", "analysis"};
    filepath += ".dag";
  }
  const unsigned int nsteps = steps.size();
  
  std::vector< std::vector<std::string> > jobnames(nsteps, std::vector<std::string>());
  std::vector< std::vector<std::string> > jobpaths(nsteps, std::vector<std::string>());
  for(auto thisStep: util::lang::indices(steps))
    {
      for(auto i: util::lang::indices(file_id))
	{
	  const std::string n = std::to_string(file_id[i]);
	  const std::string thisJobname = steps[thisStep] + "_data_" + p.showertype + "_beamen" + std::to_string(run_en_map.at(file_id[i])) + "_" + n + "_" + p.tag;
	  jobnames[thisStep].push_back(thisJobname);
	  jobpaths[thisStep].push_back(base + submission_folder + steps[thisStep] + "/" + thisJobname + ".sub");
	}
      if(thisStep == 0)
	write_dag_jobs(filepath, jobnames[thisStep], jobpaths[thisStep], std::ios_base::ate);
      else
	write_dag_jobs(filepath, jobnames[thisStep], jobpaths[thisStep], std::ios_base::app);
    }

  if(!p.last_step_only)
    {
      write_dag_hierarchy(filepath, jobnames[0], jobnames[1]);
      write_dag_repetitions(filepath, jobnames[1], 2);
    }
  write_dag_repetitions(filepath, jobnames[0], 1);

  //write individual analysis jobs which will be submitted all at once by the DAG written above
  for(auto i: util::lang::indices(file_id))
    {
      const unsigned int thisID = file_id[i];
      const unsigned int thisEnergy = run_en_map.at(thisID);
      write_submission_file(thisID, jobpaths[0][i], base, steps[0], thisEnergy, p);
      if(!p.last_step_only)
	write_submission_file(thisID, jobpaths[1][i], base, steps[1], thisEnergy, p); 
    }
}

void write_simulation(const std::string& submission_folder, const std::string& base, const DataParameters& p)
{
  constexpr unsigned int ntuples_per_energy = 5; //ntuples indexes with simulation data range from 0 to 4
  std::vector<unsigned int> energies = {{20, 30, 50, 80, 100, 120, 150, 200, 250, 300}};
  const unsigned int nenergies = energies.size();
  const unsigned int njobs = nenergies*ntuples_per_energy;

  std::string filepath;
  std::vector<std::string> steps;
  if(p.last_step_only) {
    steps = {"analysis"};
    filepath = base + "clue_" + p.datatype + "_" + p.showertype + "_" + p.tag + "_" + steps[0] + "_only.dag";
  }
  else {
    steps = {"selection", "analysis"};
    filepath = base + "clue_" + p.datatype + "_" + p.showertype + ".dag";
  }
  const unsigned int nsteps = steps.size();
  
  std::vector< std::vector<std::string> > jobnames(nsteps, std::vector<std::string>(njobs));
  std::vector< std::vector<std::string> > jobpaths(nsteps, std::vector<std::string>(njobs));
  for(auto thisStep: util::lang::indices(steps)) {
      for(auto i: util::lang::indices(energies)) {
	for(unsigned int j=0; j<ntuples_per_energy; ++j)
	  {
	    const std::string en_str = std::to_string(energies[i]);
	    const std::string idx_str = std::to_string(j);
	    const std::string this_jobname = steps[thisStep] + "_" + p.datatype + "_" + p.showertype + "_beamen" + en_str + "_" + idx_str;
	    jobnames[thisStep].at(j + i*ntuples_per_energy) = this_jobname;
	    jobpaths[thisStep].at(j + i*ntuples_per_energy) = base + submission_folder + steps[thisStep] + "/" + this_jobname + ".sub";
	  }
      }
      if(thisStep == 0)
	write_dag_jobs(filepath, jobnames[thisStep], jobpaths[thisStep], std::ios_base::ate);
      else
	write_dag_jobs(filepath, jobnames[thisStep], jobpaths[thisStep], std::ios_base::app);
    }

  if(!p.last_step_only)
    {
      write_dag_hierarchy(filepath, jobnames[0], jobnames[1]);
      write_dag_repetitions(filepath, jobnames[1], 2);
    }
  write_dag_repetitions(filepath, jobnames[0], 1);

  //write individual analysis jobs which will be submitted all at once by the DAG written above
  for(auto i: util::lang::indices(energies)) {
    for(unsigned int j=0; j<ntuples_per_energy; ++j)
      {
	write_submission_file(j, jobpaths[0][j + i*ntuples_per_energy], base, steps[0], energies[i], p);
	if(!p.last_step_only)
	  write_submission_file(j, jobpaths[1][j + i*ntuples_per_energy], base, steps[1], energies[i], p);
      }
  }
}

int main(int argc, char **argv) {
  std::unordered_map< std::string, std::vector<std::string> > valid_args;
  valid_args["--datatype"] = {"data", "sim_noproton", "sim_proton"};
  valid_args["--showertype"] = {"em", "had"};
  std::vector<std::string> free_args = {"--tag", "--w0", "--dpos"}; //any argument allowed
  std::vector<std::string> optional_args = {"--last_step_only"}; //any argument allowed
  
  int nargsmin = (valid_args.size()+free_args.size()) * 2 + 1;
  if(argc < nargsmin or argc > nargsmin+1) {
    std::cout << "You must specify the following:" << std::endl;
    for(auto& elem : valid_args) {
      std::string elem2 = elem.first;
      elem2.erase(0,2);
      std::cout << elem2 + ": ";
      print_vector_elements(elem.second);
    }
    for(std::string& elem : free_args) {
      std::string elem2 = elem;
      elem2.erase(0,2);
      std::cout << elem2 + ": required, any choice allowed" << std::endl;
    }
    std::cout << "last_step_only: optional" << std::endl;
    return 1;
  }
  for(int iarg=0; iarg<argc; ++iarg) {
    if(std::string(argv[iarg]).find("--") != std::string::npos and
       valid_args.find(std::string(argv[iarg])) == valid_args.end() and
       std::find(free_args.begin(), free_args.end(), std::string(argv[iarg])) == free_args.end() and
       std::find(optional_args.begin(), optional_args.end(), std::string(argv[iarg])) == optional_args.end())
      {
	std::cout << "The arguments currently supported are:" << std::endl;
	for(auto& elem : valid_args)
	  std::cout << elem.first << std::endl;
	for(auto& elem : free_args)
	  std::cout << elem << std::endl;
	return 1;
      }
  }

  DataParameters pars;
  std::unordered_map<std::string, std::string> chosen_args;
  for(int iarg=0; iarg<argc; ++iarg)
    {
      std::string argvstr = std::string(argv[iarg]);
      if( valid_args.find(argvstr) != valid_args.end() ) {
	if( std::find(valid_args[argvstr].begin(), valid_args[argvstr].end(), std::string(argv[iarg+1])) == valid_args[argvstr].end() )
	  {
	    std::cout << "The data type has to be one of the following: ";
	    print_vector_elements(valid_args[argvstr]);
	    return 1;
	  }
	else
	  chosen_args[argvstr] = std::string(argv[iarg+1]);
      }
      else if( std::find(free_args.begin(), free_args.end(), argvstr) != free_args.end() )
	chosen_args[argvstr] = std::string(argv[iarg+1]);
      else if(std::string(argv[iarg]) == "--last_step_only")
	pars.last_step_only = true;
    }
  
  pars.datatype = chosen_args["--datatype"];
  pars.showertype = chosen_args["--showertype"];
  pars.tag = chosen_args["--tag"];
  pars.w0 = chosen_args["--w0"];
  pars.dpos = chosen_args["--dpos"];
  
  //define common variables
  std::string cmssw_base = std::getenv("CMSSW_BASE");
  std::string condorjobs = "/src/UserCode/CondorJobs/";
  std::string submission_folder = "submission/";
  std::string condorjobs_base = cmssw_base + condorjobs;

  //create directories if required
  system( (std::string("mkdir -p ") + condorjobs_base).c_str() );
  system( (std::string("mkdir -p ") + condorjobs_base + std::string("selection/")).c_str() );
  system( (std::string("mkdir -p ") + condorjobs_base + std::string("analysis/")).c_str() );

  //write DAG files
  if(pars.datatype == "data")
    write_data(submission_folder, condorjobs_base, pars);
  else
    write_simulation(submission_folder, condorjobs_base, pars);
}
