Goal
-----------------

Assess the performance of HGCAL's clustering algorithm (CLUE) with testbeam data and simulation. 

Pipeline description
-----------------

- **1)** *selection stage*: the original NTuples are pruned, in order to keep the relevant information only

- **2)** *analysis stage*:

    - CLUE is run over the pruned NTuples
	
    - most of the quantities of interest are calculated and stored in ```csv``` and ```ROOT``` files

- **3)** *residual analysis and plotting stage*:

    - fits, histogram manipulation and dataframe operations are performed
	
    - quantities of interest are plotted using [BokehPlot](https://bitbucket.org/bfontana/bokehplot), a custom bokeh wrapper (under development)


Steps 1) and 2) were chained with a Directed Acyclic Graph (DAG) that runs within HTCondor.

Input NTuples
------------------

- **data**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/``` 

- **sim_proton**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v5/electrons/```

- **sim_noproton**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v3/electrons/```
	
Analysis type
------------------
Step #3 was divided into independent micro-analysis:

- **hit-level**: energy distributions, calibrations, responses, resolutions

- **layer-level**: CLUE densities and distances (1D and 2D), fractions of clusterized hits and energies 

- **cluster-level**: number of hits, energies, number of clusters, 'x' and 'y' position of the clusters in the detector

Scripts' description
------------------

- ```CondorJobs```: everything related to submitting jobs to the grid

    - ```bin/write_dag.cc```: creates all the required DAG submission file

    - ```selector.sh```: used by the jobs to run step #1

    - ```analyzer.sh```: used by the jobs to run step #2

    - ```clean.sh```: very simple utility that cleans the output files of the jobs once they are not needed

- ```DataProcessing```: everything related to running CLUE and extract its relevant quantities

    - ```src/CLUEAlgo.cc```, ```interface/CLUEAlgo.h``` and some other files in ```interface/```: the CLUE standalone algorithm

    - ```src/CLUEAnalysis.cc``` and ```interface/CLUEAnalysis.h```: calculation of all the quantities of interested from the results obtained by CLUE

    - ```src/selector.cc``` and ```interface/selector.h```: class that manages step #1

    - ```bin/process_data.cc```: executable that runs step #1

    - ```src/analyzer.cc``` and ```interface/analyzer.h```: class that manages step #2

    - ```bin/analyze_data.cc```: executable that runs step #2

    - ```interface/range.h```: utility that allows looping over containers by index

    - ```python/resp_res.py```: run the hit-level analysis type

    - ```python/layer_dep.py```: run the layer-level analysis type

    - ```python/cluster_dep.py```: run the cluster-level analysis type

Standard workflow
-----------------

The macros were written having a particular user in mind, but straightforward adaptations can make it work for other users as well, since the code is reasonably abstract.

If the user wants to process the ```sim_proton``` dataset, he/she would do the following:

- Produce DAG files

```bash
write_dag --datatype sim_proton
```

Alternatively, if only the analysis step is required, once can do

```bash
write_dag --datatype sim_proton --last_step_only
```

- Run the jobs (the submission files will be stored under ```CondorJobs/submission/selection/``` and ```CondorJobs/submission/analysis/```

```bash
condor_submit_dag clue_sim_noproton.dag
```

- Run the python analysis and plotting macros

```bash
python python/resp_res.py --datatype sim_proton --all    #hit level
python python/layer_dep.py --datatype sim_proton --all   #layer level
python python/cluster_dep.py --datatype sim_proton --all #cluster level
```

Please run the scripts with the ```--help``` option for choosing only specific variables.
    
Plots
-----------------
Plots should be [publicly accessible](https://bfontana.web.cern.ch/bfontana/TestBeamReconstruction/).