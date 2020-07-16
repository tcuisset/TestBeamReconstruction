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


Steps 1) and 2) were chained with a Directed Acyclic Graph (DAG) that runs n HTCondor.

Input NTuples
------------------

- **data**: 

- **sim_proton**:

- **sim_noproton**:


Analysis type
------------------
Step #3 was divided into independent micro-analysis:

- hit level

- layer level

- cluster level

Scripts list
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
    
Plots
-----------------
Plots should be [publicly accessible](https://bfontana.web.cern.ch/bfontana/TestBeamReconstruction/).