Goal
-----------------

Assess the performance of HGCAL's clustering algorithm (CLUE) with testbeam data and simulation. Working under CMSSW 11_1_0_pre2 release.

Status
-----------------

- Electromagnetic showers' studies completed and first draft of the Detector's Note ready for submission
- Now starting to replicate the same studies using hadronic showers

Pipeline description
-----------------

**1)** *selection stage*: the original NTuples are pruned, in order to keep the relevant information only

**2)** *analysis stage*:

- CLUE is run over the pruned NTuples
	
- most of the quantities of interest are calculated and stored in ```csv``` and ```ROOT``` files

**3)** *residual analysis and plotting stage*:

- fits, histogram manipulation and dataframe operations are performed
	
- quantities of interest are plotted using [BokehPlot](https://bitbucket.org/bfontana/bokehplot), a custom bokeh wrapper (under development, but already capable of doing the most common plotting operations)


Steps 1) and 2) were chained with a Directed Acyclic Graph (DAG) that runs within HTCondor.

Input NTuples
------------------

*Electromagnetic showers*

- **data**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v16/``` (HGCAL only)

- **sim_proton** (with proton contamination): ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v5/electrons/```

- **sim_noproton**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v3/electrons/```

*Hadronic showers* (there are no simulations available without proton contamination)

- **data**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/AHCAL_ntuples/v8/``` (HGCAL+AHCAL)

- **sim_proton**: ```/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/sim_ntuples/CMSSW11_0_withAHCAL_newBeamline/FTFP_BERT_EMN/v44_VtxBeam_v3/CorrectFHLay10```
	
Analysis type
------------------
Step #3 was divided into independent micro-analysis:

- **hit-level**: energy distributions, calibrations, responses, resolutions

- **layer-level**: CLUE densities and distances (1D and 2D), fractions of clusterized hits and energies 

- **cluster-level**: number of hits, energies, number of clusters, 'x' and 'y' position of the clusters in the detector

    - spatial resolution: it is possible to summarize its information across multiple tags (datasets with different conditions)

Scripts' description
------------------

- ```CondorJobs/```: everything related to submitting jobs to the grid

    - ```bin/write_dag.cc```: creates all the required DAG submission files. It takes as input the ```CondorJobs/ntuple_ids.txt``` file which lists all the run numbers available and was generated with a combination of ```ls``` (applied to the folder with the input Ntuples) and ```awk```. This file is used in combination with a ```std::map``` stored in ```CondorJobs/interface/run_en_map.h``` which pairs run numbers with incident beam energy in GeV.

    - ```launcher.sh```: used by the jobs to run step #1 and #2

    - ```clean.sh```: very simple utility that cleans the output files of the jobs once they are not needed

    - ```setup.sh```: writes a file named ```ntuple_ids.txt``` which contains the identifiers of the [data ntuples](#input-ntuples) to be considered for the electromagnetic or hadronic analysis 

- ```DataProcessing/```: everything related to running CLUE and extract its relevant quantities

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

    - ```python/summarize_tags.py```: summaryze cluster spatial resolution related quantities for different tags


Standard workflow
-----------------

The macros were written having a particular user in mind, but extremely simple and straightforward adaptations can make it work for other users as well, since the code is reasonably abstract. In particular, running everything over a new dataset should be easy.

If the user wants to process the ```sim_proton``` dataset with electromagnetic showers, he/she should do the following:

- Produce DAG files

```bash
write_dag --datatype sim_proton --showertype em --tag <anything>
```

where ```--tag``` is used for identifying a particular data production step, and should ideally indicate the conditions the data was produced; the data will be stored in a folder named after the tag. Given that the selection stage is seen as something general, the ```--tag``` only affects the analysis step.

> **_WARNING:_** If the same tag is specified more than once, the files will be written in the same folder. If the ```showertype``` and ```datatype``` are also the same, the files will be rewritten, and the old ones lost.

For hadronic showers, ```--showertype had``` is the option to use.
If only the analysis step is required, one can do

```bash
write_dag --datatype sim_proton --showertype em --tag <anything> --last_step_only
```

- Run the jobs (the submission files will be stored under ```CondorJobs/submission/selection/``` and ```CondorJobs/submission/analysis/```

```bash
condor_submit_dag CondorJobs/clue_sim_proton_em_sometag.dag
```

> **_NOTE:_** When the file has already been run Condor automatically uses its *rescue* files, *i.e.*, tries to run only the jobs that did not suceed in previous attempts. To remove all previous files, including job outputs, use ```bash CondorJobs/clean.sh``` (or manually remove the files).


- Join the output files according to their beam energy

```bash
bash DataProcessing/join_ntuples.sh --datatype sim_proton --showertype em --analysistype layerdep --tag <anything>
bash DataProcessing/join_ntuples.sh --datatype sim_proton --showertype em --analysistype clusterdep --tag <anything>
```

The outputs are currently being stored under ```/eos/user/<first username letter>/<username>/TestBeamReconstruction/job_output/```. Please create the required folders if needed. Under ```/job_output/``` the files are stored in the ```hit_dependent/```, ```layer_dependent/``` and ```cluster_dependent/``` folders.

There is no need to join the data of the **hit-level** analysis type, since they are ```.csv``` files joined by the ```pandas``` package. The two other types are instead in ```ROOT``` format and are read by ```uproot```.

*Possible improvement*: one could potentially change the way ```uproot``` reads the files so that it iterates through them (it is potentially faster). This joining step would then become unnecessary.

- Run the python analysis and plotting macros

```bash
python DataProcessing/python/resp_res.py --datatype sim_proton --showertype em --tag <anything>   #hit level
python DataProcessing/python/layer_dep.py --datatype sim_proton --showertype em --tag <anything> --all   #layer level
python DataProcessing/python/cluster_dep.py --datatype sim_proton --showertype em --tag <anything> --all #cluster level
```

If one needs to rerun the plotting stage for ```cluster_dep.py``` simply due to plot formatting or an additional cut, the option ```--use_saved_data``` can be added, effectively speeding-up the macro by reusing the previously stored dataset.

To summarize the X or Y spatial resolution information of multiple tags, calculated by the ```cluster_dep.py``` macro, an additional plotting macro is available:

```bash
python DataProcessing/python/summarize_tags.py --datatype data --showertype em --var dx
```

where the tags to be used have to be specified manually in the macro.

Please run the scripts with the ```--help``` option for more information, including running the last step only for a subset of final variables (this can be done at **layer-level** and **cluster-level**).
    
Plots
-----------------
Plots should be [publicly accessible](https://bfontana.web.cern.ch/bfontana/TestBeamReconstruction/).

Contacts
----------------
Please use CERN Phonebook's details.