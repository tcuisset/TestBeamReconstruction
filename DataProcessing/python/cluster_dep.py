import sys, os
import subprocess
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd

def graphs2d(dfs1, dfs2, axis_kwargs):
    """Plots number of cluysterized hits and clusterized energy per layer"""
    assert(len(dfs1) == len(dfs2))
    bins_hits = np.linspace(0,50,11)
    bins_en   = np.linspace(0,5000,11)

    for i,df in enumerate(dfs1):
        for ilayer in range(1,nlayers+1):
            arr = df['Nhits_layer'+str(ilayer)]
            counts, edges = np.histogram([item for items in arr for item in items], bins=bins_hits)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            bokehplot.graph(data=[same_layer_array, centers, counts], 
                            idx=i, style='square%Viridis',
                            fig_kwargs=axis_kwargs[0])
            print(i)
            print(centers, type(centers))
            print(counts, type(counts))
    bokehplot.show_frame(plot_width=300, plot_height=300)

    """
    bokehplot.add_frame("another_plot.html", nfigs=size)
    for i,df in enumerate(dfs2):
        for ilayer in range(1,nlayers+1):
            arr = df['Energy_layer'+str(ilayer)]
            counts, edges = np.histogram([item for items in arr for item in items], bins=bins_en)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            bokehplot.graph(data=[same_layer_array, centers, counts], 
                            idx=i, style='circle%Viridis',
                            fig_kwargs=axis_kwargs[1])
            print(centers)
            print(counts)
    bokehplot.show_frame(plot_width=200, plot_height=200)
    """

def main():
    usercode_path = 'src/UserCode/DataProcessing/job_output/cluster_dependent/'
    path = os.path.join(cmssw_base, usercode_path, "hadd_clusterdep.root")
    beamen_str = 'BeamEnergy'
    up_cache = {}

    file = up.open( path )
    #file.allkeys(filterclass=lambda x: issubclass(x, up.tree.TTreeMethods))
    tree = file['tree0']
    df = tree.arrays("*", outputtype=pd.DataFrame, entrystop=200, cache=up_cache)
    hits_cols = [x for x in df.columns if 'Nhits'  in x] + [beamen_str]
    en_cols   = [x for x in df.columns if 'Energy' in x]
    df_hits, df_en = df[hits_cols], df[en_cols]
    df_hits_split, df_en_split = ([] for _ in range(2))
    for en in beam_energies:
        df_hits_split.append(df_hits[ df_hits[beamen_str] == en])
        df_en_split.append(df_en[ df_en[beamen_str] == en])    

    axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusterized hits'}
    axis_kwargs_en   = {'x.axis_label': 'Layer', 'y.axis_label': 'Clusterized energy'}
    graphs2d(df_hits_split, df_en_split, (axis_kwargs_hits, axis_kwargs_en))

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    nlayers = 28
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(true_beam_energies_GeV)
    assert(len(beam_energies)==size)

    bokehplot = bkp.BokehPlot(filenames='plot_clusters.html', nfigs=size)
    main()
