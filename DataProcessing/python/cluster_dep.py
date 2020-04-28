import sys, os
import subprocess
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd

def histograms(dfs1, dfs2, axis_kwargs):
    """Plots number of cluysterized hits and clusterized energy per layer"""
    assert(len(dfs1) == len(dfs2))
    hist1, hist2 = ([] for _ in range(2))
    bins = np.linspace(0,100,11)

    for ilayer in range(1,nlayers+1):
        arr = dfs1[-1]['Nhits_layer'+str(ilayer)]
        counts, edges = np.histogram([item for items in arr for item in items], bins=bins)
        centers = (edges[:-1]+edges[1:])/2
        print(centers)
        quit()

    for i, idf in enumerate(dfs1):
        nhitsfrac = idf['nhitsfrac']
        enfrac    = idf['endfrac']
        #flattens 'idf' to one-dimension
        hist1.append( np.histogram2d(nhitsfrac, density=False, bins=bins, range=(nhitsfrac.min(),nhitsfrac.max())) )
        hist2.append( np.histogram2d(enfrac, density=False, bins=bins, range=(nhitsfrac.min(),nhitsfrac.max())) )

    """
    #number of hits
    bokehplot.histogram(data=hist1, 
                        idx=[x for x in range(size)], style='hex%Viridis',
                        fig_kwargs=axis_kwargs[0])
    #energy
    bokehplot.histogram(data=hist2, 
                        idx=[x for x in range(size,2*size)], style='hex%Viridis',
                        fig_kwargs=axis_kwargs[1])

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
    df = tree.arrays("*", outputtype=pd.DataFrame, entrystop=10, cache=up_cache)
    hits_cols = [x for x in df.columns if 'Nhits'  in x] + [beamen_str]
    en_cols   = [x for x in df.columns if 'Energy' in x]
    df_hits, df_en = df[hits_cols], df[en_cols]
    df_hits_split, df_en_split = ([] for _ in range(2))
    for en in beam_energies:
        df_hits_split.append(df_hits[ df_hits[beamen_str] == en])
        df_en_split.append(df_en[ df_en[beamen_str] == en])    

    axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusterized hits'}
    axis_kwargs_en   = {'x.axis_label': 'Layer', 'y.axis_label': 'Clusterized energy'}
    histograms(df_hits_split, df_en_split, (axis_kwargs_hits, axis_kwargs_en))

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    nlayers = 28
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(true_beam_energies_GeV)
    assert(len(beam_energies)==size)

    bokehplot = bkp.BokehPlot(filenames='plot_clusters.html', nfigs=2*size)
    main()
