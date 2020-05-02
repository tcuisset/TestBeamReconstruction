import sys
import os
import errno
import pickle
import subprocess
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd

def create_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
        
def height_for_plot_bins(bins, scale='linear'):
    nbins = len(bins)-1
    if scale == 'linear':
        factor = abs(bins[-1]-bins[0]) / nbins
        return factor * np.ones(nbins*nlayers)
    elif scale == 'log':
        centers = (bins[1:]+bins[:-1])/2
        distances = [(centers[1]-centers[0])/2]
        distances.extend( (centers[1:]-centers[:-1])/2 )
        return np.array(distances*nlayers)
    else:
        raise ValueError('height_for_plot_bins: Option not supported.')

def single_true(iterable):
    """Checks if one and only one element of iterable is True"""
    i = iter(iterable)
    return any(i) and not any(i)

def get_layer_col(df, ilayer):
    """Obtains list of columns that match a specific column name ending."""
    layer_col = df.columns.str.endswith('_layer'+str(ilayer))
    assert( single_true(layer_col) )
    return layer_col

def flatten_dataframe(df):
    flat_df = df.values.flatten()
    flat_list = [item for items in flat_df for item in items]
    return np.array(flat_list)

def graphs2d(dfs, axis_kwargs):
    """Plots number of cluysterized hits and clusterized energy per layer"""
    for i,df in enumerate(dfs):
        #determine adequate binning
        ndata = 0
        hits_max, hits_min = -np.inf, np.inf
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)
            layer_col = get_layer_col(df, ilayer)
            arr = df.loc[:, layer_col]
            if arr.size == 0:
                raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))

            arr = flatten_dataframe(arr)
            ndata += len(arr)
            m1, m2 = arr.min(), arr.max()
            if m1 < hits_min:
                hits_min = m1
            if m2 > hits_max:
                hits_max = m2
        nbins = int(15.) if hits_max>=15 else int(hits_max)
        bins = np.linspace(m1, m2, nbins+1)
        #bins = np.logspace(np.log10(m1), np.log10(m2), nbins+1)
        height_hits = height_for_plot_bins(bins, scale='linear')

        #plot data as 2d graphs
        all_counts_hits, all_layers_hits, all_centers_hits = ([] for _ in range(3))
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)
            layer_col = get_layer_col(df, ilayer)
            arr = df.loc[:, layer_col]
            arr = flatten_dataframe(arr)
            counts, edges = np.histogram(arr, bins=bins)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts_hits.extend(counts.tolist())
            all_layers_hits.extend(same_layer_array.tolist())
            all_centers_hits.extend(centers.tolist())
        print('graph')
        print(len(all_layers_hits))
        print(len(all_centers_hits))
        print(len(all_counts_hits))
        fig_kwargs = {'plot_width': 800, 'plot_height': 500,
                      't.text': 'True beam energy: {} GeV'.format(true_beam_energies_GeV[i])}
        fig_kwargs.update(axis_kwargs)
        bokehplot.graph(data=[np.array(all_layers_hits), np.array(all_centers_hits), np.array(all_counts_hits)],
                        width=np.ones((len(all_layers_hits))), height=height_hits,
                        idx=i, style='rect%Viridis', fig_kwargs=fig_kwargs)

class CacheManager:
    def __init__(self, name):
        self.name_ = name
        self.cache = {}

    def dump(self):
        obj = open(self.name_, 'wb')
        pickle.dump(self.cache, obj)
        obj.close()

    def load(self):
        try:
            obj = open(self.name_, 'rb')
            self.cache = pickle.load(obj)
            obj.close()
        except IOError:
            pass
        except EOFError:
            print('Perhaps the cache was not properly saved in a previous session?')
            raise
        return self.cache

def main():
    #load ROOT TTree
    file = up.open( data_path )
    tree = file['tree0']

    ###############################################
    ######Cluster dependent number of hits#########
    ###############################################
    #load cache
    cacheobj = CacheManager(cache_file_name_hits)
    up_cache = cacheobj.load()

    print('loading data...')
    df = tree.arrays([beamen_str, 'Nhits_layer*'], outputtype=pd.DataFrame, cache=up_cache)
    print('dumping data...')
    cacheobj.dump()
    print('done!')
    hits_cols = [x for x in df.columns if 'Nhits'  in x] + [beamen_str]
    df_hits = df[hits_cols]
    df_hits_split = []
    for en in beam_energies:
        df_hits_split.append(df_hits[ df_hits[beamen_str] == en])

    axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of hits per cluster'}
    graphs2d(df_hits_split, axis_kwargs_hits)
    print('save_frame')
    bokehplot.save_frame(plot_width=800, plot_height=500, show=False)

    ###############################################
    ######Cluster dependent energy#################
    ###############################################
    bokehplot.add_frame(os.path.join(output_html_dir, 'another_plot.html'), nfigs=size)

    #load cache
    #up_cache.clear()
    cacheobj = CacheManager(cache_file_name_en)
    up_cache = cacheobj.load()

    print('loading data...')
    df = tree.arrays([beamen_str, 'Energy_layer*'], outputtype=pd.DataFrame, cache=up_cache)
    print('dumping data...')
    cacheobj.dump()
    print('done!')
    en_cols   = [x for x in df.columns if 'Energy' in x]
    df_en = df[en_cols]
    df_en_split = []
    for en in beam_energies:
        df_en_split.append(df_en[ df_en[beamen_str] == en])    

    axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Total energy per cluster [MeV]'}
    graphs2d(df_en_split, axis_kwargs_en)
    print('save_frame')
    bokehplot.save_frame(plot_width=800, plot_height=500, show=False)

if __name__ == '__main__':
    #define local data paths and variables
    cmssw_base = '/afs/cern.ch/user/b/bfontana/CMSSW_11_1_0_pre2/'
    usercode_path = 'src/UserCode/DataProcessing/job_output/cluster_dependent/'
    data_path = os.path.join(cmssw_base, usercode_path, "hadd_clusterdep.root")
    beamen_str = 'BeamEnergy'
    cms_user = subprocess.check_output("echo $USER", shell=True).split('\n')[0]
    data_directory = 'TestBeamReconstruction'
    create_dir( os.path.join('/eos/user/', cms_user[0], cms_user, data_directory) )

    #define cache names and paths
    cache_name_hits = 'uproot_cache_hits.pickle'
    cache_file_name_hits = os.path.join('/eos/user/', cms_user[0], cms_user, data_directory, cache_name_hits)
    cache_name_en = 'uproot_cache_en.pickle'
    cache_file_name_en = os.path.join('/eos/user/', cms_user[0], cms_user, data_directory, cache_name_en)

    #define analysis constants
    nlayers = 28
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(true_beam_energies_GeV)
    assert(len(beam_energies)==size)

    #create output files with plots
    create_dir( os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory)
    output_html_file = os.path.join(output_html_dir, 'plot_clusters.html')
    bokehplot = bkp.BokehPlot(filenames=output_html_file, nfigs=size)
    
    #launch plotting code for cluster dependent quantities
    print("Input data read from {}.".format(data_path))
    print("Plots will be saved in {}.".format(output_html_file))
    main()
