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

def get_layer_col(df, starts_with, ilayer):
    """Obtains list of columns that match a specific column name ending."""
    regex = '^'+starts_with+'.*_layer'+str(ilayer)+'$'
    layer_col = df.columns.str.contains(regex, regex=True)
    assert( single_true(layer_col) )
    return layer_col

def flatten_dataframe(df):
    flat_df = df.values.flatten()
    flat_list = [item for items in flat_df for item in items]
    return np.array(flat_list)

def number_clusters_graph2d(dfs, axis_kwargs, columns_field, iframe):
    for i,df in enumerate(dfs):
        arr_number = []

        #determine adequate binning across all layers
        hits_max, hits_min = -np.inf, np.inf
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)

            layer_col = get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
            arr = df.loc[:, layer_col]
            if arr.size == 0:
                raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))
            arr_number.append( arr.apply( lambda x: len(x[columns_field+'_layer'+str(ilayer)]), axis=1, result_type='reduce') )
            del arr
            arr_number[-1] = arr_number[-1][ arr_number[-1] != 0 ] #filter out all events without clusters

            m1, m2 = arr_number[-1].min(), arr_number[-1].max()
            if m1 < hits_min:
                hits_min = m1
            if m2 > hits_max:
                hits_max = m2
        nbins = int(30.) if hits_max>30 else int(hits_max-1)
        bins = np.linspace(hits_min, hits_max, nbins+1)
        print('NBINS: ', nbins, hits_min, hits_max)
        height_hits = height_for_plot_bins(bins, scale='linear')

        #plot data as 2d graphs
        all_counts_hits, all_layers_hits, all_centers_hits = ([] for _ in range(3))
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)
            counts, edges = np.histogram(arr_number[ilayer-1], bins=bins)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts_hits.extend(counts.tolist())
            all_layers_hits.extend(same_layer_array.tolist())
            all_centers_hits.extend(centers.tolist())
        del arr_number
        print('graph')
        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[i])}
        fig_kwargs.update(axis_kwargs)
        bokehplot.graph(data=[np.array(all_layers_hits), np.array(all_centers_hits), np.array(all_counts_hits)],
                        width=np.ones((len(all_layers_hits))), height=height_hits,
                        idx=i, iframe=iframe, style='rect%Plasma', fig_kwargs=fig_kwargs)
        del all_counts_hits
        del all_layers_hits
        del all_centers_hits

def hits_and_energies_graphs2d(dfs, axis_kwargs, columns_field, iframe, weight_by_energy=False):
    """Plots number of cluysterized hits and clusterized energy per layer"""
    for i,df in enumerate(dfs):
        #determine adequate binning
        ndata = 0
        hits_max, hits_min = -np.inf, np.inf
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)

            layer_col = get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
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
        nbins = int(30.) if hits_max>30 else int(hits_max-1)
        bins = np.linspace(hits_min, hits_max, nbins+1) #bins = np.logspace(np.log10(m1), np.log10(m2), nbins+1)
        height_hits = height_for_plot_bins(bins, scale='linear')

        #plot data as 2d graphs
        all_counts_hits, all_layers_hits, all_centers_hits = ([] for _ in range(3))
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)

            layer_col = get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
            arr = flatten_dataframe( df.loc[:, layer_col] )

            if weight_by_energy:
                layer_col_weights = get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                arr_weights = flatten_dataframe( df.loc[:, layer_col_weights] )
                if arr_weights.size == 0:
                    raise ValueError('The weights dataframe {} in layer {} is empty!'.format(i, ilayer))
                counts, edges = np.histogram(arr, bins=bins, weights=arr_weights)
            else:
                counts, edges = np.histogram(arr, bins=bins)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts_hits.extend(counts.tolist())
            all_layers_hits.extend(same_layer_array.tolist())
            all_centers_hits.extend(centers.tolist())
        print('graph')
        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[i])}
        fig_kwargs.update(axis_kwargs)
        bokehplot.graph(data=[np.array(all_layers_hits), np.array(all_centers_hits), np.array(all_counts_hits)],
                        width=np.ones((len(all_layers_hits))), height=height_hits,
                        idx=i, iframe=iframe, style='rect%Plasma', fig_kwargs=fig_kwargs)
        del all_counts_hits
        del all_layers_hits
        del all_centers_hits

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
    df = tree.arrays(['*'], outputtype=pd.DataFrame, cache=up_cache)
    print('dumping data...')
    cacheobj.dump()
    print('done!')
    df_split = []
    for en in beam_energies:
        df_split.append(df[ df[beamen_str] == en])

    axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': '#hits / cluster'}
    hits_and_energies_graphs2d(df_split, axis_kwargs_hits, columns_field='Nhits', iframe=0, weight_by_energy=True)
    del df_split
    bokehplot.save_frame(iframe=0, plot_width=plot_width, plot_height=plot_height, show=False)

    ###############################################
    ######Cluster dependent energy#################
    ###############################################
    en_cols   = [x for x in df.columns if 'Energy' in x] #only need a subset of the original data
    df_en = df[en_cols]
    df_en_split = []
    for en in beam_energies:
        df_en_split.append(df_en[ df_en[beamen_str] == en])    
    del df_en #clears RAM

    axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Total energy per cluster [MeV]'}
    hits_and_energies_graphs2d(df_en_split, axis_kwargs_en, columns_field='Energy', iframe=1, weight_by_energy=False)
    bokehplot.save_frame(iframe=1, plot_width=plot_width, plot_height=plot_height, show=False)

    ###############################################
    ######Number of hits per cluster###############
    ###############################################
    axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusters'}
    number_clusters_graph2d(df_en_split, axis_kwargs_en, columns_field='Energy', iframe=2)
    bokehplot.save_frame(iframe=2, plot_width=plot_width, plot_height=plot_height, show=False)

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

    print("Input data read from {}".format(data_path))

    #create output files with plots
    create_dir( os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory)
    output_html_files = ( os.path.join(output_html_dir, 'plot_clusters_hits.html'),
                          os.path.join(output_html_dir, 'plot_clusters_energy.html'),
                          os.path.join(output_html_dir, 'plot_clusters_number.html') )
    nframes = 3
    bokehplot = bkp.BokehPlot(filenames=output_html_files, nfigs=(size, size, size), nframes=nframes)
    plot_width, plot_height = 600, 400
    main()

