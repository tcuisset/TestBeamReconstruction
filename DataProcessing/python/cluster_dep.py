import sys
import os
import errno
import pickle
import subprocess
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline

def get_mean_and_sigma(x, y=None):
    """calculate mean and std for 1D and 2D distributions"""
    if y is None:
        mean = np.mean(x)
        sigma = np.sqrt( np.mean(x - mean)**2 )
    else:
        mean_squared = np.sum(y*x**2)
        mean = np.sum(y*x)
        sumy = np.sum(y)
        mean_squared /= sumy
        mean /= sumy
        sigma = np.sqrt( mean_squared - mean**2 )
    return mean, sigma
    
def get_sigma_band(x, y=None):
    mean, sigma = get_mean_and_sigma(x, y)
    sigma /= np.sqrt(len(x))
    left, right = mean - sigma/2, mean + sigma/2
    return left, right

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

def graphs2d(dfs, axis_kwargs, columns_field, iframe, variable, weight_by_energy=False):
    """Plots cluster-related quantities per layer"""
    if variable not in ('hits', 'energy', 'number'):
        raise ValueError('graphs2d: Variable {} is not supported.'.format(variable))
    for i,df in enumerate(dfs):
        data_per_layer = []

        #determine adequate binning
        ndata = 0
        datamax, datamin = -np.inf, np.inf
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)

            layer_col = get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
            arr = df.loc[:, layer_col]
            if arr.size == 0:
                raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))

            if variable == 'number':
                data_per_layer.append( arr.apply( lambda x: len(x[columns_field+'_layer'+str(ilayer)]), axis=1, result_type='reduce') )
                data_per_layer[-1] = data_per_layer[-1][ data_per_layer[-1] != 0 ] #filter out all events without clusters
            elif variable == 'hits' or variable == 'energy':
                data_per_layer.append( flatten_dataframe(arr) )
            del arr

            ndata += len(data_per_layer[-1])
            m1, m2 = data_per_layer[-1].min(), data_per_layer[-1].max()
            if m1 < datamin:
                datamin = m1
            if m2 > datamax:
                datamax = m2
        nbins = 20 if datamax>20 else int(datamax-1)
        bins = np.linspace(datamin, datamax, nbins+1) #bins = np.logspace(np.log10(m1), np.log10(m2), nbins+1)
        height_hits = height_for_plot_bins(bins, scale='linear')

        #plot data as 2d graphs
        all_counts_hits, all_layers_hits, all_centers_hits = ([] for _ in range(3))
        means, sigmas_l, sigmas_r = ([] for _ in range(3))
        for ilayer in range(1,nlayers+1):
            print('layer ', ilayer)
            
            if weight_by_energy:
                layer_col_weights = get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                data_weights = flatten_dataframe( df.loc[:, layer_col_weights] )
                if data_weights.size == 0:
                    raise ValueError('The weights dataframe {} in layer {} is empty!'.format(i, ilayer))
                counts, edges = np.histogram(data_per_layer[ilayer-1], bins=bins, weights=data_weights)
            else:
                counts, edges = np.histogram(data_per_layer[ilayer-1], bins=bins)
            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts_hits.extend(counts.tolist())
            all_layers_hits.extend(same_layer_array.tolist())
            all_centers_hits.extend(centers.tolist())
            
            #use histogram and not original data to calculate mean and std/sqrt(n)
            #original data cannot be used for the case of weighted histograms
            mean, _ = get_mean_and_sigma(centers, counts)
            sigma_left, sigma_right = get_sigma_band(centers, counts)
            print("sigmas: ", sigma_left, mean, sigma_right)
            means.append(mean)
            sigmas_l.append(sigma_left)
            sigmas_r.append(sigma_right)

        print('graph')

        layers_x = np.arange(1,29)
        interp_thickness = 0.02
        if variable == 'hits':
            interp_degree = 1#2
            interp_smooth = len(layers_x)#len(layers_x) if i<2 else 10*len(layers_x)
        elif variable == 'energy':
            interp_degree = 1
            interp_smooth = len(layers_x)#10*len(layers_x)
        elif variable == 'number':
            interp_degree = 1
            interp_smooth = len(layers_x)/50
        layers_x_fine = np.arange(1,28,interp_thickness)
        means_sp_func = UnivariateSpline(x=layers_x, y=np.array(means), k=interp_degree, s=interp_smooth)
        means_sp = means_sp_func(layers_x_fine)
        sigmasl_sp_func = UnivariateSpline(x=layers_x, y=np.array(sigmas_l), k=interp_degree, s=interp_smooth)
        sigmasl_sp = sigmasl_sp_func(layers_x_fine)
        sigmasr_sp_func = UnivariateSpline(x=layers_x, y=np.array(sigmas_r), k=interp_degree, s=interp_smooth)
        sigmasr_sp = sigmasr_sp_func(layers_x_fine)

        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[i])}
        fig_kwargs.update(axis_kwargs)
        bokehplot.graph(data=[np.array(all_layers_hits), np.array(all_centers_hits), np.array(all_counts_hits)],
                        width=np.ones((len(all_layers_hits))), height=height_hits,
                        idx=i, iframe=iframe, style='rect%Cividis', fig_kwargs=fig_kwargs, alpha=0.6)
        bokehplot.graph(data=[layers_x, means],
                        idx=i, iframe=iframe, color='brown', style='circle', size=2, legend_label='mean')

        #pm = (u'\u00B1').encode('utf-8')
        bokehplot.graph(data=[layers_x_fine, means_sp],
                        width=interp_thickness*np.ones(len(layers_x_fine)), height=abs(sigmasr_sp-sigmasl_sp),
                        idx=i, iframe=iframe+3, color='grey', style='rect', alpha=0.6, legend_label=u'std mean error (\u03c3/\u221an)', 
                        fig_kwargs=fig_kwargs)
        bokehplot.graph(data=[layers_x_fine, means_sp],
                        idx=i, iframe=iframe+3, color='red', style='circle', size=1.5, legend_label='mean')

        del all_counts_hits
        del all_layers_hits
        del all_centers_hits
        del means
        del sigmas_l
        del sigmas_r

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
            self.cache = pickle.load(obj, encoding='bytes')
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

    #load cache
    cacheobj = CacheManager(cache_file_name_hits)
    up_cache = {}#cacheobj.load()

    print('loading data...')
    df = tree.arrays(['*'], outputtype=pd.DataFrame, cache=up_cache)
    print('dumping data...')
    cacheobj.dump()
    print('done!')

    ###############################################
    ######Cluster dependent number of hits#########
    ###############################################
    df_split = []
    for en in beam_energies:
        df_split.append(df[ df[beamen_str] == en])

    axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': '#hits / cluster'}
    graphs2d(df_split, axis_kwargs_hits, columns_field='Nhits', iframe=0, variable='hits', weight_by_energy=True)
    print('after hits_and_energies_graphs2d 1')
    del df_split
    print('save')
    bokehplot.save_frame(iframe=0, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=0, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=0, path='../../DN/figs/', mode='png')
    bokehplot.save_frame(iframe=3, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=3, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=3, path='../../DN/figs/', mode='png')

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
    graphs2d(df_en_split, axis_kwargs_en, columns_field='Energy', iframe=1, variable='energy', weight_by_energy=False)
    print('after hits_and_energies_graphs2d 2')
    print("save")
    bokehplot.save_frame(iframe=1, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=1, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=1, path='../../DN/figs/', mode='png')
    bokehplot.save_frame(iframe=4, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=4, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=4, path='../../DN/figs/', mode='png')

    ###############################################
    ######Number of clusters ######################
    ###############################################
    axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusters'}
    graphs2d(df_en_split, axis_kwargs_en, columns_field='Energy', variable='number', iframe=2, weight_by_energy=False)
    print('after numbers')
    bokehplot.save_frame(iframe=2, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=2, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=2, path='../../DN/figs/', mode='png')
    bokehplot.save_frame(iframe=5, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=5, path=cluster_dep_folder, mode='png')
    bokehplot.save_figs(iframe=5, path='../../DN/figs/', mode='png')
    
if __name__ == '__main__':
    #define local data paths and variables
    eos_base = '/eos/user/'
    cms_user = subprocess.check_output("echo $USER", shell=True, encoding='utf-8').split('\n')[0]
    data_directory = 'TestBeamReconstruction'
    data_path = os.path.join(eos_base, cms_user[0], cms_user, data_directory, "job_output/cluster_dependent/hadd_clusterdep.root")
    beamen_str = 'BeamEnergy'

    #define cache names and paths
    cache_name_hits = 'uproot_cache_hits.pickle'
    cache_file_name_hits = os.path.join(eos_base, cms_user[0], cms_user, data_directory, cache_name_hits)
    cache_name_en = 'uproot_cache_en.pickle'
    cache_file_name_en = os.path.join(eos_base, cms_user[0], cms_user, data_directory, cache_name_en)

    #define analysis constants
    nlayers = 28
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(true_beam_energies_GeV)
    assert(len(beam_energies)==size)

    print("Input data read from {}".format(data_path))

    #create output files with plots
    create_dir( os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory)
    output_html_files = ( os.path.join(output_html_dir, 'plot_clusters_hits.html'),
                          os.path.join(output_html_dir, 'plot_clusters_energy.html'),
                          os.path.join(output_html_dir, 'plot_clusters_number.html'),
                          os.path.join(output_html_dir, 'profile_clusters_hits.html'),
                          os.path.join(output_html_dir, 'profile_clusters_energy.html'),
                          os.path.join(output_html_dir, 'profile_clusters_number.html') )
    nframes = 6
    bokehplot = bkp.BokehPlot(filenames=output_html_files, nfigs=nframes*(size,), nframes=nframes)
    plot_width, plot_height = 600, 400
    cluster_dep_folder = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'layer_dep')
    create_dir( cluster_dep_folder )
    main()
