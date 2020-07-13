import sys
import os
import errno
import pickle
import subprocess
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd
import concurrent.futures
import argparse
from argparser import add_args
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
        if sumy == 0:
            mean = 0
            sigma = np.inf
        else:
            mean_squared /= sumy
            mean /= sumy
            sigma = np.sqrt( mean_squared - mean**2 )
    return mean, sigma
    
def get_sigma_band(x, y=None):
    mean, sigma = get_mean_and_sigma(x, y)
    if mean==0 and sigma==np.inf:
        left = mean
        right = mean
    else:
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

def get_layer_col(df, starts_with, ilayer=None):
    """Obtains list of columns that match a specific column name ending."""
    if ilayer is None:
        regex = '^'+starts_with
    else:
        regex = '^'+starts_with+'.*_layer'+str(ilayer)+'$'
    layer_col = df.columns.str.contains(regex, regex=True)
    if ilayer is not None:
        assert( single_true(layer_col) )
    return layer_col

def flatten_dataframe(df):
    flat_df = df.values.flatten()
    flat_list = [item for items in flat_df for item in items]
    return np.array(flat_list)

def graphs2d(dfs, axis_kwargs, columns_field, iframe, variable, weight_by_energy=False):
    """Plots cluster-related quantities per layer"""
    if variable not in ('hits', 'energy', 'number', 'pos'):
        raise ValueError('graphs2d: Variable {} is not supported.'.format(variable))
    for i,df in enumerate(dfs):
        print('Processing the {} GeV dataset...'.format(beam_energies[i]))
        data_per_layer = []
        data_per_layer_cut = []

        datamax, datamin = -np.inf, np.inf
        for ilayer in range(1,nlayers+1):
            #determine adequate binning
            layer_col = get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
            arr = df.loc[:, layer_col]
            if arr.size == 0:
                raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))

            energy_lower_cut = 0#1500
            if variable == 'number':
                data_per_layer.append( arr.apply( lambda x: len(x[columns_field+'_layer'+str(ilayer)]), axis=1, result_type='reduce') )
                data_per_layer[-1] = data_per_layer[-1][ data_per_layer[-1] != 0 ] #filter out all events without clusters
                data_per_layer_cut.append( arr.apply( lambda x: len(x[columns_field+'_layer'+str(ilayer)][ x[columns_field+'_layer'+str(ilayer)] > energy_lower_cut ]), axis=1, result_type='reduce') )
                data_per_layer_cut[-1] = data_per_layer[-1][ data_per_layer[-1] != 0 ] #filter out all events without clusters
            elif variable == 'energy':
                arr = flatten_dataframe(arr)
                data_per_layer.append( arr )
                arr = arr[ arr > energy_lower_cut ]
                data_per_layer_cut.append( arr )
            elif variable == 'hits' or variable == 'pos':
                extra_layer_col = get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                extra_arr = df.loc[:, extra_layer_col]
                arr = flatten_dataframe(arr)
                data_per_layer.append( arr )
                arr = arr[ extra_arr > energy_lower_cut ]
                data_per_layer_cut.append( arr )
            del arr

            if len(data_per_layer[-1]) != 0:
                m1, m2 = data_per_layer[-1].min(), data_per_layer[-1].max()
                if m1 < datamin:
                    datamin = m1
                if m2 > datamax:
                    datamax = m2
        if variable == 'pos':
            nbins = 50 
        else:
            nbins = 20 if datamax>20 else int(datamax-1)
        bins = np.linspace(datamin, datamax, nbins+1)
        height_hits = height_for_plot_bins(bins, scale='linear')
        all_counts, all_layers, all_centers = ([] for _ in range(3))
        means, sigmas_l, sigmas_r = ([] for _ in range(3))

        #plot data as 2d graphs
        for ilayer in range(1,nlayers+1):                
            if weight_by_energy:
                layer_col_weights = get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                data_weights = flatten_dataframe( df.loc[:, layer_col_weights] )
                if data_weights.size == 0:
                    raise ValueError('The weights dataframe {} in layer {} is empty!'.format(i, ilayer))
                counts, edges = np.histogram(data_per_layer[ilayer-1], bins=bins, weights=data_weights)
                counts_cut, _ = np.histogram(data_per_layer_cut[ilayer-1], bins=bins, weights=data_weights)
            else:
                counts, edges = np.histogram(data_per_layer[ilayer-1], bins=bins)
                counts_cut, _ = np.histogram(data_per_layer_cut[ilayer-1], bins=bins)

            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts.extend(counts.tolist())
            all_counts_cut.extend(counts_cut.tolist())
            all_layers.extend(same_layer_array.tolist())
            all_centers.extend(centers.tolist())
            
            #use histogram and not original data to calculate mean and std/sqrt(n)
            #original data cannot be used for the case of weighted histograms
            mean, _ = get_mean_and_sigma(centers, counts)
            sigma_left, sigma_right = get_sigma_band(centers, counts)
            means.append(mean)
            sigmas_l.append(sigma_left)
            sigmas_r.append(sigma_right)

        del data_per_layer
        
        layers_x = np.arange(1,29)
        interp_thickness = 0.02
        if variable == 'hits':
            interp_degree = 1#2
            interp_smooth = len(layers_x)#len(layers_x) if i<2 else 10*len(layers_x)
        elif variable == 'energy':
            interp_degree = 1
            interp_smooth = len(layers_x)#10*len(layers_x)
        elif variable == 'number' or variable == 'pos':
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
        bokehplot.graph(data=[np.array(all_layers), np.array(all_centers), np.array(all_counts)],
                        width=np.ones((len(all_layers))), height=height_hits,
                        idx=i, iframe=iframe, style='rect%Cividis', fig_kwargs=fig_kwargs, alpha=0.6)
        bokehplot.graph(data=[layers_x, means],
                        idx=i, iframe=iframe, color='brown', style='circle', size=2, legend_label='mean')
        #with cut
        bokehplot.graph(data=[np.array(all_layers), np.array(all_centers), np.array(all_counts_cuts)],
                        width=np.ones((len(all_layers))), height=height_hits,
                        idx=i+size, iframe=iframe, style='rect%Cividis', fig_kwargs=fig_kwargs, alpha=0.6)

        bokehplot.graph(data=[layers_x_fine, means_sp],
                        width=interp_thickness*np.ones(len(layers_x_fine)), height=abs(sigmasr_sp-sigmasl_sp),
                        idx=i, iframe=iframe+3, color='grey', style='rect', alpha=0.6, legend_label=u'std mean error (\u03c3/\u221an)', 
                        fig_kwargs=fig_kwargs)  #pm = (u'\u00B1').encode('utf-8')
        bokehplot.graph(data=[layers_x_fine, means_sp],
                        idx=i, iframe=iframe+3, color='red', style='circle', size=1.5, legend_label='mean')

        del all_counts
        del all_counts_cut
        del all_layers
        del all_centers
        del means
        del sigmas_l
        del sigmas_r

class CacheManager:
    def __init__(self, name):
        self.name_ = name
        self.cache = up.ArrayCache("1 GB")

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

def save_plots(frames_to_save):
    mode = 'png'
    second_folder = '../../DN/figs/'
    #bokehplot.save_frame(iframe=0, plot_width=plot_width, plot_height=plot_height, show=False)
    #bokehplot.save_frame(iframe=3, plot_width=plot_width, plot_height=plot_height, show=False)
    bokehplot.save_figs(iframe=frames_to_save[0], path=cluster_dep_folder, mode=mode)
    for i in range(size):
        src_file = os.path.join( cluster_dep_folder, os.path.splitext( os.path.basename(output_html_files[frames_to_save[0]]) )[0] + '_' + str(i) + '.' + mode )
        dst_file = os.path.join( second_folder, os.path.splitext( os.path.basename(output_html_files[frames_to_save[0]]) )[0] + '_' + str(i) + '.' + mode )
        subprocess.run('cp '+src_file+' '+dst_file, shell=True, stderr=subprocess.PIPE)
        print('copying from '+src_file+' to '+dst_file+'...')
    bokehplot.save_figs(iframe=frames_to_save[1], path=cluster_dep_folder, mode=mode)
    for i in range(size):
        src_file = os.path.join( cluster_dep_folder, os.path.splitext( os.path.basename(output_html_files[frames_to_save[1]]) )[0] + '_' + str(i) + '.' + mode )
        dst_file = os.path.join( second_folder, os.path.splitext( os.path.basename(output_html_files[frames_to_save[1]]) )[0] + '_' + str(i) + '.' + mode )
        subprocess.run('cp '+src_file+' '+dst_file, shell=True, stderr=subprocess.PIPE)
        print('copying from '+src_file+' to '+dst_file+'...')

def main():
    #load ROOT TTree
    file = up.open( data_path )
    tree = file['tree0']
    executor = concurrent.futures.ThreadPoolExecutor() #executor for parallel data loading
    up_cache = {}
    frame_shift = int(nframes/2)
    
    ######Cluster dependent number of hits#########
    if FLAGS.hits or FLAGS.all:
        print('loading hits data...')
        df_hits = tree.arrays([beamen_str, 'Nhits*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)
        print('done!')

        df_hits_split = []
        for en in beam_energies:
            df_hits_split.append(df_hits[ df_hits[beamen_str] == en])
        del df_hits

        axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': '#hits / cluster'}
        graphs2d(df_hits_split, axis_kwargs_hits, columns_field='Nhits', iframe=0, variable='hits', weight_by_energy=True)
        del df_hits_split
        print('saving the plots...')
        save_plots(frames_to_save=(0,0+frame_shift))
        print('done!')

    ######Cluster dependent energy#################
    if FLAGS.energies or FLAGS.all:
        print('loading energy data...')
        df_energy = tree.arrays([beamen_str, 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)
        print('done!')

        df_energy_split = []
        for en in beam_energies:
            df_energy_split.append(df_energy[ df_energy[beamen_str] == en])    
        del df_energy #clears RAM

        axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Total energy per cluster [MeV]'}
        graphs2d(df_energy_split, axis_kwargs_en, columns_field='Energy', iframe=1, variable='energy', weight_by_energy=False)
        del df_energy_split
        print('saving the plots...')
        save_plots(frames_to_save=(1,1+frame_shift))
        print('done!')
        
    ######Number of clusters######################
    if FLAGS.numbers or FLAGS.all:
        print('loading number data...')
        df_number = tree.arrays([beamen_str, 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True) #'Energy*' is used by len() function to get the number of clusters
        print('done!')

        df_number_split = []
        for en in beam_energies:
            df_number_split.append(df_number[ df_number[beamen_str] == en])
        del df_number

        axis_kwargs_numbers = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusters'}
        graphs2d(df_number_split, axis_kwargs_numbers, columns_field='Energy', variable='number', iframe=2, weight_by_energy=False)
        del df_number_split
        print('saving the plots...')
        save_plots(frames_to_save=(2,2+frame_shift))
        print('done!')

    ######Custer X positions######################
    if FLAGS.posx or FLAGS.all:
        print('loading x positions data...')
        df_posx = tree.arrays([beamen_str, 'X*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)
        print('done!')

        df_posx_split = []
        for en in beam_energies:
            df_posx_split.append(df_posx[ df_posx[beamen_str] == en])
        del df_posx

        axis_kwargs_posx = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' X position"}
        graphs2d(df_posx_split, axis_kwargs_posx, columns_field='X', variable='pos', iframe=3, weight_by_energy=False)
        del df_posx_split
        print('saving the plots...')
        save_plots(frames_to_save=(3,3+frame_shift))
        print('done!')

    ######Custer Y positions######################
    if FLAGS.posy or FLAGS.all:
        print('loading y positions data...')
        df_posy = tree.arrays([beamen_str, 'Y*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)
        print('done!')

        df_posy_split = []
        for en in beam_energies:
            df_posy_split.append(df_posy[ df_posy[beamen_str] == en])
        del df_posy

        axis_kwargs_posy = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' Y position"}
        graphs2d(df_posy_split, axis_kwargs_posy, columns_field='Y', variable='pos', iframe=4, weight_by_energy=False)
        del df_posy_split
        print('saving the plots...')
        save_plots(frames_to_save=(4,4+frame_shift))
        print('done!')

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
    cache_name_energy = 'uproot_cache_energy.pickle'
    cache_file_name_energy = os.path.join(eos_base, cms_user[0], cms_user, data_directory, cache_name_energy)
    cache_name_number = 'uproot_cache_number.pickle'
    cache_file_name_number = os.path.join(eos_base, cms_user[0], cms_user, data_directory, cache_name_number)

    #define analysis constants
    nlayers = 28
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = 2*len(true_beam_energies_GeV)
    assert(len(beam_energies)==size/2)

    print("Input data read from: {}".format(data_path))

    #create output files with plots
    create_dir( os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory)
    output_html_files = ( os.path.join(output_html_dir, 'plot_clusters_hits.html'),
                          os.path.join(output_html_dir, 'plot_clusters_energy_nocut.html'),
                          os.path.join(output_html_dir, 'plot_clusters_number_nocut.html'),
                          os.path.join(output_html_dir, 'plot_clusters_posx_nocut.html'),
                          os.path.join(output_html_dir, 'plot_clusters_posy_nocut.html'),
                          os.path.join(output_html_dir, 'profile_clusters_hits.html'),
                          os.path.join(output_html_dir, 'profile_clusters_energy_nocut.html'),
                          os.path.join(output_html_dir, 'profile_clusters_number_nocut.html'),
                          os.path.join(output_html_dir, 'profile_clusters_posx_nocut.html'),
                          os.path.join(output_html_dir, 'profile_clusters_posy_nocut.html') )
    nframes = len(output_html_files)
    bokehplot = bkp.BokehPlot(filenames=output_html_files, nfigs=nframes*(size,), nframes=nframes)
    plot_width, plot_height = 600, 400
    cluster_dep_folder = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'cluster_dep')
    create_dir( cluster_dep_folder )

    parser = argparse.ArgumentParser()
    FLAGS, _ = add_args(parser, 'clusters')
    main()
