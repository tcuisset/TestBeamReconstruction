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
import argparse
from argparser import add_args
import utils

def graphs_single(df, columns_field, iframe, do_2D=False, energy_index=None):
    """Plots layer-related quantities; `df` refers to one single beam energy"""
    if do_2D == False and energy_index == None:
        raise ValueError('graphs_single: When doing 1D plots, please specify the energy you want to plot.')

    #determine adequate binning
    data_per_layer = []
    if columns_field == 'Densities':
        data_per_layer_seeds = []
        data_per_layer_noseeds = []
        
    datamax, datamin = -np.inf, np.inf
    nhits_allowed = 1
    for ilayer in range(1,nlayers+1):
        layer_col = utils.get_layer_col(df, starts_with=columns_field, ilayer=ilayer)
        arr = df.loc[:, layer_col]
        if arr.size == 0:
            raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))
        arr = utils.flatten_dataframe(arr)
        if columns_field == 'Densities':
            layer_col_en = utils.get_layer_col(df, starts_with='Energies', ilayer=ilayer)
            arr_en = utils.flatten_dataframe( df.loc[:, layer_col_en] )
            layer_col_seeds = utils.get_layer_col(df, starts_with='isSeed', ilayer=ilayer)
            arr_seeds = utils.flatten_dataframe( df.loc[:, layer_col_seeds] )
            layer_col_sizes = utils.get_layer_col(df, starts_with='ClustSize', ilayer=ilayer)
            arr_sizes = utils.flatten_dataframe( df.loc[:, layer_col_sizes] )
            arr_seeds_only = arr[ (arr_en > energy_cut) & (arr_seeds == True) & (arr_sizes >= nhits_allowed) ] #select only seeds with energy > 1000MeV
            arr_no_seeds = arr[ (arr_en > energy_cut) & (arr_seeds == False) & (arr_sizes >= nhits_allowed) ] #select only non-seeds with energy > 1000MeV
            arr = arr[ (arr_en > energy_cut) & (arr_sizes >= nhits_allowed) ] #select only the hits with energy > 1000MeV
            data_per_layer_seeds.append( arr_seeds_only )
            data_per_layer_noseeds.append( arr_no_seeds )
        elif columns_field == 'Distances':
            layer_col_seeds = utils.get_layer_col(df, starts_with='isSeed', ilayer=ilayer)
            arr_seeds = df.loc[:, layer_col_seeds]
            arr_seeds = utils.flatten_dataframe(arr_seeds)
            arr[ arr_seeds == True ] = 0.5 #select CLUE seeds
        data_per_layer.append( arr )
        
        m1, m2 = data_per_layer[-1].min(), data_per_layer[-1].max()
        if m1 < datamin:
            datamin = m1
        if m2 > datamax:
            datamax = m2

    if columns_field == 'Nhits' or columns_field == 'Energy':
        nbins = 60
    else:
        nbins = 20 if datamax>20 else int(datamax-1)
    bins = np.linspace(datamin, datamax, nbins+1) #bins = np.logspace(np.log10(m1), np.log10(m2), nbins+1)
    height = utils.calculate_rect_side_for_plot_bins(bins, scale='linear', nlayers=nlayers)

    #plot data as 2d graphs
    all_counts, all_layers, all_centers = ([] for _ in range(3))
    means, sigmas_l, sigmas_r = ([] for _ in range(3))
    for ilayer in range(1,nlayers+1):
        counts, edges = np.histogram(data_per_layer[ilayer-1], bins=bins)        
        centers = (edges[:-1]+edges[1:])/2
        same_layer_array = ilayer*np.ones(len(centers))
        all_counts.extend(counts.tolist())
        all_layers.extend(same_layer_array.tolist())
        all_centers.extend(centers.tolist())

        #use histogram and not original data to calculate mean and std/sqrt(n)
        #original data cannot be used for the case of weighted histograms
        mean, _ = utils.get_mean_and_sigma(centers, counts)
        sigma_left, sigma_right = utils.get_sigma_band(centers, counts)
        means.append(mean)
        sigmas_l.append(sigma_left)
        sigmas_r.append(sigma_right)

    if columns_field == 'Densities':
        xlabel = 'Density [MeV]'
    elif columns_field == 'Distances':
        xlabel = 'Distance to nearest higher density [cm]'
    elif columns_field == 'Nhits':
        xlabel = 'Fraction of clusterized hits'
    elif columns_field == 'Energy':
        xlabel = 'Fraction of clusterized energy'

    if do_2D:
        print('plotting 2D (i.e., as a function of the layers)...')
        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]),
                      'x.axis_label': 'Layer', 'y.axis_label': xlabel}
        bokehplot.graph(data=[np.array(all_layers), np.array(all_centers), np.array(all_counts)],
                        width=np.ones((len(all_layers))), height=height,
                        idx=energy_index, iframe=iframe, style='rect%Cividis', fig_kwargs=fig_kwargs, alpha=0.6)
        bokehplot.graph(data=[np.arange(1,29), np.array(means)],
                        idx=energy_index, iframe=iframe, color='red', style='circle', size=2, legend_label='mean')

    else:
        print('plotting 1D (i.e., counts in the Y axis) ...')
        fig_kwargs_1D = {'plot_width': plot_width, 'plot_height': plot_height,
                         't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]),
                         'x.axis_label': xlabel, 'y.axis_label': 'Counts'}
        hist, hist2, hist3, peaks, peaks2, peaks3, indexes, leg_labels, leg_labels2, leg_labels3 = ([] for _ in range(10))
        hist_min, hist_max, norm, hist_bins = 1000, 5500, False, 100
        for ilayer in range(nlayers):
            if columns_field == 'Densities':
                fig_kwargs_1D['t.text'] = 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]) + " | E{(hit) > 1000MeV "

                hist_tmp = np.histogram(data_per_layer[ilayer], bins=hist_bins, density=norm, range=(hist_min,hist_max))
                hist_seeds_tmp = np.histogram(data_per_layer_seeds[ilayer], bins=hist_bins, density=norm, range=(hist_min,hist_max))
                hist_noseeds_tmp = np.histogram(data_per_layer_noseeds[ilayer], bins=hist_bins, density=norm, range=(hist_min,hist_max))

                if all(hist_tmp[0] == 0) or all(np.isnan(hist_tmp[0])):
                    hist.append( [np.zeros((hist_tmp[1].size-1)), hist_tmp[1]] )
                    peaks.append( 0. )
                else:
                    hist.append( hist_tmp )
                    peaks.append( utils.peak_abciss(hist_tmp[0], hist_tmp[1]) )

                if all(hist_seeds_tmp[0] == 0) or all(np.isnan(hist_seeds_tmp[0])):
                    hist2.append( [np.zeros((hist_seeds_tmp[1].size-1)), hist_seeds_tmp[1]] )
                    peaks2.append( 0. )
                else:
                    hist2.append( hist_seeds_tmp )
                    peaks2.append( utils.peak_abciss(hist_seeds_tmp[0], hist_seeds_tmp[1]) )

                if all(hist_noseeds_tmp[0] == 0) or all(np.isnan(hist_noseeds_tmp[0])):
                    hist3.append( [np.zeros((hist_noseeds_tmp[1].size-1)), hist_noseeds_tmp[1]] )
                    peaks3.append( 0. )
                else:
                    hist3.append( hist_noseeds_tmp )
                    peaks3.append( utils.peak_abciss(hist_noseeds_tmp[0], hist_noseeds_tmp[1]) )
                
            else:
                hist.append( np.histogram(data_per_layer[ilayer], bins=100, density=True) )
            indexes.append( ilayer )
            leg_labels.append( 'Layer ' + str(ilayer+1) + ', all hits' )
            leg_labels2.append( 'Layer ' + str(ilayer+1) + ', seeds only')
            leg_labels3.append( 'Layer ' + str(ilayer+1) + ', no seeds')

        bokehplot.histogram(data=hist, idx=indexes, legend_label=leg_labels, iframe=iframe, style='step', fig_kwargs=fig_kwargs_1D)
        if columns_field == 'Densities':
            bokehplot.histogram(data=hist2, idx=indexes, legend_label=leg_labels2,
                                iframe=iframe, style='step', color='red', fig_kwargs=fig_kwargs_1D)
            bokehplot.histogram(data=hist3, idx=indexes, legend_label=leg_labels3,
                                iframe=iframe, style='step', color='purple', fig_kwargs=fig_kwargs_1D)
            bokehplot.graph(data=[np.arange(1,29), np.array(peaks)], idx=nlayers, legend_label='All hits', iframe=iframe, style='circle', line=True,
                            fig_kwargs={'t.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]),
                                        'x.axis_label': 'Layer', 'y.axis_label': 'Maximum density [MeV]', 'l.location': 'bottom_left'})
            bokehplot.graph(data=[np.arange(1,29), np.array(peaks2)], idx=nlayers, legend_label='Seeds only', iframe=iframe, style='square', color='red', line=True)
            bokehplot.graph(data=[np.arange(1,29), np.array(peaks3)], idx=nlayers, legend_label='No seeds', iframe=iframe, style='triangle', color='purple', line=True)

        #plotting additional vertical lines
        if columns_field == 'Densities':
            pass
            #bokehplot.line(x=[[sigmaNoiseTimesKappa,sigmaNoiseTimesKappa] for _ in range(nlayers)], 
            #               y=[[0,hist[i][0].max()] for i in range(nlayers)], 
            #               idx=indexes, iframe=iframe, color='orange', legend_label='approximate kappa cut')
        elif columns_field == 'Distances':
            bokehplot.line(x=[[1.3,1.3] for _ in range(nlayers)], 
                           y=[[0,hist[i][0].max()] for i in range(nlayers)], 
                           idx=indexes, iframe=iframe, color='orange', legend_label='1.3cm')

    del all_counts
    del all_layers
    del all_centers
    del means
    del sigmas_l
    del sigmas_r

def graphs_double(df, mode, iframe, energy_index=2): #default to energy_index=2 (=50GeV)
    """Plots density vs distance or posx vs poxy (the latter weighted by the energy density). The variables are defined according to CLUE."""
    if mode == 'dens_dist':
        columns_fields = ('Densities', 'Distances', 'isSeed')
        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]),
                      'x.axis_label': 'Density [MeV]', 'y.axis_label': 'Distance to nearest higher density [cm]'}
        nbins_max_x = 35
        nbins_max_y = 8
    elif mode == 'pos':
        columns_fields = ('PosX', 'PosY', 'Densities')
        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'Beam energy: {} GeV'.format(true_beam_energies_GeV[energy_index]),
                      'x.axis_label': 'X hit position [cm]', 'y.axis_label': 'Y hit position [cm]'}
        nbins_max_x = 30
        nbins_max_y = 30
    else:
        raise ValueError('Wrong mode!')
    len_col_fields = len(columns_fields)
    nbins = (nbins_max_x, nbins_max_y)

    for ilayer in range(1,nlayers+1):
        datamax, datamin = [-np.inf,-np.inf], [np.inf,np.inf]
        layer_col = tuple( utils.get_layer_col(df, starts_with=columns_fields[x], ilayer=ilayer) for x in range(len_col_fields) )
        arr = [ df.loc[:, layer_col[x]] for x in range(len_col_fields) ]
        for x in range(len(columns_fields)):
            if arr[x].size == 0:
                raise ValueError('The dataframe in layer {} of variable {} is empty!'.format(ilayer, x))
        arr = [utils.flatten_dataframe(arr[x]) for x in range(len_col_fields) ]
        if mode == 'dens_dist':
            arr_selection_seeds = (arr[2] == True) #seeds
            arr[1][arr_selection_seeds] = 0.5
            data_weights = arr[0] ** 6
        elif mode == 'pos':
            arr[0] = arr[0][ arr[2] > energy_cut ]
            arr[1] = arr[1][ arr[2] > energy_cut ]
            data_weights = arr[2][ arr[2] > energy_cut ]

        m1 = (arr[0].min(), arr[1].min()) 
        m2 = (arr[0].max(), arr[1].max()) 
        if m1[0] < datamin[0]:
            datamin[0] = m1[0]
        if m2[0] > datamax[0]:
            datamax[0] = m2[0]
        if m1[1] < datamin[1]:
            datamin[1] = m1[1]
        if m2[1] > datamax[1]:
            datamax[1] = m2[1]

        if mode == 'dens_dist':
            bins = ( np.linspace(datamin[0], datamax[0], nbins[0]+1),
                     np.linspace(datamin[1], datamax[1], nbins[1]+1) )
        else:
            limit = 5
            bins = ( np.linspace(-limit, limit, nbins[0]+1),
                     np.linspace(-limit, limit, nbins[1]+1) )

        counts, xedges, yedges = np.histogram2d(x=arr[0], y=arr[1], bins=bins, weights=data_weights, density=False)
        bokehplot.histogram(data=(counts, xedges, yedges),
                            idx=ilayer-1, iframe=iframe, style='quad%blues', fig_kwargs=fig_kwargs,
                            continuum_value=0, continuum_color='white')
        layer_label = 'Layer '+str(ilayer)
        label_props = {'text_font_size': '10pt', 'x_units': 'screen', 'y_units': 'screen',
                       'border_line_color': 'black', 'background_fill_color': 'LightCyan'}
        bokehplot.label(layer_label, idx=ilayer-1, iframe=iframe, x=340, y=290, **label_props)
        #bokehplot.box(x=[sigmaNoiseTimesKappa,xmaxlimit], y=[1.3, datamax[1]], idx=ilayer-1, color='red', line_width=3)

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
            print('WARNING: The cache was reset.')
            self.cache = {} 
        return self.cache

def main():

    chosen_energy = 50 #density vs densities plots will only refer to this energy (one plot per layer)
    assert(chosen_energy in beam_energies)
    
    for iEn,thisEn in enumerate(beam_energies):
        
        #load ROOT TTree for a specific energy
        file = up.open( data_paths[iEn] )
        tree = file['tree0']

        #load cache for a specific energy
        cacheobj = CacheManager( cache_file_names[iEn] )
        up_cache = cacheobj.load()

        if FLAGS.densities_distances:
            if thisEn == chosen_energy:
                print('Loading density and distance data for {}GeV...'.format(beam_energies[iEn]))
                df = tree.arrays(['Distances*', 'Densities*', 'isSeed*', 'Energies*'], outputtype=pd.DataFrame, cache=up_cache)
                index = beam_energies.index(chosen_energy)
                graphs_double(df, mode='dens_dist', iframe=output_html_files_map['densities_distances'][1], energy_index=index )

        if FLAGS.densities or FLAGS.densities_2D:
            if thisEn == chosen_energy:
                if not FLAGS.densities_distances: #avoid loading the same data again
                    print('Loading density data for {}GeV...'.format(beam_energies[iEn]))
                    df = tree.arrays(['Densities*', 'Energies*', 'isSeed*', 'ClustSize*'], outputtype=pd.DataFrame, cache=up_cache)
                index = beam_energies.index(chosen_energy)
                if FLAGS.densities:
                    graphs_single(df, columns_field='Densities', energy_index=index, iframe=output_html_files_map['densities'][1], do_2D=False)
            if FLAGS.densities_2D:
                graphs_single(df, columns_field='Densities', energy_index=index, iframe=output_html_files_map['densities_2D'][1], do_2D=True)
                        
        if FLAGS.distances or FLAGS.distances_2D:
            if thisEn == chosen_energy:
                if not FLAGS.densities_distances: #avoid loading the same data again
                    print('Loading distance data for {}GeV...'.format(beam_energies[iEn]))
                    df = tree.arrays(['Distances*', 'isSeed*'], outputtype=pd.DataFrame, cache=up_cache)
                index = beam_energies.index(chosen_energy)
                if FLAGS.distances:
                    graphs_single(df, columns_field='Distances', energy_index=index, iframe=output_html_files_map['distances'][1], do_2D=False)
            if FLAGS.distances_2D:
                graphs_single(df, columns_field='Distances', energy_index=index, iframe=output_html_files_map['distances_2D'][1], do_2D=True)

        if FLAGS.densities or FLAGS.distances or FLAGS.densities_2D or FLAGS.distances_2D or FLAGS.densities_distances:
            if thisEn == chosen_energy:
                del df

        if FLAGS.hits_fraction:
            print('Loading hits fraction data for {}GeV...'.format(beam_energies[iEn]))
            df = tree.arrays(['Nhits*'], outputtype=pd.DataFrame, cache=up_cache)
            graphs_single(df, columns_field='Nhits', iframe=output_html_files_map['hits_fraction'][1], energy_index=iEn, do_2D=True)

        if FLAGS.energy_fraction:
            print('Loading energy fraction data for {}GeV...'.format(beam_energies[iEn]))
            df = tree.arrays(['EnergyFrac*'], outputtype=pd.DataFrame, cache=up_cache)
            graphs_single(df, columns_field='Energy', iframe=output_html_files_map['energy_fraction'][1], energy_index=iEn, do_2D=True)

        if FLAGS.posx_posy:
            if thisEn == chosen_energy:
                print('Loading X and Y position data for {}GeV...'.format(beam_energies[iEn]))
                df = tree.arrays(['Densities*', 'PosX*', 'PosY*'], outputtype=pd.DataFrame, cache=up_cache)
                index = beam_energies.index(chosen_energy)
                graphs_double(df, mode='pos', iframe=output_html_files_map['posx_posy'][1], energy_index=index )


        cacheobj.dump() #dump cache for specific energy
    
    print('Saving BokehPlot frames...')
    layer_dep_folder = os.path.join(eos_base, cms_user[0], cms_user, 'www', analysis_directory, 'layer_dep', FLAGS.datatype, FLAGS.showertype)
    utils.create_dir( layer_dep_folder )
    presentation_path = os.path.join(home, release, 'DN/figs', 'layer_dep', FLAGS.datatype, FLAGS.showertype)
    utils.create_dir( presentation_path )
    bokehplot.save_frames(plot_width=plot_width, plot_height=plot_height, show=False)
    for iframe in range(nframes):
        utils.create_dir( layer_dep_folder )
        bokehplot.save_figs(iframe=iframe, path=layer_dep_folder, mode='png')
        #bokehplot.save_figs(iframe=iframe, path=presentation_path, mode='png')

if __name__ == '__main__':
    #define analysis constants
    nlayers = 28 if FLAGS.showertype=='em' else 40
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    energy_cut = 1000
    sigmaNoiseTimesKappa = 9 * 10. / 6.
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(true_beam_energies_GeV)
    assert(len(beam_energies)==size)

    #define parser for user input arguments
    parser = argparse.ArgumentParser()
    FLAGS, _ = add_args(parser, 'layers')
    utils.input_sanity_checks(FLAGS, sys.argv)
    
    #define local data paths and variables
    eos_base = '/eos/user/'
    cms_user = subprocess.check_output("echo $USER", shell=True, encoding='utf-8').split('\n')[0]
    release = 'CMSSW_11_1_0_pre2/src/'
    home = subprocess.check_output(b'echo $HOME', shell=True, encoding='utf-8').split('\n')[0]
    analysis_directory = 'TestBeamReconstruction/'
    data_directory = 'job_output/layer_dependent/'
    data_path_start = os.path.join(eos_base, cms_user[0], cms_user, analysis_directory, data_directory)
    data_paths = [os.path.join(data_path_start, 'hadd_layerdep_' + FLAGS.datatype + '_' + FLAGS.showertype + '_beamen' + str(x) + '.root') for x in beam_energies]

    #define cache names and paths
    cache_file_name_start = os.path.join(eos_base, cms_user[0], cms_user, analysis_directory)
    cache_file_names = [os.path.join(cache_file_name_start, 'uproot_cache_layerdep_beamen' + str(x) + '.root') for x in beam_energies]

    print("Input data:")
    for x in data_paths:
        print(x)

    #create output files with plots
    output_html_dir = os.path.join(eos_base, cms_user[0], cms_user, 'www', analysis_directory, 'layer_dep', FLAGS.datatype, FLAGS.showertype)
    utils.create_dir( output_html_dir )
    plot_width, plot_height = 600, 400
    
    #the keys are the attributes of FLAGS (except 'all')
    #the values are: 1) the name of all potential bokehplot frames, 2) number of bokehplot figures in each frame       
    outlambda = lambda x: os.path.join(output_html_dir, FLAGS.datatype + '_' + FLAGS.showertype + x)
    output_html_potential_files_map = { 'densities':           ( outlambda('_densities_1D.html'),       nlayers + 1 ), 
                                        #one additional plot to show the max vs. layer 
                                        'distances':           ( outlambda('_distances_1D.html'),       nlayers), 
                                        'densities_distances': ( outlambda('_dens_vs_dist.html'),       nlayers ),
                                        'posx_posy':           ( outlambda('_posx_vs_posy.html'),       nlayers ),
                                        'densities_2D':        ( outlambda('_densities_2D.html'),       size    ),
                                        'distances_2D':        ( outlambda('_distances_2D.html'),       size    ),
                                        'hits_fraction':       ( outlambda('_hits_fraction_2D.html'),   size    ),
                                        'energy_fraction':     ( outlambda('_energy_fraction_2D.html'), size    )  }

    #select only the frames requested by the user
    #the values will become: 1) the name of all potential bokehplot frames, 2) a frame identifier, 3) number of bokehplot figures in each frame       
    counter = 0
    output_html_files_map = dict()
    for k,tup in output_html_potential_files_map.items():
        if getattr(FLAGS,k) == True:
            output_html_files_map.update({k: (tup[0],counter,tup[1])})
            counter += 1

    output_html_files_list = [tup[0] for k,tup in output_html_files_map.items()]
    nfigs = [tup[2] for k,tup in output_html_files_map.items()]
    nframes = len(output_html_files_list)
    bokehplot = bkp.BokehPlot(filenames=output_html_files_list, nfigs=nfigs, nframes=nframes)
    main()
