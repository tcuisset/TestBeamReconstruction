import sys
import os
import pickle
import subprocess
import warnings
import bokehplot as bkp
import uproot as up
import numpy as np
import pandas as pd
import concurrent.futures
import argparse
from tqdm import tqdm
from argparser import add_args
from scipy.interpolate import UnivariateSpline
import utils
from utils import CacheManager

def graphs_per_energy(df, axis_kwargs, columns_field, iframe, variable, energy_index, weight_by_energy=False):
    """Plots cluster-related quantities per layer"""
    if variable not in ('hits', 'energy', 'number', 'pos'):
        raise ValueError('graphs_per_energy: Variable {} is not supported.'.format(variable))

    dl = [ [] for _ in range(ncuts) ] #data per layer
    
    datamax, datamin = [ -np.inf for _ in range(ncuts) ], [ np.inf for _ in range(ncuts) ] #the second list is for the results with cuts
    for ilayer in range(1,nlayers+1):
        #determine adequate binning
        col_arr = df.loc[ :, utils.get_layer_col(df, starts_with=columns_field, ilayer=ilayer) ]
        if col_arr.size == 0:
            raise ValueError('The dataframe {} in layer {} is empty!'.format(i, ilayer))

        for thisCut in range(ncuts):
            if variable == 'number':
                dl[thisCut].append( col_arr.apply( lambda x: len( x[columns_field+'_layer'+str(ilayer)][ x[columns_field+'_layer'+str(ilayer)] > energy_cuts[thisCut] ] ), axis=1, result_type='reduce') )
                dl[thisCut][-1] = dl[thisCut][-1][ dl[thisCut][-1] != 0 ] #filter out all events without clusters
            elif variable == 'energy':
                arr = utils.flatten_dataframe(col_arr)
                arr = arr[ arr > energy_cuts[thisCut] ]
                dl[thisCut].append( arr )
            elif variable == 'hits' or variable == 'pos':
                extra_layer_col = utils.get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                extra_arr = df.loc[:, extra_layer_col]
                extra_arr = utils.flatten_dataframe(extra_arr)
                arr = utils.flatten_dataframe(col_arr)[ extra_arr > energy_cuts[thisCut] ]
                dl[thisCut].append( arr )

            if len(dl[thisCut][-1]) != 0:
                m1, m2 = dl[thisCut][-1].min(), dl[thisCut][-1].max()
            if m1 < datamin[thisCut]:
                datamin[thisCut] = m1
            if m2 > datamax[thisCut]:
                datamax[thisCut] = m2 if variable != 'energy' else 30000 #fixed y axis facilitates the comparison between energies
                
    if variable == 'pos':
        nbins = tuple( 50 for x in range(ncuts) )
    else:
        nbins = tuple( 20 if datamax[x]>20 else int(datamax[x]-1) for x in range(ncuts) )
    bins = tuple( np.linspace(datamin[x], datamax[x], nbins[x]+1) for x in range(ncuts) )
    height_hits = tuple( utils.height_for_plot_bins(bins[x], scale='linear', nlayers=nlayers) for x in range(ncuts) )
    all_counts, all_layers, all_centers = ([[] for _ in range(ncuts)] for _ in range(3))
    means, sigmas_l, sigmas_r = ( [[] for _ in range(ncuts)] for _ in range(3) ) #the second list is for the results with cuts

    #plot data as 2d graphs
    for ilayer in range(1,nlayers+1):
        for thisCut in range(ncuts):
            if weight_by_energy:
                layer_col_weights = utils.get_layer_col(df, starts_with='Energy', ilayer=ilayer)
                weights_array = df.loc[:, layer_col_weights]
                data_weights = utils.flatten_dataframe( weights_array )
                data_weights = data_weights[ data_weights > energy_cuts[thisCut] ]
                if data_weights.size == 0:
                    if thisCut==0:
                        warnings.warn('The weights dataframe with energy cut {} in layer {} is empty!'.format(energy_cuts[thisCut], ilayer))
                    counts, edges = np.zeros(nbins[thisCut]), np.linspace(datamin[thisCut],datamax[thisCut],nbins[thisCut]+1)
                else:
                    counts, edges = np.histogram(dl[thisCut][ilayer-1], bins=bins[thisCut], weights=data_weights)
            else:
                counts, edges = np.histogram(dl[thisCut][ilayer-1], bins=bins[thisCut], range=(-1.,1.))

            centers = (edges[:-1]+edges[1:])/2
            assert(len(counts) == len(centers))
            same_layer_array = ilayer*np.ones(len(centers))
            all_counts[thisCut].extend(counts.tolist())
            all_layers[thisCut].extend(same_layer_array.tolist())
            all_centers[thisCut].extend(centers.tolist())

            #use histogram and not original data to calculate mean and std/sqrt(n)
            #original data cannot be used for the case of weighted histograms
            mean, _ = utils.get_mean_and_sigma(centers, counts)
            sigma_left, sigma_right = utils.get_sigma_band(centers, counts)
            means[thisCut].append(mean)
            sigmas_l[thisCut].append(sigma_left)
            sigmas_r[thisCut].append(sigma_right)

    dl.clear()
    del dl

    layers_x = np.arange(1,nlayers+1)
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

    means_sp_func   = tuple( UnivariateSpline(x=layers_x, y=np.array(means[x]), k=interp_degree, s=interp_smooth) for x in range(ncuts) )
    means_sp        = tuple( means_sp_func[x](layers_x_fine) for x in range(ncuts) )
    sigmasl_sp_func = tuple( UnivariateSpline(x=layers_x, y=np.array(sigmas_l[x]), k=interp_degree, s=interp_smooth) for x in range(ncuts) )
    sigmasl_sp      = tuple( sigmasl_sp_func[x](layers_x_fine) for x in range(ncuts) )
    sigmasr_sp_func = tuple( UnivariateSpline(x=layers_x, y=np.array(sigmas_r[x]), k=interp_degree, s=interp_smooth) for x in range(ncuts) )
    sigmasr_sp      = tuple( sigmasr_sp_func[x](layers_x_fine) for x in range(ncuts) )

    fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height}
    fig_kwargs.update(axis_kwargs)

    for thisCut in range(ncuts):
        tmp_fig_kwargs = { 't.text': '{} GeV beam energy | E(cluster) > {} MeV'.format(true_beam_energies_GeV[energy_index], energy_cuts[thisCut]) }
        tmp_fig_kwargs.update(fig_kwargs)
        bokehplot.graph(data=[np.array(all_layers[thisCut]), np.array(all_centers[thisCut]), np.array(all_counts[thisCut])],
                        width=np.ones((len(all_layers[thisCut]))), height=height_hits[thisCut],
                        idx=energy_index + thisCut*size_shift, iframe=iframe, style='rect%Cividis', fig_kwargs=tmp_fig_kwargs, alpha=0.6)
        bokehplot.graph(data=[layers_x, means[thisCut]],
                        idx=energy_index + thisCut*size_shift, iframe=iframe, color='brown', style='circle', size=2, legend_label='mean')

        bokehplot.graph(data=[layers_x_fine, means_sp[thisCut]],
                        width=interp_thickness*np.ones(len(layers_x_fine)), height=abs(sigmasr_sp[thisCut]-sigmasl_sp[thisCut]),
                        idx=energy_index + thisCut*size_shift, iframe=iframe+frame_shift, color='grey', style='rect', alpha=0.6, legend_label=u'std mean error (\u03c3/\u221an)', 
                        fig_kwargs=tmp_fig_kwargs)  #pm = (u'\u00B1').encode('utf-8')
        bokehplot.graph(data=[layers_x_fine, means_sp[thisCut]],
                        idx=energy_index + thisCut*size_shift, iframe=iframe+frame_shift, color='red', style='circle', size=1.5, legend_label='mean')

def graphs_per_layer(tree, cache, axis_kwargs, iframe, variable, energy_index=2, weight_by_energy=False):
    """Plots cluster-related quantities per layer"""
    if variable not in ('pos', 'spatialres_x', 'spatialres_y', 'spatialres_xy'):
        raise ValueError('graphs_per_layer: Variable {} is not supported.'.format(variable))

    hdf5_name_ = 'cluster_' + FLAGS.showertype + '_' + FLAGS.datatype + '_' + str(FLAGS.chosen_energy) + 'GeV_' + FLAGS.tag
    hdf5_name = os.path.join(data_path_start, hdf5_name_ + '.h5' )
    print(hdf5_name)
    compression_level = 9

    executor = concurrent.futures.ThreadPoolExecutor()
    nbins = 200 if variable == 'pos' else 100 if variable == 'spatialres_xy' else 500
    limit_up, limit_down = 4, 2
    fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height}
    fig_kwargs.update(axis_kwargs)

    #store data for faster retrieval at a later time
    if not FLAGS.use_saved_data:
        emptyfracs = []

        for ilayer in tqdm(range(1,nlayers+1)):
            nx = 'X_layer'+str(ilayer)
            ny = 'Y_layer'+str(ilayer) 
            ndx = 'dX_layer'+str(ilayer)
            ndy = 'dY_layer'+str(ilayer)
            ne = 'Energy_layer'+str(ilayer)
            nn = 'Nhits_layer'+str(ilayer)

            df_ = tree.arrays([nx, ny, ndx, ndy, ne, nn], outputtype=pd.DataFrame, flatten=True, executor=executor, blocking=True, cache=cache)
            nhits_selection = (df_[nn] >= nhits_min) & (df_[nn] < nhits_max)
            df_ = df_[nhits_selection]
            df_layer = df_
            df_layer.to_hdf(hdf5_name, key='layer'+str(ilayer), complevel=compression_level, mode='r+' if ilayer != 1 else 'w')
            
            df_layer_out = df_[ nhits_selection & (df_[nx] < -98) & (df_[ny] < -98)]
            emptyfracs.append(df_layer_out.shape[0] / df_.shape[0])

        emptyfracs = pd.Series(emptyfracs)
        emptyfracs.to_hdf(hdf5_name, key='fracs', complevel=compression_level, mode='r+')


    print('Converting data to plot format...')
    emptyfracs = pd.read_hdf(hdf5_name, 'fracs')
    bias, ebias, res, eres = ([] for _ in range(4))
    for ilayer in tqdm(range(1,nlayers+1)):
        df_layer = pd.read_hdf(hdf5_name, 'layer'+str(ilayer))
        nx = 'X_layer'+str(ilayer) 
        ny = 'Y_layer'+str(ilayer) 
        ndx = 'dX_layer'+str(ilayer)
        ndy = 'dY_layer'+str(ilayer)
        ne = 'Energy_layer'+str(ilayer)
        nn = 'Nhits_layer'+str(ilayer)

        for thisCut in range(ncuts):
            if 'spatialres' in variable and thisCut != 0:
                continue

            if variable == 'pos':
                text = '{} GeV beam energy | E(cluster) > {} MeV | Layer {}'.format(true_beam_energies_GeV[energy_index], energy_cuts[thisCut], ilayer)
            else:
                text = '{} GeV beam energy | Layer {}'.format(true_beam_energies_GeV[energy_index], ilayer)
            tmp_fig_kwargs = {'t.text': text}
            tmp_fig_kwargs.update(fig_kwargs)

            df_layer_cut = df_layer[ df_layer[ne] > energy_cuts[thisCut] ]
            #df_layer_cut = df_layer_cut[ (df_layer_cut[nn] >= 5) & (df_layer_cut[nx]<-2.3) ]
            df_layer_cut = df_layer_cut[ (df_layer_cut[nn] >= 5) ]

            if variable == 'pos':
                if weight_by_energy:
                    bokehplot.histogram( np.histogram2d( x=df_layer_cut[nx].to_numpy(), y=df_layer_cut[ny].to_numpy(), bins=nbins,
                                                         weights=df_layer_cut[ne].to_numpy(), range=[[-limit_up,limit_down],[-limit_down,limit_up]] ),
                                         idx=ilayer-1+thisCut*nlayers, iframe=iframe, style='quad%Cividis', fig_kwargs=tmp_fig_kwargs)
                else:
                    bokehplot.histogram( np.histogram2d(x=df_layer_cut[nx].to_numpy(), y=df_layer_cut[ny].to_numpy(), bins=nbins, range=[[-limit_up,limit_down],[-limit_down,limit_up]]),
                                         idx=ilayer-1+thisCut*nlayers, iframe=iframe, style='quad%Cividis', fig_kwargs=tmp_fig_kwargs)

            elif variable == 'spatialres_x':
                bokehplot.histogram( np.histogram(df_layer_cut[ndx].to_numpy(), bins=nbins, range=[-1.,1.]), 
                                     idx=ilayer-1, iframe=iframe, color='green', fig_kwargs=tmp_fig_kwargs)
                common_args = {'pdf': 'gaus', 'line_width': 2.5, 'alpha': 0.8}
                coeff, var = bokehplot.fit(p0=[2000, 0.,  .2], idx=ilayer-1, iframe=iframe, obj_idx=0, color='black', line_dash='dashed', debug=False, **common_args)

            elif variable == 'spatialres_y':
                bokehplot.histogram( np.histogram(df_layer_cut[ndy].to_numpy(), bins=nbins, range=[-1.,1.]), 
                                     idx=ilayer-1+thisCut*nlayers, iframe=iframe, color='blue', fig_kwargs=tmp_fig_kwargs )
                common_args = {'pdf': 'gaus', 'line_width': 2.5, 'alpha': 0.8}
                coeff, var = bokehplot.fit(p0=[2000, 0.,  .2], idx=ilayer-1, iframe=iframe, obj_idx=0, color='black', line_dash='dashed', debug=False, **common_args)
            elif variable == 'spatialres_xy':
                bokehplot.histogram( np.histogram2d(x=df_layer_cut[ndx].to_numpy(), y=df_layer_cut[ndy].to_numpy(), bins=nbins, range=[[-1.,1.],[-1.,1.]]), idx=ilayer-1, iframe=iframe, style='quad%CET_L17', fig_kwargs=tmp_fig_kwargs )
                kwargs_other = {'x.axis_label': "#hits / cluster", 'y.axis_label': "hit energy [MeV]",
                                't.text': text, 'plot_width': plot_width, 'plot_height': plot_height }
                bokehplot.histogram( np.histogram2d(x=df_layer_cut[nn].to_numpy(), y=df_layer_cut[ne].to_numpy(), bins=(15,100), range=[[0.1,15.],[0.1,200]]), idx=ilayer-1+nlayers, iframe=iframe, style='quad%CET_L17', scale='log', fig_kwargs=kwargs_other )

            if variable == 'spatialres_x' or variable == 'spatialres_y':
                bias.append( coeff[1] )
                res.append( coeff[2] )
                err = np.sqrt(np.diag(var))
                ebias.append( round(err[1],3) )
                eres.append( round(err[2],3) )

        del df_layer

    if variable == 'pos':
        emptyfrac_kwargs = {'t.text': '{} GeV beam energy'.format(true_beam_energies_GeV[energy_index]),
                            'x.axis_label': 'Layer', 'y.axis_label': 'Fraction of clusters were the position measurement failed',
                            'plot_width': plot_width, 'plot_height': plot_height}
        bokehplot.histogram( (np.array(emptyfracs), np.arange(.5,nlayers+1)), 
                             idx=ncuts*nlayers, iframe=iframe, color='orange',
                             fig_kwargs=emptyfrac_kwargs)
    elif variable == 'spatialres_x' or variable == 'spatialres_y':
        extra = variable[-1]
        hdf5_name_summ = os.path.join(data_path_start, hdf5_name_ + '_d' + variable[-1] + '_summary.h5' )
        pd.Series(bias).to_hdf(hdf5_name_summ, key='bias'+extra, complevel=compression_level, mode='w')
        pd.Series(ebias).to_hdf(hdf5_name_summ, key='ebias'+extra, complevel=compression_level, mode='r+')
        pd.Series(res).to_hdf(hdf5_name_summ, key='res'+extra, complevel=compression_level, mode='r+')
        pd.Series(eres).to_hdf(hdf5_name_summ, key='eres'+extra, complevel=compression_level, mode='r+')

def save_plots(frame_key):
    mode = 'png'
    presentation_folder = os.path.join(home, release, 'DN/figs', 'cluster_dep', FLAGS.datatype, FLAGS.showertype)
    utils.create_dir( presentation_folder )

    #save frames
    bokehplot.save_frame(iframe=output_html_files_map[frame_key][1], plot_width=plot_width, plot_height=plot_height, show=False)
    #if version2 in output_html_files_map.keys():
    #    bokehplot.save_frame(iframe=output_html_files_map[frame_key+version2][1], plot_width=plot_width, plot_height=plot_height, show=False)

    #save figs
    bokehplot.save_figs(iframe=output_html_files_map[frame_key][1], path=cluster_dep_folder, mode=mode)
    for i in range(size):
        src_file = os.path.join( cluster_dep_folder, os.path.splitext( os.path.basename(output_html_files_map[frame_key][0]) )[0] + '_' + str(i) + '.' + mode )
        dst_file = os.path.join( presentation_folder, os.path.splitext( os.path.basename(output_html_files_map[frame_key][0]) )[0] + '_' + str(i) + '.' + mode )
        print('copying from '+src_file+' to '+dst_file+'...')
        subprocess.run('cp '+src_file+' '+dst_file, shell=True, stderr=subprocess.PIPE)

    #save version2 figs
    if version2 in output_html_files_map.keys():
        bokehplot.save_figs(iframe=output_html_files_map[frame_key+version2][1], path=cluster_dep_folder, mode=mode)
        for i in range(size):
            src_file = os.path.join( cluster_dep_folder, os.path.splitext( os.path.basename(output_html_files_map[frame_key + version2][0]) )[0] + '_' + str(i) + '.' + mode )
            dst_file = os.path.join( presentation_folder, os.path.splitext( os.path.basename(output_html_files_map[frame_key + version2][0]) )[0] + '_' + str(i) + '.' + mode )
            print('copying from '+src_file+' to '+dst_file+'...')
            subprocess.run('cp '+src_file+' '+dst_file, shell=True, stderr=subprocess.PIPE)

def main():
    executor = concurrent.futures.ThreadPoolExecutor() #executor for parallel data loading

    for iEn,thisEn in enumerate(beam_energies):

        #check the file exists
        if not os.path.isfile(data_paths[iEn]):
            warnings.warn('Missing dataset for {} GeV.'.format(thisEn))
            continue

        #load ROOT TTree
        file = up.open( data_paths[iEn] )
        tree = file['tree0']
        
        #load cache for a specific energy
        cacheobj = CacheManager( cache_file_names[iEn] )
        up_cache = cacheobj.load()

        ######Cluster dependent number of hits#########
        if FLAGS.hits:
            print('Loading hits data for {}GeV...'.format(beam_energies[iEn]))
            df_hits = tree.arrays(['Nhits*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)
            axis_kwargs_hits = {'x.axis_label': 'Layer', 'y.axis_label': '#hits / cluster'}
            graphs_per_energy(df_hits, axis_kwargs_hits, columns_field='Nhits', iframe=output_html_files_map['hits'][1], variable='hits', energy_index=iEn, weight_by_energy=True)
            del df_hits

        ######Cluster dependent energy#################
        if FLAGS.energies:
            print('Loading energy data for {}GeV...'.format(beam_energies[iEn]))
            df_energy = tree.arrays(['Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)

            axis_kwargs_en = {'x.axis_label': 'Layer', 'y.axis_label': 'Total energy per cluster [MeV]'}
            graphs_per_energy(df_energy, axis_kwargs_en, columns_field='Energy', iframe=output_html_files_map['energies'][1], variable='energy', energy_index=iEn, weight_by_energy=False)
            del df_energy

        ######Number of clusters######################
        if FLAGS.numbers:
            print('Loading number data for {}GeV...'.format(beam_energies[iEn]))
            df_number = tree.arrays(['Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True) #'Energy*' is used by len() function to get the number of clusters

            axis_kwargs_numbers = {'x.axis_label': 'Layer', 'y.axis_label': 'Number of clusters'}
            graphs_per_energy(df_number, axis_kwargs_numbers, columns_field='Energy', variable='number', iframe=output_html_files_map['numbers'][1], energy_index=iEn, weight_by_energy=False)
            del df_number

        ######Custer X positions######################
        if FLAGS.posx:
            print('Loading and storing X positions data for {}GeV...'.format(beam_energies[iEn]))
            df_posx = tree.arrays(['X*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)

            axis_kwargs_posx = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' X position"}
            graphs_per_energy(df_posx, axis_kwargs_posx, columns_field='X', variable='pos', iframe=output_html_files_map['posx'][1], energy_index=iEn, weight_by_energy=True)
            del df_posx

        ######Custer Y positions######################
        if FLAGS.posy:
            print('Loading and storing Y positions data for {}GeV...'.format(beam_energies[iEn]))
            df_posy = tree.arrays(['Y*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)

            axis_kwargs_posy = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' Y position"}
            graphs_per_energy(df_posy, axis_kwargs_posy, columns_field='Y', variable='pos', iframe=output_html_files_map['posy'][1], energy_index=iEn, weight_by_energy=True)
            del df_posy

        ######Custer X vs Cluster Y positions######################
        if FLAGS.posx_posy and thisEn == FLAGS.chosen_energy:
            print('Loading and storing X and Y positions data for {}GeV...'.format(beam_energies[iEn]))

            axis_kwargs_posxy = {'x.axis_label': "Clusters' X position", 'y.axis_label': "Clusters' Y position"}
            graphs_per_layer(tree, up_cache, axis_kwargs_posxy, variable='pos', iframe=output_html_files_map['posx_posy'][1], energy_index=iEn, weight_by_energy=True)

        ######Custer 1D X spatial resolution######################
        if FLAGS.dx and thisEn == FLAGS.chosen_energy:
            print('Loading and storing X spatial resolution data for {}GeV...'.format(beam_energies[iEn]))

            axis_kwargs_dx = {'x.axis_label': "Clusters' X spatial resolution [cm]", 'y.axis_label': "Counts"}
            graphs_per_layer(tree, up_cache, axis_kwargs_dx, variable='spatialres_x', iframe=output_html_files_map['dx'][1], energy_index=iEn, weight_by_energy=False)

        ######Custer 1D X spatial resolution######################
        if FLAGS.dy and thisEn == FLAGS.chosen_energy:
            print('Loading and storing X spatial resolution data for {}GeV...'.format(beam_energies[iEn]))

            axis_kwargs_dy = {'x.axis_label': "Clusters' Y spatial resolution [cm]", 'y.axis_label': "Counts"}
            graphs_per_layer(tree, up_cache, axis_kwargs_dy, variable='spatialres_y', iframe=output_html_files_map['dy'][1], energy_index=iEn, weight_by_energy=False)

        ######Custer 1D X spatial resolution######################
        if FLAGS.dx_dy and thisEn == FLAGS.chosen_energy:
            print('Loading and storing X and Y spatial resolution data for {}GeV...'.format(beam_energies[iEn]))

            axis_kwargs_dx = {'x.axis_label': "Clusters' X spatial resolution [cm]", 'y.axis_label': "Clusters' Y spatial resolution [cm]"}
            graphs_per_layer(tree, up_cache, axis_kwargs_dx, variable='spatialres_xy', iframe=output_html_files_map['dx_dy'][1], energy_index=iEn, weight_by_energy=False)

        ######Custer 2D X spatial resolution######################
        if FLAGS.dx_2D:
            print('Loading and storing 2D X spatial resolution data for {}GeV...'.format(beam_energies[iEn]))
            df_dx2D = tree.arrays(['dX*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)

            axis_kwargs_dx2D = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' X spatial resolution"}
            graphs_per_energy(df_dx2D, axis_kwargs_dx2D, columns_field='dX', variable='pos', iframe=output_html_files_map['dx_2D'][1], energy_index=iEn, weight_by_energy=False)
            del df_dx2D

        ######Custer Y spatial resolution######################
        if FLAGS.dy_2D:
            print('Loading and storing Y spatial resolution data for {}GeV...'.format(beam_energies[iEn]))
            df_dy2D = tree.arrays(['dY*', 'Energy*'], outputtype=pd.DataFrame, cache=up_cache, executor=executor, blocking=True)

            axis_kwargs_dy2D = {'x.axis_label': 'Layer', 'y.axis_label': "Clusters' Y spatial resolution"}
            graphs_per_energy(df_dy2D, axis_kwargs_dy2D, columns_field='dY', variable='pos', iframe=output_html_files_map['dy_2D'][1], energy_index=iEn, weight_by_energy=False)
            del df_dy2D

    print('Saving the plots...')
    for k in output_html_files_map.keys():
        if version2 not in k:
            save_plots(frame_key=k)

if __name__ == '__main__':
    #define parser for user input arguments
    parser = argparse.ArgumentParser()
    FLAGS, _ = add_args(parser, 'clusters')
    utils.input_sanity_checks(FLAGS, sys.argv)

    #define analysis constants
    nlayers = 28 if FLAGS.showertype=='em' else 40
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    energy_cuts = (0,1000)
    ncuts = len(energy_cuts)
    size = ncuts*len(true_beam_energies_GeV)
    size_shift = len(true_beam_energies_GeV)
    version2 = '_prof'
    nhits_min, nhits_max = 1, 10000
    assert(FLAGS.chosen_energy in beam_energies)

    #define local data paths and variables
    eos_base = '/eos/user/'
    cms_user = subprocess.check_output("echo $USER", shell=True, encoding='utf-8').split('\n')[0]
    release = 'CMSSW_11_1_0_pre2/src/' #subprocess.check_output(b'echo $CMSSW_VERSION', shell=True, encoding='utf-8').split('\n')[0] + '/src/'
    home = subprocess.check_output(b'echo $HOME', shell=True, encoding='utf-8').split('\n')[0]
    data_directory = 'TestBeamReconstruction'
    data_path_start = os.path.join(utils.data_path, FLAGS.tag, 'cluster_dependent')
    data_paths = [os.path.join(data_path_start, 'hadd_clusterdep_' + FLAGS.datatype + '_' + FLAGS.showertype + '_beamen' + str(x) + '.root') for x in beam_energies]
        
    #define cache names and paths
    cache_file_name_start = utils.data_path
    cache_file_names = [os.path.join(cache_file_name_start, 'uproot_cache_clusterdep_beamen' + str(x) + '.root') for x in beam_energies]

    utils.print_input_data(data_paths)

    #create output files with plots
    utils.create_dir( os.path.join(utils.plot_html_path, 'cluster_dep', FLAGS.datatype, FLAGS.showertype, FLAGS.tag) )
    output_html_dir = os.path.join(utils.plot_html_path, 'cluster_dep', FLAGS.datatype, FLAGS.showertype, FLAGS.tag)
    outlambda = lambda x: os.path.join(output_html_dir, FLAGS.datatype + '_' + FLAGS.showertype + '_' + str(FLAGS.chosen_energy) + x)
    output_html_files_potential_map = { 'hits':                ( outlambda('_plot_clusters_hits.html'),         size ),
                                        'energies':            ( outlambda('_plot_clusters_energy.html'),       size ),
                                        'numbers':             ( outlambda('_plot_clusters_number.html'),       size ),
                                        'posx':                ( outlambda('_plot_clusters_posx.html'),         size ),
                                        'posy':                ( outlambda('_plot_clusters_posy.html'),         size ),
                                        'hits' + version2:     ( outlambda(version2 + '_clusters_hits.html'),   size ),
                                        'energies' + version2: ( outlambda(version2 + '_clusters_energy.html'), size ),
                                        'numbers' + version2:  ( outlambda(version2 + '_clusters_number.html'), size ),
                                        'posx' + version2:     ( outlambda(version2 + '_clusters_posx.html'),   size ),
                                        'posy' + version2:     ( outlambda(version2 + '_clusters_posy.html'),   size ),
                                        'posx_posy':           ( outlambda('_plot_clusters_posx_vs_posy.html'.format(FLAGS.chosen_energy)), ncuts*nlayers + 1 ),
                                        'dx':                  ( outlambda('_plot_clusters_dx_nowindow.html'),  nlayers),
                                        'dy':                  ( outlambda('_plot_clusters_dy.html'),           nlayers),
                                        'dx_dy':               ( outlambda('_plot_clusters_dxdy.html'),         2*nlayers ),
                                        'dx_2D':               ( outlambda('_plot_clusters_dx2D.html'),         size ),
                                        'dy_2D':               ( outlambda('_plot_clusters_dy2D.html'),         size )}
    counter = 0
    output_html_files_map = dict()
    for k,tup in output_html_files_potential_map.items():
        if version2 not in k:
            if getattr(FLAGS,k) == True:
                output_html_files_map.update({k: (tup[0],counter,tup[1])})
                counter += 1
        else:
            if getattr(FLAGS,k[:-5]) == True:
                output_html_files_map.update({k: (tup[0],counter,tup[1])})
                counter += 1
                
    output_html_files_list = [tup[0] for k,tup in output_html_files_map.items()]
    nfigs = [tup[2] for k,tup in output_html_files_map.items()]
    nframes = len(output_html_files_list)
    frame_shift = int(nframes/2)
    
    bokehplot = bkp.BokehPlot(filenames=output_html_files_list, nfigs=nfigs, nframes=nframes)
    plot_width, plot_height = 600, 400
    cluster_dep_folder = os.path.join(utils.plot_html_path, 'cluster_dep', FLAGS.datatype, FLAGS.showertype, FLAGS.tag)
    utils.create_dir( cluster_dep_folder )

    main()
