import sys
import os
import warnings
import subprocess
import errno
import glob
import argparse
import pandas as pd
import numpy as np
import h5py
from argparser import AddArgs
from functools import reduce
import bokehplot as bkp
from bokeh.models import Range1d
import utils

class ProcessData:
    @staticmethod
    def join(glob_path):
        dfs = []
        for j,i in enumerate(glob.glob(glob_path)):
            df_tmp = pd.read_csv(i)
            dfs.append(df_tmp)
        red = reduce(lambda left,right: pd.merge(left, right, left_index=True, right_index=True, how='outer')
                            .fillna(-1.), dfs)
        df_beamen = []
        for b in beam_energies:
            #energy sum columns that correspond to a specific beam energy
            #columns with the same energy are put into the same histogram
            cols = [x for x in red.columns if ( 'ensum' in x and red.iat[1,red.columns.get_loc(x)+1]==b ) ]
            df_beamen.append( red.loc[:, cols].stack() if cols != [] else None )
        return df_beamen

    @staticmethod
    def scale_energy(dfs, scale_factor, shift_factor):
        dfs_scaled = []
        for i,df in enumerate(dfs): #each dataframe is a pandas Series
            if df is not None:
                df_tmp = df.apply(lambda x: scale_factor[i]*x + shift_factor[i] if x!= -1 else x) #perform no scaling and shifting when the value is -1 (fillna used before)
                dfs_scaled.append(df_tmp)
            else:
                dfs_scaled.append(None)
        return dfs_scaled

class HandleHistograms:
    @staticmethod
    def create(dfs, bins, iframe):
        hist = []
        axis_kwargs = {'x.axis_label': 'Total RecHit energy per event [MeV]', 'y.axis_label': 'Counts'}
        for i,idf in enumerate(dfs):
            if idf is not None:
                idf = idf[idf>-1.]
                h = np.histogram(idf, density=False, bins=bins[i], range=(idf.min(),idf.max())) #flattens 'df' to one-dimension
                bokehplot.histogram(data=h, idx=i, iframe=iframe, style='step',
                                    legend_label=str(true_beam_energies_GeV[i])+' GeV', line_color=line_colors, alpha=0.5,
                                    fig_kwargs=axis_kwargs)
                hist.append(h)
            else:
                hist.append(None)

        return hist

    @staticmethod
    def fit(hist, parameters, ranges, iframe):
        sigma_units_left, sigma_units_right = 1, 2.5
        common_args = {'pdf': 'gaus', 'line_width': 2.5, 'alpha': 0.8}

        means = []
        means_err = []
        responses = []
        responses_err = []
        resolutions = []
        resolutions_err = []
        for i in range(len(true_beam_energies_GeV)):
            #First fit
            if hist[i] is None:
                print('WARNING: Missing dataset for {} GeV.'.format(beam_energies[i]))
                means.append( 0. )
                means_err.append( 0. )
                responses.append( 0. ) 
                responses_err.append( 0. )
                resolutions.append( 0. )
                resolutions_err.append( 0. )
                continue;
            elif len(hist[i][0])==0 or len(hist[i][1])==0: #no data => no fit!
                print('WARNING: There appears to be no data available for {} GeV in fit #1.'.format(beam_energies[i]))
                hist[i] = None

            if hist[i] is not None:
                try:
                    coeff, _ = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=0, iframe=iframe, 
                                             color=line_colors[i], line_dash='dashed', **common_args)
                    fit_bounds = (coeff[1]-sigma_units_left*coeff[2], coeff[1]+sigma_units_right*coeff[2])
                except (RuntimeError, ValueError) as e:
                    fit_bounds = (hist[i][1][0], hist[i][1][-1])
            else:
                fit_bounds = (hist[i][1][0], hist[i][1][-1])

            #Create and draw amputed histogram 
            #Detail: the last item of the histogram with the counts is removed since it refers to the events between the last
            #edge selected and the next one. Nedges = Nbinswithcounts + 1
            selection = (hist[i][1]>fit_bounds[0]) & (hist[i][1]<fit_bounds[1])
            indexes_selected = np.nonzero(selection) #remove last edge to match with number of histogram entries (#bins)

            hist_modified = [hist[i][0][indexes_selected[0][:-1]], hist[i][1][selection]]
            if hist[i] is None or len(hist_modified[0])==1: #no data => no fit!
                hist[i] = None
                print('WARNING: There appears to be no data available for {} GeV in fit#2'.format(beam_energies[i]))
            else:
                bokehplot.histogram(data=hist_modified, idx=i, iframe=iframe, style='step',
                                    legend_label=[str(x)+' GeV' for x in true_beam_energies_GeV], line_color=line_colors[i], 
                                    #alpha=1., fig_kwargs={'x_range': ranges[i], 'y_range': Range1d(-30.,900)})
                                    alpha=1., fig_kwargs={'x_range': ranges[i]})

            #Second fit
            if hist[i] is not None:
                try:
                    coeff, var = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=1, iframe=iframe, 
                                               color=line_colors[i], **common_args)
                    err = np.sqrt(np.diag(var))
                    err1 = round(err[1],2)
                    err2 = round(err[2],2)
                except (RuntimeError, ValueError):
                    coeff = [0., 0., 0.]
                    err1 = 0
                    err2 = 0
                mean_label = 'mean='+str(round(coeff[1],2))+'+-'+str(err1)+' MeV'
                sigma_label = 'sigma='+str(round(coeff[2],2))+'+-'+str(err2)+' MeV'
                font_size = {'text_font_size': '10pt', 'x_units': 'screen', 'y_units': 'screen'}
                bokehplot.label(mean_label,  idx=i, iframe=iframe, x=10, y=310, **font_size)
                bokehplot.label(sigma_label, idx=i, iframe=iframe, x=10, y=290, **font_size)
            else: #no data => no fit!
                coeff = [0., 0., 0.]
                err1 = 0
                err2 = 0
            means.append( coeff[1] )
            means_err.append( err1 )
            responses.append( (coeff[1]-true_beam_energies_MeV[i]) / (true_beam_energies_MeV[i]) ) 
            responses_err.append( err1 / true_beam_energies_MeV[i] )
            resolutions.append( coeff[2] / true_beam_energies_MeV[i] )
            resolutions_err.append( err2 / true_beam_energies_MeV[i] )
        return means, means_err, responses, responses_err, resolutions, resolutions_err

def create_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def response_and_resolution_graphs(resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2, frameid):
    """Plots responses and resolutions with their errors"""
    c_resp1, c_eresp1, c_en = clean_zeros(resp1, eresp1)
    c_resp2, c_eresp2, _ = clean_zeros(resp2, eresp2)
    c_res1, c_eres1, _ = clean_zeros(res1, eres1)
    c_res2, c_eres2, _ = clean_zeros(res2, eres2)

    if FLAGS.showertype == 'em':
        locations = {'data': 'bottom_right', 'sim_proton': 'top_right', 'sim_noproton': 'bottom_right'}
    elif FLAGS.showertype == 'had':
        locations = {'data': 'bottom_right', 'sim_proton': 'bottom_right'}
    axis_kwargs1 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response (E/True - 1)', 'l.location': locations[FLAGS.datatype]}
    bokehplot.graph(idx=1, iframe=frameid, data=[c_en,c_resp1], 
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))],[c_eresp1/2,c_eresp1/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs1)
    bokehplot.graph(idx=1, iframe=frameid, data=[c_en,c_resp2], 
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))],
                            [c_eresp2/2,c_eresp2/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')

    #'t.text': 'Differences after original RecHits calibration',
    axis_kwargs2 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response difference', 'l.location': 'bottom_right'}
    bokehplot.graph(idx=2, iframe=frameid, data=[c_en, c_resp2-c_resp1], 
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))],
                            [np.sqrt(np.power(c_eresp1,2)+np.power(c_eresp2,2)/2), np.sqrt(np.power(c_eresp1,2)+np.power(c_eresp2,2)/2)]],
                    style='square', line=True, color='blue', legend_label="Reconstructable - Clusterized", fig_kwargs=axis_kwargs2)
    #fig = bokehplot.get_figure(idx=3, iframe=frameid)
    #fig.legend.location = 'bottom_right'

    axis_kwargs3 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': u" Fractional energy Resolution (\u03c3 / E)"}
    bokehplot.graph(idx=0, iframe=frameid, data=[c_en,c_res1], 
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))],[c_eres1/2,c_eres1/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs3)
    bokehplot.graph(idx=0, iframe=frameid, data=[c_en,c_res2], 
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))], [c_eres2/2,c_eres2/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')


def clean_zeros(mean, emean):
    mean, emean = np.array(mean), np.array(emean)
    idx1, idx2 = mean!=0, emean!=0
    assert(np.array_equal(idx1,idx2))
    mean, emean = mean[idx1], emean[idx1]
    en = np.array(true_beam_energies_GeV)[idx1]
    return mean, emean, en

def linear_fit_graph(mean, emean, idx, iframe):
    axis_kwargs = {'t.text': 'Original RecHits calibration' if idx==0 else 'Clusterized RecHits calibration', 
                   'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Reconstructed energy [GeV]'}
    c_mean, c_emean, c_en = clean_zeros(mean, emean)
    bokehplot.graph(idx=idx, iframe=iframe, 
                    data=[c_en, c_mean/1000], #convert to GeV
                    errors=[[np.zeros(len(c_en)),np.zeros(len(c_en))], [c_emean/2,c_emean]],
                    style='circle', line=False, color='black', fig_kwargs=axis_kwargs)
    if len(c_mean) > 1: #at least two points are required for the fit
        coeff, var = bokehplot.fit(pdf='linear', idx=idx, iframe=iframe, p0=([1.,0.]), color='red', 
                                   legend_label='y(x) = m*x + b', fig_kwargs={'l.location': 'bottom_right'})
        err = np.sqrt(np.diag(var))
        err1 = round(err[0],2)
        err2 = round(err[1],2)
        pm = u'\u00B1'
        m_label = 'm = '+str(round(coeff[0],2))+pm+str(err1)
        b_label = 'b = '+str(round(coeff[1],2))+pm+str(err2)+' GeV'
        font_size = {'text_font_size': '9pt', 'x_units': 'data', 'y_units': 'data'}
        bokehplot.label(m_label, idx=idx, iframe=iframe, x=15, y=280, **font_size)
        bokehplot.label(b_label, idx=idx, iframe=iframe, x=15, y=260, **font_size)
    else:
        coeff = [1, 1]
    return coeff[0], coeff[1]

def analyze_data():
    #files with sum of rechit energy
    usercode_path = 'src/UserCode/DataProcessing/job_output'

    #difference due to historic reasons; this will have to be removed if the analysis step is rerun
    path = os.path.join(eos_base, cms_user[0], cms_user, data_directory, FLAGS.tag,
                        'hit_dependent/outEcut_' + FLAGS.datatype + "_" + FLAGS.showertype)

    if FLAGS.showertype == 'em':
        bins = (1000, 1800, 4200, 5000, 5000, 4200, 5700, 5500, 5500, 500)
        histo_ranges1 = (Range1d(0,30000), Range1d(0, 40000), Range1d(27000, 58000), Range1d(52000, 94000), 
                         Range1d(64000, 120000), Range1d(88000,135000), Range1d(120000,165000), Range1d(170000, 220000), 
                         Range1d(200000,280000), Range1d(240000,320000))
        histo_ranges2 = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                         Range1d(64000, 120000), Range1d(88000,135000), Range1d(120000,165000), Range1d(160000, 220000), 
                         Range1d(200000,280000), Range1d(240000,315000))
        pars1 = ([750, 20000.,  2000.], #20GeV
                 [750, 30000.,  1200.], #30GeV
                 [750, 50000.,  2100.], #50GeV
                 [750, 80000.,  2500.], #80GeV
                 [750, 100000., 3000.], #100GeV
                 [750, 110000., 2800.], #120GeV
                 [750, 150000., 3000.], #150GeV
                 [750, 190000., 3500.], #200GeV
                 [750, 245000., 4000.], #250GeV
                 [750, 290000., 4500.]) #300GeV
        
    elif FLAGS.showertype == 'had':
        bins = (1000, 1000, 1000, 1000, 1000, 1000, 1000, 750, 600, 500)
        if FLAGS.datatype == 'data':
            histo_ranges1 = (Range1d(0,30000), Range1d(2000, 40000), Range1d(20000, 70000), Range1d(40000, 105000), 
                             Range1d(50000, 130000), Range1d(60000,150000), Range1d(120000,165000), Range1d(100000, 250000), 
                             Range1d(120000,330000), Range1d(180000,380000))
        else:
            histo_ranges1 = (Range1d(0,30000), Range1d(2000, 40000), Range1d(20000, 58000), Range1d(40000, 94000), 
                             Range1d(58000, 117000), Range1d(70000,135000), Range1d(120000,165000), Range1d(170000, 220000), 
                             Range1d(200000,280000), Range1d(260000,350000))
        histo_ranges2 = histo_ranges1
        pars1 = ([750, 15000.,  2000.], #20GeV
                 [750, 24000.,  1200.], #30GeV
                 [750, 40000.,  2100.], #50GeV
                 [750, 75000.,  2500.], #80GeV
                 [750, 95000., 3000.], #100GeV
                 [750, 115000., 2800.], #120GeV
                 [750, 150000., 3000.], #150GeV
                 [750, 195000., 3500.], #200GeV
                 [750, 250000., 4000.], #250GeV
                 [750, 310000., 4500.]) #300GeV

    data1 = ProcessData.join(path + '*_noclusters.csv')

    hist1 = HandleHistograms.create(data1, bins, iframe=0)
    mean1, emean1, _, _, _, _ = HandleHistograms.fit(hist1, pars1, histo_ranges1, iframe=0)

    last_frame_id = bokehplot.get_nframes() - 1
    calibration_slope1, calibration_shift1 = linear_fit_graph(mean1, emean1, idx=0, iframe=last_frame_id-1)
    if calibration_slope1 == 0:
        print('WARNING: slope #1 is zero!')
        calibration_slope1 = 1;
    correction_value1 = 1 / calibration_slope1

    data1_corrected = ProcessData.scale_energy(data1, correction_value1*np.ones(size), -calibration_shift1*np.ones(size))
    hist1_corrected = HandleHistograms.create(data1_corrected, bins, iframe=1)
    pars1_corrected = tuple([x[0], x[1]/calibration_slope1, x[2]] for x in pars1)
    histo_ranges1_corrected = tuple(Range1d(x.start/calibration_slope1, x.end/calibration_slope1) for x in histo_ranges1)
    _, _, resp1, eresp1, res1, eres1 = HandleHistograms.fit(hist1_corrected, pars1_corrected, histo_ranges1_corrected, iframe=1)

    #files with sum of clusterized rechit energy
    pars2 = ([750, 18000.,  2000.], #20GeV
             [750, 25000.,  1200.], #30GeV
             [750, 43000.,  2000.], #50GeV
             [750, 77000.,  2500.], #80GeV
             [750, 98000.,  2700.], #100GeV
             [750, 109000., 2750.], #120GeV
             [750, 140000., 3500.], #150GeV
             [750, 185000., 3500.], #200GeV
             [750, 240000., 4000.], #250GeV
             [750, 290000., 4500.]) #300GeV 
    data2 = ProcessData.join(path + '*[0-9].csv')
    hist2 = HandleHistograms.create(data2, bins, iframe=2)
    mean2, emean2, _, _, _, _ = HandleHistograms.fit(hist2, pars2, histo_ranges2, iframe=2)
    
    data2_corrected = ProcessData.scale_energy(data2, correction_value1*np.ones(size), -calibration_shift1*np.ones(size))
    hist2_corrected = HandleHistograms.create(data2_corrected, bins, iframe=3)
    pars2_corrected = tuple([x[0], x[1]/calibration_slope1, x[2]] for x in pars2)
    histo_ranges2_corrected = tuple(Range1d(x.start/calibration_slope1, x.end/calibration_slope1) for x in histo_ranges2)
    _, _, resp2, eresp2, _, _ = HandleHistograms.fit(hist2_corrected, pars2_corrected, histo_ranges2_corrected, iframe=3)

    calibration_slope2, calibration_shift2 = linear_fit_graph(mean2, emean2, idx=1, iframe=last_frame_id-1)
    if calibration_slope2 == 0:
        print('WARNING: slope #2 is zero!')
        calibration_slope2 = 1
    correction_value2 = 1 / calibration_slope2

    #if FLAGS.datatype == 'data': #REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #    correction_value2 = 1
    #    calibration_slope2 = 1
    
    data2_corrected2 = ProcessData.scale_energy(data2, correction_value2*np.ones(size), -calibration_shift2*np.ones(size))
    hist2_corrected2 = HandleHistograms.create(data2_corrected2, bins, iframe=4)
    pars2_corrected2 = tuple([x[0], x[1]/calibration_slope2, x[2]] for x in pars2)
    histo_ranges2_corrected2 = tuple(Range1d(x.start/calibration_slope2, x.end/calibration_slope2) for x in histo_ranges2)
    _, _, _, _, res2, eres2 = HandleHistograms.fit(hist2_corrected2, pars2_corrected2, histo_ranges2_corrected2, iframe=4)

    save_folder = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'resp_res', FLAGS.datatype, FLAGS.showertype, FLAGS.tag)
    utils.create_dir( save_folder )
    presentation_path = os.path.join(home, release, 'DN/figs', 'resp_res', FLAGS.datatype)
    print(presentation_path)
    utils.create_dir( presentation_path )
    print(save_folder)
    print(presentation_path)
    for i in range(len(nfigs)-1):
       bokehplot.save_frame(iframe=i, plot_width=plot_width, plot_height=plot_width, show=False)
       bokehplot.save_figs(iframe=i, path=save_folder, mode='png')
       bokehplot.save_figs(iframe=i, path=presentation_path, mode='png')

    #Write data in the HDF5 format
    variables_to_store = [resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2]
    with h5py.File(h5filename, 'w') as hf:
        for name,var in zip(variables_created,variables_to_store):
            hf.create_dataset(name, data=var)

def final_plots():
    variables_stored = []
    hf = h5py.File(h5filename, 'r') 
    for name in variables_created:
        variables_stored.append( hf.get(name) )

    save_folder = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'resp_res', FLAGS.datatype, FLAGS.tag)
    utils.create_dir( save_folder )
    presentation_path = os.path.join(home, release, 'DN/figs', 'resp_res', FLAGS.datatype)
    utils.create_dir( presentation_path )
    frameid = bokehplot.get_nframes()-1
    response_and_resolution_graphs(*variables_stored, frameid=frameid)

    bokehplot.save_frame(iframe=frameid, plot_width=plot_width, plot_height=plot_width, nrows=1, ncols=3, show=False)
    bokehplot.save_figs(iframe=frameid, path=save_folder, mode='png')
    bokehplot.save_figs(iframe=frameid, path=presentation_path, mode='png')

    hf.close()

def main():
    if FLAGS.analyze_only: #analysis and some plotting (store data into HDF5)
        analyze_data()
    elif FLAGS.plot_only: #final resolution and response plots (read data from HDF5)
        final_plots()
    else:
        analyze_data()
        final_plots()

if __name__ == '__main__':
    FLAGS, _ = AddArgs().resp_res()
    utils.input_sanity_checks(FLAGS, sys.argv)

    #beam energies used
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(beam_energies)
    assert(size==len(true_beam_energies_GeV))
    plot_width, plot_height = 400, 400

    #local CMSSW variables and paths
    eos_base = '/eos/user'
    cms_user = subprocess.check_output(b'echo $USER', shell=True, encoding='utf-8').split('\n')[0]
    release = 'TestBeamAnalysis/src/' #'/' + subprocess.check_output(b'echo $CMSSW_VERSION', shell=True, encoding='utf-8').split('\n')[0] + '/src/'
    home = subprocess.check_output(b'echo $HOME', shell=True, encoding='utf-8').split('\n')[0]
    data_directory = 'TestBeamReconstruction'

    output_html_dir = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'resp_res', FLAGS.datatype, FLAGS.showertype, FLAGS.tag)
    outlambda = lambda x: os.path.join(output_html_dir, FLAGS.datatype + '_' + FLAGS.showertype + '_' + x)
    output_html_files = ( outlambda('pure_rechit_energy_Ecut.html'),
                          outlambda('pure_rechit_energy_Ecut_scaled.html'),
                          outlambda('clusterized_rechit_energy_Ecut.html'),
                          outlambda('clusterized_rechit_energy_Ecut_scaled_with_original_calibration.html'),
                          outlambda('clusterized_rechit_energy_Ecut_scaled_with_clusterized_calibration.html'),
                          outlambda('linear_regressions.html'),
                          outlambda('responses_and_resolutions.html') 
    )
    nfigs = [size, size, size, size, size, 2, 3]
    bokehplot = bkp.BokehPlot(filenames=output_html_files, nframes=len(nfigs), nfigs=nfigs,
                              nwidgets=[0 for _ in range(len(nfigs))])
    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']

    #HDF5 data file related variables
    h5filename = os.path.join( eos_base, cms_user[0], cms_user, data_directory, FLAGS.tag, FLAGS.datatype + "_" + FLAGS.showertype + '_' + os.path.splitext( os.path.basename(__file__) )[0] + '.h5')
    variables_created = ('resp1', 'eresp1', 'res1', 'eres1', 'resp2', 'eresp2', 'res2', 'eres2')

    main()
