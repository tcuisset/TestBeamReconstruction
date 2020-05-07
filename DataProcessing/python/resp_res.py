import sys
import os
import subprocess
import errno
import glob
import pandas as pd
import numpy as np
from functools import reduce
import bokehplot as bkp
from bokeh.models import Range1d

class ProcessData:
    @staticmethod
    def join(glob_path):
        dfs = []
        for j,i in enumerate(glob.glob(glob_path)):
            df_tmp = pd.read_csv(i)
            dfs.append(df_tmp)
        df_reduced = reduce(lambda left,right: pd.merge(left, right, left_index=True, right_index=True, how='outer')
                            .fillna(-1.), dfs)

        df_beamen = []
        for b in beam_energies:
            #energy sum columns that correspond to a specific beam energy
            #columns with the same energy are put into the same histogram
            cols = [x for x in df_reduced.columns if ( 'ensum' in x and df_reduced.at[1,'beamen'+x[-3:]]==b ) ]
            df_beamen.append( df_reduced.loc[:, cols].stack() if cols != [] else None )
        return df_beamen

    @staticmethod
    def scale_energy(dfs, scale_factor, shift_factor):
        dfs_scaled = []
        for i,df in enumerate(dfs): #each dataframe is a pandas Series
            df_tmp = df.apply(lambda x: scale_factor[i]*x + shift_factor[i] if x!= -1 else x) #perform no scaling and shifting when the value is -1 (fillna used before)
            dfs_scaled.append(df_tmp)
        return dfs_scaled

class HandleHistograms:
    @staticmethod
    def create(dfs, bins, iframe):
        hist = []
        for i,idf in enumerate(dfs):
            if idf is not None:
                idf = idf[idf>-1]
                hist.append( np.histogram(idf, density=False, bins=bins[i], range=(idf.min(),idf.max())) ) #flattens 'df' to one-dimension
                
        axis_kwargs = {'x.axis_label': 'Total RecHit energy per event [MeV]', 'y.axis_label': 'Counts'}
        bokehplot.histogram(data=hist, idx=[x for x in range(len(true_beam_energies_GeV))], iframe=iframe, style='step',
                            legend_label=[str(x)+' GeV' for x in true_beam_energies_GeV], fill_color='white', line_color=line_colors, alpha=0.5,
                            fig_kwargs=axis_kwargs)
        return hist

    @staticmethod
    def fit(hist, parameters, ranges, iframe):
        sigma_units_left, sigma_units_right = 1, 2.5
        common_args = {'pdf':'gaus', 'line_width':2.5, 'alpha':0.8}

        means = []
        means_err = []
        responses = []
        responses_err = []
        resolutions = []
        resolutions_err = []
        for i in range(len(true_beam_energies_GeV)):
            #First fit
            coeff, _ = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=0, iframe=iframe, color=line_colors[i], **common_args)
            fit_bounds = (coeff[1]-sigma_units_left*coeff[2], coeff[1]+sigma_units_right*coeff[2])
            
            #Create and draw amputed histogram 
            #Detail: the last item of the histogram with the counts is removed since it refers to the events between the last
            #edge selected and the next one. Nedges = Nbinswithcounts + 1
            selection = (hist[i][1]>fit_bounds[0]) & (hist[i][1]<fit_bounds[1])
            indexes_selected = np.nonzero(selection) #return indexes of the non-zero (True) elements
            hist_modified = [hist[i][0][indexes_selected][:-1], hist[i][1][selection]]
            bokehplot.histogram(data=hist_modified, idx=i, iframe=iframe, style='step',
                                legend_label=[str(x)+' GeV' for x in true_beam_energies_GeV], fill_color='white', line_color=line_colors[i], 
                                alpha=1., fig_kwargs={'x_range': ranges[i], 'y_range': Range1d(-30.,850)})

            #Second fit
            coeff, var = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=1, iframe=iframe, color=line_colors[i], **common_args)
            err = np.sqrt(np.diag(var))
            err1 = round(err[1],2)
            err2 = round(err[2],2)
            mean_label = 'mean='+str(round(coeff[1],2))+'+-'+str(err1)+' MeV'
            sigma_label = 'sigma='+str(round(coeff[2],2))+'+-'+str(err2)+' MeV'
            font_size = {'text_font_size': '10pt', 'x_units': 'screen'}
            bokehplot.label(mean_label,  idx=i, iframe=iframe, x=10, y=600, **font_size)
            bokehplot.label(sigma_label, idx=i, iframe=iframe, x=10, y=550, **font_size)
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
    axis_kwargs1 = {'t.text': 'Responses after original RecHits calibration',
                    'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response (E/True - 1)'}
    bokehplot.graph(idx=3, iframe=frameid, data=[np.array(true_beam_energies_GeV),np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eresp1)/2,np.array(eresp1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs1)
    bokehplot.graph(idx=3, iframe=frameid, data=[np.array(true_beam_energies_GeV),np.array(resp2)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eresp2)/2,np.array(eresp2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')

    axis_kwargs2 = {'t.text': 'Differences after original RecHits calibration',
                    'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response difference', 'l.location': 'bottom_right'}
    bokehplot.graph(idx=4, iframe=frameid, data=[np.array(true_beam_energies_GeV), np.array(resp2)-np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 ), 
                             np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 )]],
                    style='square', line=True, color='blue', legend_label="Reconstructable - Clusterized", fig_kwargs=axis_kwargs2)
    #fig = bokehplot.get_figure(idx=3, iframe=frameid)
    #fig.legend.location = 'bottom_right'

    axis_kwargs3 = {'t.text': 'Resolutions after original and clusterized RecHits calibrations',
                    'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': u" Fractional energy Resolution (\u03c3 / E)"}
    bokehplot.graph(idx=2, iframe=frameid, data=[np.array(true_beam_energies_GeV),np.array(res1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eres1)/2,np.array(eres1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs3)
    bokehplot.graph(idx=2, iframe=frameid, data=[np.array(true_beam_energies_GeV),np.array(res2)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eres2)/2,np.array(eres2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')


def linear_fit_graph(mean, emean, idx, iframe):
    axis_kwargs = {'t.text': 'Original RecHits calibration' if idx==0 else 'Clusterized RecHits calibration', 
                   'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Reconstructed energy [GeV]'}
    bokehplot.graph(idx=idx, iframe=iframe, 
                    data=[np.array(true_beam_energies_GeV),np.array(mean)/1000], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(emean)/2,np.array(emean)/2]],
                    style='circle', line=False, color='black', fig_kwargs=axis_kwargs)
    coeff, var = bokehplot.fit(pdf='linear', idx=idx, iframe=iframe, p0=([1.,0.]), color='red', 
                               legend_label='y(x) = m*x + b', fig_kwargs={'l.location': 'bottom_right'})
    err = np.sqrt(np.diag(var))
    err1 = round(err[0],2)
    err2 = round(err[1],2)
    pm = (u'\u00B1').encode('utf-8')
    m_label = 'm = '+str(round(coeff[0],2))+pm+str(err1)
    b_label = 'b = '+str(round(coeff[1],2))+pm+str(err2)+' GeV'
    font_size = {'text_font_size': '9pt', 'x_units': 'data'}
    bokehplot.label(m_label, idx=idx, iframe=iframe, x=15, y=280, **font_size)
    bokehplot.label(b_label, idx=idx, iframe=iframe, x=15, y=260, **font_size)
    return coeff[0], coeff[1]

def main():
    #files with sum of rechit energy
    usercode_path = 'src/UserCode/DataProcessing/job_output'
    path_start = '' if ecut_str == '' else 'Ecut'
    path = os.path.join(cmssw_base, usercode_path, 'out'+path_start+'_')
    bins = (1000, 1800, 4200, 5000, 5000, 4200, 5700, 5500, 5500, 500)
    histo_ranges1 = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                    Range1d(64000, 120000), Range1d(88000,130000), Range1d(120000,165000), Range1d(160000, 210000), 
                    Range1d(210000,250000), Range1d(240000,295000))
    histo_ranges2 = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                    Range1d(64000, 120000), Range1d(88000,130000), Range1d(120000,165000), Range1d(160000, 210000), 
                    Range1d(210000,250000), Range1d(240000,295000))

    pars1 = ([750, 20000.,  2000.], #20GeV
             [750, 30000.,  1200.], #30GeV
             [750, 50000.,  2100.], #50GeV
             [750, 80000.,  2500.], #80GeV
             [750, 100000., 3000.], #100GeV
             [750, 110000., 2800.], #120GeV
             [750, 150000., 3000.], #150GeV
             [750, 190000., 3500.], #200GeV
             [750, 215000., 4000.], #250GeV
             [750, 280000., 4500.]) #300GeV 
    data1 = ProcessData.join(path + '*_noclusters.csv')
    hist1 = HandleHistograms.create(data1, bins, iframe=0)
    mean1, emean1, _, _, _, _ = HandleHistograms.fit(hist1, pars1, histo_ranges1, iframe=0)

    last_frame_id = bokehplot.get_nframes() - 1
    calibration_slope1, calibration_shift1 = linear_fit_graph(mean1, emean1, idx=0, iframe=last_frame_id)
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
             [750, 220000., 4000.], #250GeV
             [750, 270000., 4500.]) #300GeV 
    data2 = ProcessData.join(path + '*[0-9][0-9][0-9].csv')
    hist2 = HandleHistograms.create(data2, bins, iframe=2)
    mean2, emean2, _, _, _, _ = HandleHistograms.fit(hist2, pars2, histo_ranges2, iframe=2)

    data2_corrected = ProcessData.scale_energy(data2, correction_value1*np.ones(size), -calibration_shift1*np.ones(size))
    hist2_corrected = HandleHistograms.create(data2_corrected, bins, iframe=3)
    pars2_corrected = tuple([x[0], x[1]/calibration_slope1, x[2]] for x in pars2)
    histo_ranges2_corrected = tuple(Range1d(x.start/calibration_slope1, x.end/calibration_slope1) for x in histo_ranges2)
    _, _, resp2, eresp2, _, _ = HandleHistograms.fit(hist2_corrected, pars2_corrected, histo_ranges2_corrected, iframe=3)

    calibration_slope2, calibration_shift2 = linear_fit_graph(mean2, emean2, idx=1, iframe=last_frame_id)
    correction_value2 = 1 / calibration_slope2
    data2_corrected2 = ProcessData.scale_energy(data2, correction_value2*np.ones(size), -calibration_shift2*np.ones(size))
    hist2_corrected2 = HandleHistograms.create(data2_corrected2, bins, iframe=4)
    pars2_corrected2 = tuple([x[0], x[1]/calibration_slope2, x[2]] for x in pars2)
    histo_ranges2_corrected2 = tuple(Range1d(x.start/calibration_slope2, x.end/calibration_slope2) for x in histo_ranges2)
    _, _, _, _, res2, eres2 = HandleHistograms.fit(hist2_corrected2, pars2_corrected2, histo_ranges2_corrected2, iframe=4)

    response_and_resolution_graphs(resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2, last_frame_id)
    for i in range(len(nfigs)):
        if i==len(nfigs)-1:
            bokehplot.save_frame(iframe=i, plot_width=400, plot_height=400, nrows=2, ncols=3, show=False)
        else:
            bokehplot.save_frame(iframe=i, plot_width=350, plot_height=350, show=False)

if __name__ == '__main__':
    ecut_str = '' if len(sys.argv)==1 else sys.argv[1]

    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    size = len(beam_energies)
    assert(size==len(true_beam_energies_GeV))

    cms_user = subprocess.check_output("echo $USER", shell=True).split('\n')[0]
    data_directory = 'TestBeamReconstruction'
    create_dir( os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory)
    output_html_files = ( os.path.join(output_html_dir, 'pure_rechit_energy'+ecut_str+'.html'),
                          os.path.join(output_html_dir, 'pure_rechit_energy'+ecut_str+'_scaled.html'),
                          os.path.join(output_html_dir, 'clusterized_rechit_energy'+ecut_str+'.html'),
                          os.path.join(output_html_dir, 'clusterized_rechit_energy'+ecut_str+'_scaled_with_original_calibration.html'),
                          os.path.join(output_html_dir, 'clusterized_rechit_energy'+ecut_str+'_scaled_with_clusterized_calibration.html'),
                          os.path.join(output_html_dir, 'summary_graphs'+ecut_str+'.html') 
    )
    nfigs = (size, size, size, size, size, 5)
    bokehplot = bkp.BokehPlot(filenames=output_html_files, nframes=len(nfigs), nfigs=nfigs)
    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
    main()


"""
Errors for 'resp1 / resp2':
[np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 ), 
np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 )]],
"""
