import sys
import os
import subprocess
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
    def shift_energy(dfs, means):
        assert(len(means) == len(true_beam_energies_GeV))
        assert(len(means) == len(dfs))
        dfs_shifted = []
        for i,df in enumerate(dfs): #each dataframe is a pandas Series
            shift_factor = true_beam_energies_MeV[i]/means[i]
            df_tmp = df.apply(lambda x: shift_factor*x if x!= -1 else x) #perform no shifting when the value is -1 (fillna used before)
            dfs_shifted.append(df_tmp)
        return dfs_shifted

class HandleHistograms:
    @staticmethod
    def create(dfs, bins):
        hist = []
        for i,idf in enumerate(dfs):
            if idf is not None:
                idf = idf[idf>-1]
                hist.append( np.histogram(idf, density=False, bins=bins[i], range=(idf.min(),idf.max())) ) #flattens 'df' to one-dimension
                
        axis_kwargs = {'x.axis_label': 'Total RecHit energy per event [MeV]', 'y.axis_label': 'Counts'}
        bokehplot.histogram(data=hist, idx=[x for x in range(len(true_beam_energies_GeV))], style='step',
                            legend_label=[str(x)+' GeV' for x in true_beam_energies_GeV], fill_color='white', line_color=line_colors, alpha=0.5,
                            fig_kwargs=axis_kwargs)
        return hist

    @staticmethod
    def fit(hist, parameters, ranges):
        sigma_units = 3
        common_args = {'pdf':'gaus', 'line_width':2.5, 'alpha':0.8}

        means = []
        means_err = []
        responses = []
        responses_err = []
        resolutions = []
        resolutions_err = []
        for i in range(len(true_beam_energies_GeV)):
            #First fit
            coeff, _ = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=0, color=line_colors[i], **common_args)
            fit_bounds = (coeff[1]-sigma_units*coeff[2], coeff[1]+sigma_units*coeff[2])
            
            #Create and draw amputed histogram 
            #Detail: the last item of the histogram with the counts is removed since it refers to the events between the last
            #edge selected and the next one. Nedges = Nbinswithcounts + 1
            selection = (hist[i][1]>fit_bounds[0]) & (hist[i][1]<fit_bounds[1])
            indexes_selected = np.nonzero(selection) #return indexes of the non-zero (True) elements
            hist_modified = [hist[i][0][indexes_selected][:-1], hist[i][1][selection]]
            bokehplot.histogram(data=hist_modified, idx=i, style='step',
                                legend_label=[str(x)+' GeV' for x in true_beam_energies_GeV], fill_color='white', line_color=line_colors[i], 
                                alpha=1., fig_kwargs={'x_range': ranges[i]})

            #Second fit
            coeff, var = bokehplot.fit(p0=parameters[i], idx=i, obj_idx=1, color=line_colors[i], **common_args)
            err = np.sqrt(np.diag(var))
            err1 = round(err[1],2)
            err2 = round(err[2],2)
            mean_label = 'mean='+str(round(coeff[1],2))+'+-'+str(err1)+' MeV'
            sigma_label = 'sigma='+str(round(coeff[2],2))+'+-'+str(err2)+' MeV'
            font_size = {'text_font_size': '8pt'}
            bokehplot.label(mean_label,  idx=i, x=10, y=320, **font_size)
            bokehplot.label(sigma_label, idx=i, x=10, y=305, **font_size)
            means.append( coeff[1] )
            means_err.append( err1 )
            responses.append( (coeff[1]-true_beam_energies_MeV[i]) / (true_beam_energies_MeV[i]) ) 
            resolutions.append( coeff[2] / true_beam_energies_MeV[i] )
            responses_err.append( err1 / true_beam_energies_MeV[i] )
            resolutions.append( coeff[2] / true_beam_energies_MeV[i] )
            resolutions_err.append( err2 / true_beam_energies_MeV[i] )
        bokehplot.show_frame(plot_width=300, plot_height=300)
        return means, means_err, responses, responses_err, resolutions, resolutions_err

            
def final_graphs(resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2, frameid):
    """Plots responses and resolutions with their errors"""
    axis_kwargs1 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response (E/True - 1)'}
    bokehplot.graph(data=[np.array(true_beam_energies_GeV),np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eresp1)/2,np.array(eresp1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs1)
    bokehplot.graph(data=[np.array(true_beam_energies_GeV),np.array(resp2)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eresp2)/2,np.array(eresp2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')

    axis_kwargs2 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response difference'}
    bokehplot.graph(idx=1, data=[np.array(true_beam_energies_GeV), np.array(resp2)-np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 ), 
                             np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 )]],
                    style='square', line=True, color='blue', legend_label="Reconstructable - Clusterized", fig_kwargs=axis_kwargs2)
    fig = bokehplot.get_figure(idx=1, iframe=frameid)
    fig.legend.location = 'bottom_right'

    axis_kwargs3 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': u" Resolution (\u03c3 / E) [MeV]"}
    bokehplot.graph(idx=2, data=[np.array(true_beam_energies_GeV),np.array(res1)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eres1)/2,np.array(eres1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs3)
    bokehplot.graph(idx=2, data=[np.array(true_beam_energies_GeV),np.array(res2)], 
                    errors=[[np.zeros(len(true_beam_energies_GeV)),np.zeros(len(true_beam_energies_GeV))],
                            [np.array(eres2)/2,np.array(eres2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')
    bokehplot.show_frame(plot_width=300, plot_height=300)

def main():
    #files with sum of rechit energy
    usercode_path = 'src/UserCode/DataProcessing/job_output'
    path = os.path.join(cmssw_base, usercode_path, 'out_')
    average_shift = 1.05
    bins = (1000, 1800, 4200, 5000, 5000, 4200, 5700, 5500, 5500, 500)
    histo_ranges1 = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                    Range1d(64000, 120000), Range1d(88000,130000), Range1d(120000,165000), Range1d(160000, 210000), 
                    Range1d(210000,250000), Range1d(240000,295000))
    histo_ranges2 = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                    Range1d(64000, 120000), Range1d(88000,130000), Range1d(120000,165000), Range1d(160000, 210000), 
                    Range1d(210000,250000), Range1d(240000,295000))
    histo_ranges1_shifted = tuple(Range1d(average_shift*x.start, average_shift*x.end) for x in histo_ranges1)
    histo_ranges2_shifted = tuple(Range1d(average_shift*x.start, average_shift*x.end) for x in histo_ranges2)

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
    hist1 = HandleHistograms.create(data1, bins)
    mean1, emean1, resp1, eresp1, _, _ = HandleHistograms.fit(hist1, pars1, histo_ranges1)

    #files with sum of clusterized rechit energy
    bokehplot.add_frame('clusterized_rechit_energy.html', nfigs=len(true_beam_energies_GeV))
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
    hist2 = HandleHistograms.create(data2, bins)
    mean2, emean2, resp2, eresp2, _, _ = HandleHistograms.fit(hist2, pars2, histo_ranges2)

    #shift the hits to calculate resolution
    bokehplot.add_frame('pure_rechit_energy_shifted.html', nfigs=len(true_beam_energies_GeV))
    pars1_shifted = tuple([x[0], x[1]*average_shift, x[2]] for x in pars1)
    data1_shifted = ProcessData.shift_energy(data1, mean1)
    hist1_shifted = HandleHistograms.create(data1_shifted, bins)
    _, _, _, _, res1_shifted, eres1_shifted = HandleHistograms.fit(hist1_shifted, pars1_shifted, histo_ranges1_shifted)

    bokehplot.add_frame('clusterized_rechit_energy_shifted.html', nfigs=len(true_beam_energies_GeV))
    pars2_shifted = tuple([x[0], x[1]*average_shift, x[2]] for x in pars2)
    data2_shifted = ProcessData.shift_energy(data2, mean2)
    hist2_shifted = HandleHistograms.create(data2_shifted, bins)
    _, _, _, _, res2_shifted, eres2_shifted = HandleHistograms.fit(hist2_shifted, pars2_shifted, histo_ranges2_shifted)

    bokehplot.add_frame('response_and_resolution.html', nfigs=3)
    last_frame_id = bokehplot.get_nframes() - 1
    final_graphs(resp1, eresp1, res1_shifted, eres1_shifted, resp2, eresp2, res2_shifted, eres2_shifted, last_frame_id)

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    assert(len(beam_energies)==len(true_beam_energies_GeV))

    bokehplot = bkp.BokehPlot(filenames='pure_rechit_energy.html', nfigs=len(true_beam_energies_GeV))
    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
    main()


"""
Errors for 'resp1 / resp2':
[np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 ), 
np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 )]],
"""
