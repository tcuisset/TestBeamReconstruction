import sys
import os
import subprocess
import glob
import pandas as pd
import numpy as np
from functools import reduce
import bokehplot as bkp
from bokeh.models import Range1d

def create_histograms(path, clusters_bool):
    dfs = []
    end_string = '*[0-9][0-9][0-9].csv' if clusters_bool else '*_noclusters.csv'
    for j,i in enumerate(glob.glob(path + end_string)):
        df_tmp = pd.read_csv(i)
        dfs.append(df_tmp)
    #merge all columns and set values to zero whenever some DataFrames have more rows than others
    df = reduce(lambda left,right: pd.merge(left, right, left_index=True, right_index=True, how='outer').fillna(-1.), dfs)

    df_beamen = []
    for b in beam_energies:
        #energy sum columns that correspond to a specific beam energy
        #columns with the same energy are put into the same histogram
        cols = [x for x in df.columns if ( 'ensum' in x and df.at[1,'beamen'+x[-3:]]==b ) ]
        df_beamen.append( df.loc[:, cols].stack() if cols != [] else None )

    bins = (1000, 1800, 4200, 5000, 5000, 4200, 5700, 5500, 5500, 500)
    hist = []
    for i,idf in enumerate(df_beamen):
        if idf is not None:
            idf = idf[idf>-1]
            range_min = -1. if idf.min() <= -1. else idf.min()
            hist.append( np.histogram(idf, density=False, bins=bins[i], range=(range_min,idf.max())) ) #flattens 'df' to one-dimension

    axis_kwargs = {'x.axis_label': 'Total RecHit energy per event [MeV]', 'y.axis_label': 'Counts'}
    bokehplot.histogram(data=hist, idx=[x for x in range(len(true_beam_energies))], style='step',
                        legend_label=[str(x)+' GeV' for x in true_beam_energies], fill_color='white', line_color=line_colors, alpha=0.5,
                        fig_kwargs=axis_kwargs)
    return hist

def fits(hist, clusters_bool):
    sigma_units = 3
    histo_ranges = (Range1d(0,30000), Range1d(11000, 35000), Range1d(27000, 58000), Range1d(52000, 94000), 
                    Range1d(64000, 120000), Range1d(88000,130000), Range1d(120000,165000), Range1d(160000, 210000), 
                    Range1d(210000,250000), Range1d(240000,295000))
    common_args = {'pdf':'gaus', 'line_width':2.5, 'alpha':0.8}
    p0_parameters = ( ([7500,  20000.,  2000.], #20GeV
                       [14000,  30000.,  1200.], #30GeV
                       [35000, 50000.,  2100.], #50GeV
                       [40000, 80000.,  2500.], #80GeV
                       [30000, 100000.,  3000.], #100GeV
                       [25000, 120000.,  2800.], #120GeV
                       [40000, 150000., 3000.], #150GeV
                       [40000, 200000., 3500.], #200GeV
                       [32000, 215000., 4000.], #250GeV
                       [2500,  280000., 4500.]) #300GeV 
                  ) if not clusters_bool else (
                      ([7000,  18000.,  2000.], #20GeV
                       [14000, 25000.,  1200.], #30GeV
                       [35000, 43000.,  2000.], #50GeV
                       [40000, 77000.,  2500.], #80GeV
                       [30000, 98000., 2700.], #100GeV
                       [25000, 109000., 2750.], #120GeV
                       [40000, 140000., 3500.], #150GeV
                       [35000, 180000., 3500.], #200GeV
                       [32000, 220000., 4000.], #250GeV
                       [2300,  270000., 4500.]) #300GeV 
                  )

    responses = []
    responses_err = []
    resolutions = []
    resolutions_err = []
    for i in range(len(true_beam_energies)):
        #First fit
        coeff, _ = bokehplot.fit(p0=p0_parameters[i], idx=i, obj_idx=0, color=line_colors[i], **common_args)
        fit_bounds = (coeff[1]-sigma_units*coeff[2], coeff[1]+sigma_units*coeff[2])

        #Create and draw amputed histogram 
        #Detail: the last item of the histogram with the counts is removed since it refers to the events between the last
        #edge selected and the next one. Nedges = Nbinswithcounts + 1
        selection = (hist[i][1]>fit_bounds[0]) & (hist[i][1]<fit_bounds[1])
        indexes_selected = np.nonzero(selection) #return indexes of the non-zero (True) elements
        hist_modified = [hist[i][0][indexes_selected][:-1], hist[i][1][selection]]
        bokehplot.histogram(data=hist_modified, idx=i, style='step',
                            legend_label=[str(x)+' GeV' for x in true_beam_energies], fill_color='white', line_color=line_colors[i], 
                            alpha=1., fig_kwargs={'x_range': histo_ranges[i]})

        #Second fit
        coeff, var = bokehplot.fit(p0=p0_parameters[i], idx=i, obj_idx=1, color=line_colors[i], **common_args)
        err = np.sqrt(np.diag(var))
        err1 = round(err[1],2)
        err2 = round(err[2],2)
        mean_label = 'mean='+str(round(coeff[1],2))+'+-'+str(err1)+' MeV'
        sigma_label = 'sigma='+str(round(coeff[2],2))+'+-'+str(err2)+' MeV'
        font_size = {'text_font_size': '8pt'}
        bokehplot.label(mean_label,  idx=i, x=10, y=320, **font_size)
        bokehplot.label(sigma_label, idx=i, x=10, y=305, **font_size)
        responses.append( (coeff[1]-true_beam_energies_MeV[i]) / (true_beam_energies_MeV[i]) ) 
        resolutions.append( coeff[2] / true_beam_energies_MeV[i] )
        responses_err.append( err1 / true_beam_energies_MeV[i] )
        resolutions_err.append( err2 / true_beam_energies_MeV[i] )
    bokehplot.show_frame(plot_width=300, plot_height=300)
    return responses, responses_err, resolutions, resolutions_err

def final_graphs(resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2):
    """Plots responses and resolutions with their errors"""
    axis_kwargs1 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response (E/True - 1)'}
    bokehplot.graph(data=[np.array(true_beam_energies),np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies)),np.zeros(len(true_beam_energies))],
                            [np.array(eresp1)/2,np.array(eresp1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs1)
    bokehplot.graph(data=[np.array(true_beam_energies),np.array(resp2)], 
                    errors=[[np.zeros(len(true_beam_energies)),np.zeros(len(true_beam_energies))],
                            [np.array(eresp2)/2,np.array(eresp2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')

    axis_kwargs2 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Response difference'}
    bokehplot.graph(idx=1, data=[np.array(true_beam_energies), np.array(resp2)-np.array(resp1)], 
                    errors=[[np.zeros(len(true_beam_energies)),np.zeros(len(true_beam_energies))],
                            [np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 ), 
                             np.sqrt( ( np.power(np.array(eresp1),2)+np.power(np.array(eresp2),2) ) / 2 )]],
                    style='square', line=True, color='blue', legend_label="Reconstructable - Clusterized", fig_kwargs=axis_kwargs2)
    fig = bokehplot.get_figure(idx=1, iframe=2)
    fig.legend.location = 'bottom_right'

    axis_kwargs3 = {'x.axis_label': 'Beam energy [GeV]', 'y.axis_label': 'Resolution [MeV]'}
    bokehplot.graph(idx=2, data=[np.array(true_beam_energies),np.array(res1)], 
                    errors=[[np.zeros(len(true_beam_energies)),np.zeros(len(true_beam_energies))],
                            [np.array(eres1)/2,np.array(eres1)/2]],
                    style='square', line=True, color='green', legend_label='Reconstructable', fig_kwargs=axis_kwargs3)
    bokehplot.graph(idx=2, data=[np.array(true_beam_energies),np.array(res2)], 
                    errors=[[np.zeros(len(true_beam_energies)),np.zeros(len(true_beam_energies))],
                            [np.array(eres2)/2,np.array(eres2)/2]],
                    style='triangle', line=True, color='orange', legend_label='Clusterized hits')
    bokehplot.show_frame(plot_width=300, plot_height=300)

def main():
    #files with sum of rechit energy
    usercode_path = 'src/UserCode/DataProcessing/job_output'
    path = os.path.join(cmssw_base, usercode_path, 'out_')
    clusters_bool = False
    histograms1 = create_histograms(path, clusters_bool)
    resp1, eresp1, res1, eres1 = fits(histograms1, clusters_bool)
    bokehplot.add_frame('clusterized_rechit_energy.html', nfigs=len(true_beam_energies))

    #files with sum of clusterized rechit energy
    clusters_bool = True
    histograms2 = create_histograms(path, clusters_bool)
    resp2, eresp2, res2, eres2 = fits(histograms2, clusters_bool)
    bokehplot.add_frame('response_and_resolution.html', nfigs=3)
    final_graphs(resp1, eresp1, res1, eres1, resp2, eresp2, res2, eres2)

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies)
    assert(len(beam_energies)==len(true_beam_energies))

    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
    bokehplot = bkp.BokehPlot(filenames='pure_rechit_energy.html', nfigs=len(true_beam_energies))
    main()


"""
Erros for 'resp1 / resp2':
[np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 ), 
np.sqrt( ( np.power(np.array(eresp2),2)/np.power(np.array(resp1),2) + np.power(np.array(eresp1)*np.array(resp2),2)/np.power(resp1,4) ) / 2 )]],
"""
