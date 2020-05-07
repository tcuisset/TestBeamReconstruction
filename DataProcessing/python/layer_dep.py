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

class ProcessData:
    @staticmethod
    def join(glob_path, field):
        dfs = [None for _ in range(len(beam_energies))]
        col_indexes = [0]
        col_indexes.extend([x for x in range(1,2*nlayers+1,2)]) if 'nhits' in field else col_indexes.extend([x for x in range(0, 2*nlayers+1, 2)])
        print('grouping dataframes...')
        for j,i in enumerate(glob.glob(glob_path)):
            print('{}/111 dataframes... \r'.format(j+1))
            df_tmp = pd.read_csv(i, usecols=col_indexes)
            tuple_index = beam_energies.index( df_tmp['beamen'].iloc[0] )
            df_tmp = df_tmp.iloc[:,1:] #remove beam energy column; we do not need it anymore!

            #rename columns for pd.concat() to work (remove ntuple code)
            df_tmp.columns = [x[:-4] for x in df_tmp.columns] 

            if dfs[tuple_index] is None:
                dfs[tuple_index] = df_tmp
            else:
                dfs[tuple_index] = pd.concat([dfs[tuple_index], df_tmp], ignore_index=True, sort=True)
        return dfs
            
def histograms(dfs, field):
    """Plots responses and resolutions with their errors"""
    xlabelextra = ' of clusterized hits' if field=='nhitsfrac_' else ' of clusterized energy'
    axis_kwargs1 = {'y.axis_label': 'Fraction'+xlabelextra, 'x.axis_label': 'Layer'}
    nbins = 30
    bins = np.linspace(0,1,nbins+1)
    height_hits = height_for_plot_bins(bins, scale='linear')
    for i,idf in enumerate(dfs):
        #bokehplot.add_frame(os.path.join(output_html_dir, 'plot_layers_'+field+str(i)+'.html'), nfigs=1)
        xvalues, yvalues, counts = ([] for _ in range(3))
        for ilayer in range(nlayers):
            frac = idf[field+'layer'+str(ilayer)]
            hist = np.histogram(frac, density=False, bins=bins, range=(frac.min(),frac.max()))
            ncounts = len(hist[0])
            thislayer = ilayer + 1
            bincenters = (hist[1][:-1]+hist[1][1:])/2
            xvalues.extend( thislayer*np.ones(len(bincenters)) ) #layer values
            yvalues.extend( bincenters ) #fraction values
            counts.extend( hist[0] ) #count values

        fig_kwargs = {'plot_width': plot_width, 'plot_height': plot_height,
                      't.text': 'True beam energy: {} GeV'.format(true_beam_energies_GeV[i])}
        fig_kwargs.update(axis_kwargs1)
        print('{}/{} 2d graphs...'.format(i+1,size))
        bokehplot.graph(data=[np.array(xvalues), np.array(yvalues), np.array(counts)],
                        width=np.ones((len(xvalues))), height=height_hits,
                        idx=i, style='rect%Viridis', fig_kwargs=fig_kwargs)
        #bokehplot.save_frame()

def create_dir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

def main():
    usercode_path = 'src/UserCode/DataProcessing/job_output/layer_dependent/'
    path = os.path.join(cmssw_base, usercode_path, 'outEcut_*_layerdep.csv')

    data_hits = ProcessData.join(path, 'nhitsfrac_')
    histograms(data_hits, 'nhitsfrac_')
    bokehplot.save_frame(show=False)

    bokehplot.add_frame(os.path.join(output_html_dir, 'plot_layers_energies.html'), nfigs=size)
    data_en = ProcessData.join(path, 'enfrac_')
    histograms(data_en, 'enfrac_')
    bokehplot.save_frame(show=False)

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    cms_user = subprocess.check_output("echo $USER", shell=True).split('\n')[0]
    nlayers = 28

    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    size = len(true_beam_energies_GeV)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    assert(len(beam_energies)==len(true_beam_energies_GeV))

    data_directory = 'TestBeamReconstruction'
    create_dir( os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory) )
    output_html_dir = os.path.join('/eos/user/', cms_user[0], cms_user, 'www', data_directory)
    output_html_file = os.path.join(output_html_dir, 'plot_layers_hits.html')
    bokehplot = bkp.BokehPlot(filenames=output_html_file, nfigs=size)
    plot_width, plot_height = 600, 300
    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
    main()
