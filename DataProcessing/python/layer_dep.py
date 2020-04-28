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
            cols = [x for x in df_reduced.columns if ( ('nhitsfrac' in x or 'enfrac' in x) and df_reduced.at[1,'beamen'+x[-3:]]==b ) ]
            df_beamen.append( df_reduced.loc[:, cols].stack() if cols != [] else None )
        return df_beamen
            
def histograms(dfs):
    """Plots responses and resolutions with their errors"""
    axis_kwargs1 = {'x.axis_label': 'Fraction', 'y.axis_label': 'Counts'}
    hist1, hist2 = ([] for _ in range(2))
    bins = np.linspace(0,1,11)
    for i,idf in enumerate(dfs):
        nhitsfrac = idf['nhitsfrac']
        enfrac    = idf['endfrac']
        #flattens 'idf' to one-dimension
        hist1.append( np.histogram(nhitsfrac, density=False, bins=bins, range=(nhitsfrac.min(),nhitsfrac.max())) )

    bokehplot.histogram(data=hist1, idx=[x for x in range(len(true_beam_energies_GeV))], style='step',
                        legend_label=['fraction of clusterized hits'], line_color='blue',
                        fig_kwargs=axis_kwargs)
    bokehplot.histogram(data=hist2, idx=[x for x in range(len(true_beam_energies_GeV))], style='step',
                        legend_label=['fraction of clusterized energy'], line_color='red',
                        fig_kwargs=axis_kwargs)
    font_size = {'text_font_size': '13pt', 'text_color': 'green'}
    bokehplot.label([str(x)+' GeV' for x in true_beam_energies_GeV], 
                    idx=[x for x in range(len(true_beam_energies_GeV))], 
                    x=10, y=320, **font_size)

    bokehplot.show_frame(plot_width=300, plot_height=300)

def main():
    usercode_path = 'src/UserCode/DataProcessing/job_output/layer_dependent/'
    path = os.path.join(cmssw_base, usercode_path, 'out_')

    data1 = ProcessData.join(path + '*layerdep.csv')
    histograms(data1)

if __name__ == '__main__':
    cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
    beam_energies = (20,30,50,80,100,120,150,200,250,300)
    true_beam_energies_GeV = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_beam_energies_MeV = tuple(x*1000 for x in true_beam_energies_GeV)
    assert(len(beam_energies)==len(true_beam_energies_GeV))

    bokehplot = bkp.BokehPlot(filenames='plot_layers.html', nfigs=len(true_beam_energies_GeV))
    line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
    main()
