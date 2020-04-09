import sys
import os
import subprocess
import glob
import pandas as pd
import numpy as np
from functools import reduce
import bokehplot as bkp

cmssw_base = subprocess.check_output("echo $CMSSW_BASE", shell=True).split('\n')[0]
usercode_path = 'src/UserCode/DataProcessing/job_output_old'
path = os.path.join(cmssw_base, usercode_path, 'out_')
dfs = []
for j,i in enumerate(glob.glob(path + '*.csv')):
    df_tmp = pd.read_csv(i)
    dfs.append(df_tmp)
#merge all columns and set values to zero whenever some DataFrames have more rows than others
df = reduce(lambda left,right: pd.merge(left, right, left_index=True, right_index=True, how='outer').fillna(-1.), dfs)

beam_energies = (20,30,50,80,100,120,150,200,250,300)
df_beamen = []
for b in beam_energies:
    #energy sum columns that correspond to a specific beam energy
    #columns with the same energy are put into the same histogram
    cols = [x for x in df.columns if ( 'ensum' in x and df.at[1,'beamen'+x[-3:]]==b ) ]
    df_beamen.append( df.loc[:, cols].stack() if cols != [] else None )

hist = []
for idf in df_beamen:
    if idf is not None:
        idf = idf[idf>-1]
        range_min = -1. if idf.min() <= -1. else idf.min()
        hist.append( np.histogram(idf, density=False, bins=100, range=(range_min,idf.max())) ) #flattens 'df' to one-dimension

b = bkp.BokehPlot(filenames='200GeV_test.html', nfigs=len(beam_energies))
line_colors = ['black', 'blue', 'green', 'red', 'orange', 'purple', 'greenyellow', 'brown', 'pink', 'grey']
axis_kwargs = {'x.axis_label': 'Total RecHit energy per event [MeV]', 'y.axis_label': 'Counts'}
b.histogram(data=hist, idx=[x for x in range(len(beam_energies))], 
            legend_label=[str(x)+' GeV' for x in beam_energies], fill_color='white', line_color=line_colors,
            fig_kwargs=axis_kwargs)
common_args = {'pdf':'gaus', 'line_width':2, 'alpha':0.7, 'legend_label':'gaussian'}
p0_parameters = ([10000, 1600.,  300.], #20GeV
                 [15000, 2300.,  300.], #30GeV
                 [30000, 4000.,  300.], #50GeV
                 [30000, 6500.,  300.], #80GeV
                 [30000, 8000.,  300.], #100GeV
                 [25000, 9500.,  300.], #120GeV
                 [40000, 12000., 300.], #150GeV
                 [40000, 16000., 300.], #200GeV
                 [35000, 19500., 300.], #250GeV
                 [2500,  23500., 300.]) #300GeV

resolutions = []
responses = []
for i in range(len(beam_energies)):
    coeff, var = b.fit(p0=p0_parameters[i], idx=i, color=line_colors[i], **common_args)
    err = np.sqrt(np.diag(var))
    mean_label = 'mean='+str(round(coeff[1],1))+'+-'+str(round(err[1],1))+' MeV'
    sigma_label = 'sigma='+str(round(coeff[2],1))+'+-'+str(round(err[2],1))+' MeV'
    font_size = {'text_font_size': '8pt'}
    b.label(mean_label,  idx=i, x=10, y=320, **font_size)
    b.label(sigma_label, idx=i, x=10, y=305, **font_size)
    responses.append( (coeff[1]-beam_energies[i])/beam_energies[i] )
    resolutions.append( coeff[2] )
b.show_frame(plot_width=300, plot_height=300)




import numpy as np
import matplotlib.pyplot as plt
plt.subplot(2, 1, 1)
plt.plot([x for x in range(len(beam_energies))], responses, 'o-')
plt.title('Plots')
plt.ylabel('Response')
plt.subplot(2, 1, 2)
plt.plot([x for x in range(len(beam_energies))], resolutions, '.-')
plt.ylabel('Resolution')
plt.savefig('test.png')
#b.add_figure()
