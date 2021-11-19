import os
import numpy as np
import subprocess
import bokehplot as bkp
import pandas as pd
import argparse
from argparser import add_args

def main():
    eos_base = '/eos/user/'
    cms_user = subprocess.check_output("echo $USER", shell=True, encoding='utf-8').split('\n')[0]
    data_directory = 'TestBeamReconstruction'
    output_html_dir = os.path.join(eos_base, cms_user[0], cms_user, 'www', data_directory, 'cluster_dep', FLAGS.datatype, FLAGS.showertype)
    nlayers = 28 if FLAGS.showertype == 'em' else 40

    outlambda = lambda x: os.path.join(output_html_dir, FLAGS.datatype + '_' + FLAGS.showertype + '_' + str(FLAGS.chosen_energy) + x)
    if FLAGS.var == 'dx':
        out = outlambda('_clusters_dx_summary.html')
        label = 'x'
    elif FLAGS.var == 'dy':
        out = outlambda('_clusters_dy_summary.html')
        label = 'y'

    bokehplot = bkp.BokehPlot(filenames=out, nfigs=2, nframes=1)
    plot_width, plot_height = 600, 400
    true_energies = (20,30,49.99,79.93,99.83,119.65,149.14,197.32,243.61,287.18)
    true_en = true_energies[min( range(len(true_energies)), key=lambda i: abs(true_energies[i]-FLAGS.chosen_energy) )]
    figkwargs = {'t.text': '{} GeV beam energy'.format(true_en), 'x.axis_label': 'Layer',
                 'plot_width': plot_width, 'plot_height': plot_height,
                 'l.location': 'top_left', 'l.click_policy': 'hide', 'l.label_text_font_size': '7pt'}
    compression_level = 9

    hdf5_path_start = os.path.join(eos_base, cms_user[0], cms_user, data_directory)

    tags = ['W2p9_dpos1p3', 'W2p3_dpos1p3', 'W4p0_dpos1p3', 'W5p0_dpos1p3', 'W2p9_dpos3p4', 'W2p9_dpos999', 'W5p0_dpos999']
    colors = ['green', 'orange', 'blue', 'purple', 'turquoise', 'red', 'grey']
    assert(len(colors) == len(tags))
    for i,t in enumerate(tags):
        n = os.path.join( hdf5_path_start, t, 'cluster_dependent',
                          'cluster_' + FLAGS.showertype + '_' + FLAGS.datatype + '_' + str(FLAGS.chosen_energy) + 'GeV_' +
                          t + '_' + FLAGS.var + '_summary.h5' )

        bias = pd.read_hdf(n,  'bias' + label)
        ebias = pd.read_hdf(n, 'ebias'+ label)
        res = pd.read_hdf(n,   'res'  + label)
        eres = pd.read_hdf(n,  'eres' + label)

        if i==0:
            bkw = dict(figkwargs, **{'y.axis_label': 'Position ' + label.upper() + ' bias [cm]'})
            rsw = dict(figkwargs, **{'y.axis_label': 'Position ' + label.upper() + ' resolution [cm]'})
        else:
            bkw, rsw = None, None

        #tupdate = 'W0='+t[2:5].replace('p', '.')+', d='+t[-3:].replace('p', '.')
        tupdate = 'W='+t[1:4].replace('p', '.')+', d='+t[-3:].replace('p', '.')
        bokehplot.graph(data=[np.arange(1,nlayers+1), np.array(bias)],
                        errors=[[np.zeros(nlayers),np.zeros(nlayers)],[np.array(ebias)/2, np.array(ebias)/2]],
                        idx=0, iframe=0, line=True, color=colors[i], legend_label=tupdate,
                        fig_kwargs=bkw)
        bokehplot.graph(data=[np.arange(1,nlayers+1), np.array(res)],
                        errors=[[np.zeros(nlayers),np.zeros(nlayers)],[np.array(eres)/2, np.array(eres)/2]],
                        idx=1, iframe=0, line=True, color=colors[i], legend_label=tupdate,
                        fig_kwargs=rsw)

    bokehplot.save_frame(plot_width=plot_width, plot_height=plot_height, show=False)

if __name__ == '__main__':
    #define parser for user input arguments
    parser = argparse.ArgumentParser()
    FLAGS, _ = add_args(parser, 'summary')
    main()
