# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/09/17
content:    Script for figure 2.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset
from singlecell.utils.distplot import distplot


# Classes / functions
class Fig1(Dataset):
    pass


# Script
if __name__ == '__main__':

    ds = Fig1(
            counts_table='dengue',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)

    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)
    ds.counts.normalize('counts_per_million', inplace=True)

    n = ds.samplesheet['numberDengueReads'].astype(int)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)
    ds.counts.log(inplace=True)
    dst = ds.split(phenotypes=('time', 'MOI'))
    gene = 'DDIT3'
    gid = ds.featuresheet.index[ds.featuresheet['GeneName'] == gene][0]
    colors = sns.color_palette(n_colors=3)
    fig, axs = plt.subplots(
         1, 2,
         sharey=True,
         figsize=(6.4, 3.1))
    for it, time in enumerate(('4', '12', '24', '48')):
        for im, moi in enumerate(('0', '1', '10')):
            dsti = dst[(time, moi)]

            x = [
                np.maximum(-0.5, np.log10(ds.counts.pseudocount + dsti.samplesheet['virus_reads_per_million'].values)),
                np.maximum(-0.5, dsti.counts.loc[gid].values),
                ]

            # Singular value for totally uninfected cells
            if x[0].max() == -0.5:
                x[0] = x[0] + 0.1 * np.random.rand(x[1].shape[0])

            for ico, (ax, co) in enumerate(zip(axs, x)):
                if False:#it == 0 and ico == 1:
                    label = moi
                else:
                    label = None

                distplot(
                        co,
                        hist=False,
                        rug=False,
                        kde=True,
                        ax=ax,
                        color=colors[im],
                        label=label,
                        bottom=4 - it)

            ## Plot an arrow at the arithmetic mean
            #vm = np.maximum(-0.5, np.log10(ds.counts.pseudocount + dsti.samplesheet['virus_reads_per_million'].values.mean()))
            #axs[0].arrow(
            #        vm, 4.9 - it, 0, -0.9,
            #        lw=3,
            #        head_width=0.1,
            #        head_length=0.1,
            #        length_includes_head=True,
            #        color=colors[im])

            ## Plot another arrow at the median
            #vg = np.maximum(-0.5, np.median(np.log10(ds.counts.pseudocount + dsti.samplesheet['virus_reads_per_million'].values)))
            #axs[0].arrow(
            #        vg, 4.9 - it, 0, -0.9,
            #        lw=3,
            #        ls='--',
            #        head_width=0.1,
            #        head_length=0.1,
            #        length_includes_head=True,
            #        color=colors[im])

    for iax, ax in enumerate(axs):
        ax.set_yticks((1, 2, 3, 4))
        ax.set_yticklabels(('4', '12', '24', '48')[::-1])
        ax.set_xticks((-0.5, 0, 1, 2, 3, 4, 5))
        ax.set_xticklabels(('$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'))
        ax.set_ylim(0.5, 6.5)
        if iax == 1:
            ax.set_xlim(-0.5, 4)
            ax.set_xlabel(gene+' counts per million transcripts')
        else:
            ax.set_xlim(-0.5, 5.5)
            ax.set_xlabel('Virus per million transcripts')
            ax.set_ylabel('Time since infection [hrs]')
            ax.yaxis.set_label_coords(-0.15, 0.35)
        ax.grid(True)

    #axs[1].legend(loc='upper right', title='MOI')
    axs[0].annotate(
            '',
            xy=(-0.12, 0),
            xycoords='axes fraction',
            xytext=(-0.12, 0.7),
            arrowprops=dict(arrowstyle="->", color='k'))

    plt.tight_layout(rect=(0.01, 0, 1, 1), w_pad=0)

#    print('Plot number of infected cells')
#    dq = Dataset(
#            counts_table=None,
#            samplesheet='quantified',
#            featuresheet=None,
#            )
#    dq.samplesheet.rename(columns={
#        'ratio [pmol/mg]': 'ratio',
#        'qPCR_ACTB [pM]': 'ACTB',
#        'time [h]': 'time'},
#        inplace=True)
#
#    # Count cells for each condition (table 1)
#
#    dq.query_samples_by_metadata('experiment == "10017006"', inplace=True)
#    dq.samplesheet.loc[dq.samplesheet['ratio'] == '', 'ratio'] = 0
#    dq.samplesheet.loc[:, 'ratio'] = dq.samplesheet.loc[:, 'ratio'].astype(float)
#    dq.samplesheet.loc[dq.samplesheet['ACTB'] == '', 'ACTB'] = 1
#    dq.samplesheet.loc[:, 'ACTB'] = dq.samplesheet.loc[:, 'ACTB'].astype(float)
#    dq.query_samples_by_metadata('ACTB >= 0.01', inplace=True)
#    dq.samplesheet['infected'] = dq.samplesheet['ratio'] > 1e-6
#
#    n = (dq.samplesheet[['time', 'MOI', 'infected']]
#           .groupby(['time', 'MOI'])
#           .mean()
#           .loc[:, 'infected']
#           .unstack('MOI'))
#    n.loc['0'] = 0
#    n = n.loc[['0', '4', '12', '24', '48']]
#    fig, ax = plt.subplots(1, 1, figsize=(3.6, 3.1))
#    x = [0, 4, 12, 24, 48]
#    for moi in ('0', '1', '10'):
#        ax.plot(
#            x, n.loc[:, moi],
#            '-o',
#            lw=2,
#            label=moi)
#
#    ax.set_xlabel('Time since infection [hrs]')
#    ax.set_ylabel('Fraction of cells infected')
#    ax.legend(loc='upper left', title='MOI')
#    ax.grid(True)
#
#    plt.tight_layout()

    plt.ion()
    plt.show()
