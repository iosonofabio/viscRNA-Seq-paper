# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/09/17
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


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('--virus', choices=['dengue', 'Zika'], default='dengue',
                        help='Virus to look at')
    args = parser.parse_args()
    virus = args.virus

    print('Load dataset')
    ds = Dataset(
            counts_table=virus.lower(),
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)

    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

    print('Normalize')
    ds.counts.normalize('counts_per_million', inplace=True)

    print('Add virus reads')
    n = ds.samplesheet['number{:}Reads'.format(virus.capitalize())].astype(int)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)

    print('Log counts')
    ds.counts.log(inplace=True)

    print('Get correlations')
    co = ds.correlation.correlate_features_phenotypes(
            phenotypes='virus_reads_per_million',
            fillna=0).fillna(0)

    # Plot ER-stress/proapoptotic genes
    gnames = ['BBC3', 'TRIB3', 'TNFRSF10B', 'EIF2A',
              'DDIT3', 'PPP1R15A', 'ATF4', 'EDEM1',
              'XBP1', 'ATF6', 'ATF3', 'ATG5',
              'CASP3', 'CASP9', 'CASP4', 'CASP6',
              ]
    dsd = ds.copy()
    gids = [ds.featuresheet.index[ds.featuresheet['GeneName'] == gname][0] for gname in gnames]
    dsd.counts = dsd.counts.loc[gids]

    dsh = dsd.copy()
    if virus == 'dengue':
        dsh.query_samples_by_metadata('virus_reads_per_million > 1e3', inplace=True)
        dsh.query_samples_by_metadata('MOI != "0"', inplace=True)
    else:
        dsh.query_samples_by_metadata('virus_reads_per_million > 1e1', inplace=True)
        dsh.query_samples_by_metadata('MOI == "1"', inplace=True)
    fig, axs = plt.subplots(
            sharex=True, sharey=True,
            nrows=4, ncols=4, figsize=(12, 10))
    axs = axs.ravel()
    coh = dsh.correlation.correlate_features_phenotypes(['virus_reads_per_million'], fillna=0)
    gids_sorted = sorted(gids, key=lambda x: coh.loc[x, 'virus_reads_per_million'], reverse=True)
    x = dsh.samplesheet.loc[:, 'virus_reads_per_million'].values
    for gid, ax in zip(gids_sorted, axs):
        gname = dsh.featuresheet.loc[gid, 'GeneName']
        y = dsh.counts.loc[gid].values
        ax.scatter(
                x, y,
                alpha=0.5,
                s=10,
                label='$\\rho = {:.2f}$'.format(coh.loc[gid, 'virus_reads_per_million']),
                zorder=10)
        ax.set_title(gname)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.grid(True)
        ax.set_ylim(-1.1, 5.5)
        ax.set_yticks([-1, 0, 2, 4])
        ax.set_yticklabels(['$0$', '$1$', '$10^2$', '$10^4$'])
        if virus == 'dengue':
            ax.set_xlim(3e3, 10**(5.5))
        else:
            ax.set_xlim(3e1, 10**(4.5))
        ax.set_xscale('log')
        ax.legend(loc='best')
    fig.text(0.5, 0.02,
             virus.capitalize()+' amount per million transcripts',
             ha='center')
    fig.text(0.02, 0.5,
             'Gene counts per million transcripts',
             ha='right', va='center', rotation=90)
    plt.tight_layout(rect=(0.02, 0.03, 1, 1), w_pad=1)

    plt.ion()
    plt.show()
