# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/09/17
content:    Script for studying bystander effects.
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

    ds = Dataset(
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

    # Only select cells without virus
    ds.query_samples_by_metadata('virus_reads_per_million < 0.1', inplace=True)

    # Check table with number of cells
    table = (ds.samplesheet
               .groupby(['time', 'MOI'])
               .count()
               .iloc[:, 0]
               .unstack()
               .fillna(0)
               .astype(int)
               .loc[['4', '12', '24', '48']])

    print('Selecting only early 2 time points')
    # The rest has too few uninfected cells
    ds.query_samples_by_metadata('time in ["4", "12"]', inplace=True)
    dsm = ds.split('MOI')
    ks = dsm['0'].compare(dsm['1'])['P-value']

    # Get the top hits for GO analysis
    hits = ds.featuresheet.loc[ks.nsmallest(100).index, 'GeneName'].values
    with open('../tables/bystander_top_100.tsv', 'w') as f:
        f.write('\n'.join(hits))

    sys.exit()

    # Bonferroni correction
    ks = np.minimum(1, ks * len(ks))

    # Print cumulative histogram of P-values
    x = ks.sort_values().values
    y = 1.0 - np.linspace(0, 1, len(x))
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(x, y, lw=2, color='darkred')
    ax.set_xlabel('KS P-value')
    ax.set_ylabel('Cumulative distribution')
    ax.set_xscale('log')

    plt.tight_layout()

    # Plot the top hits together with a control gene
    gcontrols = ['HM13', 'SQSTM1', 'PSMB2']
    gidcontrols = [ds.featuresheet.index[(ds.featuresheet['GeneName'] == gc)][0] for gc in gcontrols]
    gids = list(ks.nsmallest(3).index) + gidcontrols
    gnames = ds.featuresheet.loc[gids, 'GeneName'].values
    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 6), sharex=True, sharey=True)
    axs = axs.ravel()
    for i, (ax, gid, gname) in enumerate(zip(axs, gids, gnames)):
        x0 = np.sort(10**(dsm['0'].counts.loc[gid].values))
        y0 = 1.0 - np.linspace(0, 1, len(x0))
        x1 = np.sort(10**(dsm['1'].counts.loc[gid].values))
        y1 = 1.0 - np.linspace(0, 1, len(x1))
        ax.plot(x0, y0, lw=2, color='steelblue', label='control')
        ax.plot(x1, y1, lw=2, color='darkred', label='bystander')
        ax.set_title(gname)
        ax.grid(True)
        ax.set_xscale('log')
        pval = ks.loc[gid]
        ax.text(0.05, 0.06, 'P = {:.2f}'.format(pval),
                fontsize=9,
                ha='left', va='bottom', transform=ax.transAxes)

        if i == 0:
            ax.legend(loc='upper left', title='Cell condition')

        if i == 4:
            ax.set_xlabel('Counts per million transcripts')

    fig.text(0.02, 0.5, 'Cumulative distribution', rotation=90, va='center')

    plt.tight_layout(rect=(0.04, 0, 1, 1), w_pad=0.5)

    plt.ion()
    plt.show()
