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


def load_cell_cycle_table():
    fn = '/home/fabio/university/postdoc/virus_singlecell/data/cell_cycle/cr201684x7.tsv'
    tb = pd.read_csv(fn, sep='\t', index_col=0)
    tb.index.name = 'EnsemblID'
    tb.rename(columns={
        'Symbol': 'GeneName',
        'Phase ': 'Phase',
        'Core 67': 'Core',
        }, inplace=True)
    return tb


# Script
if __name__ == '__main__':

    print('Load dataset')
    ds = Dataset(
            counts_table='dengue',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)
    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    ds.samplesheet.loc[:, 'time'] = pd.Categorical(
            ds.samplesheet.loc[:, 'time'].astype(int),
            )
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

    print('Normalize')
    ds.counts.normalize('counts_per_million', inplace=True)

    print('Add virus reads')
    n = ds.samplesheet['numberDengueReads'].astype(int)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)
    ds.samplesheet['log_virus_reads_per_million'] = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])

    #print('Log counts')
    #ds.counts.log(inplace=True)

    print('Get cell cycle genes (Core 67)')
    cc = load_cell_cycle_table()[['GeneName', 'Periodic Rank', 'Phase', 'Core']]
    cc.query('Core != "No"', inplace=True)

    print('Hierarchical clustering of cells based on cell cycle and virus')
    dsn = ds.copy()
    dsn.counts = dsn.counts.loc[cc.index]
    for col in ['Periodic Rank', 'Phase']:
        dsn.featuresheet.loc[:, col] = cc.loc[
                dsn.featuresheet.index, col]
    dsn.rename(axis='features', column='GeneName', inplace=True)

    print('Log counts')
    dsn.counts.log(inplace=True)
    dsn.featuresheet.loc[:, 'Mean'] = dsn.counts.mean(axis=1)

    # Only keep decently expressed genes
    dsn.query_features_by_metadata('Mean > 1', inplace=True)

    hier = dsn.cluster.hierarchical(
            axis='samples',
            #phenotypes=['virus_reads_per_million'],
            optimal_ordering=True,
            log_features=False,
            )
    #dsn.counts = dsn.counts.loc[:, hier['leaves']]

    #hier_fea = dsn.cluster.hierarchical(
    #        axis='features',
    #        optimal_ordering=True,
    #        log_features=False,
    #        )

    dsn.plot.clustermap(
            cluster_samples=hier['linkage'],
            cluster_features=False,
            labels_samples=False,
            annotate_features={
                'Phase': 'Set1',
                'Mean': 'viridis',
                },
            annotate_samples={
                'log_virus_reads_per_million': 'viridis',
                'time': 'Set1',
                },
            z_score=0,
            cmap='viridis',
            colorbars=True,
            vmin=-2,
            vmax=+2,
            )

    #plt.subplots_adjust(left=-0.15, right=0.8, bottom=0.06, top=1.2)

    print('t-SNE of the cell cycle genes')
    vs = dsn.dimensionality.tsne(perplexity=40)
    fig, axs = plt.subplots(
            1, 3, figsize=(11, 4),
            sharex=True, sharey=True,
            )
    for ax, cby in zip(axs, ['log_virus_reads_per_million', 'MCM6', 'CCNB1']):
        dsn.plot.scatter_reduced_samples(
                vs,
                color_by=cby,
                ax=ax,
                alpha=0.5,
                s=10,
                )
        ax.set_title(cby)

    print('Print again, only early time point')
    dsn4 = dsn.query_samples_by_metadata('time == 4')
    vs = dsn4.dimensionality.tsne(perplexity=40)
    fig, axs = plt.subplots(
            1, 3, figsize=(11, 4),
            sharex=True, sharey=True,
            )
    for ax, cby in zip(axs, ['log_virus_reads_per_million', 'MCM6', 'CCNB1']):
        dsn.plot.scatter_reduced_samples(
                vs,
                color_by=cby,
                ax=ax,
                alpha=0.5,
                s=10,
                )
        ax.set_title(cby)


    plt.ion()
    plt.show()
