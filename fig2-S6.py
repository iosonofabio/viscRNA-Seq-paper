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
from singlet.counts_table import CountsTable


# Classes / functions
class CompareOldNewDengue(Dataset):
    pass


# Script
if __name__ == '__main__':

    ## Check ERCCs and others
    #ct = CountsTable.from_tablename('old_dengue')
    ## LOOK OK

    print('Load dataset')
    ds = CompareOldNewDengue(
            counts_table='both_dengue',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)

    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)

    # Number of cells
    nc = (ds.samplesheet
            .groupby(['experiment', 'time', 'MOI'])
            .count()
            .iloc[:, 0]
            .unstack('experiment')
            .sort_index(axis=0, level=('time', 'MOI'))
            )
    print(nc)

    print('Calculate coverage')
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

    print('Normalize counts')
    ds.counts.normalize('counts_per_million', inplace=True)

    print('Get and normalize number of virus reads')
    n = ds.samplesheet['numberDengueReads'].astype(int)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)

    print('Select genes with decent expression')
    indg = ((ds.counts >= 50).sum(axis=1) >= 20)
    ds.counts = ds.counts.loc[indg]

    print('Split experiments')
    dspl = ds.split(phenotypes='experiment')

    print('Get correlations by experiment')
    cos = {}
    for exp, dsi in dspl.items():
        co = dsi.correlation.correlate_features_phenotypes(
                phenotypes='virus_reads_per_million',
                fillna=0).fillna(0)
        cos[exp] = co
    cos = pd.concat([cos['10017003'], cos['10017006']], axis=1)
    cos.columns = ['10017003', '10017006']

    print('Translate relevant correlations')
    cosg = cos.copy()
    ind = (np.abs(cosg) > 0.3).any(axis=1)
    indi = ind[ind].index
    cosg = cosg.loc[ind]
    cosg.index = pd.Index(
            ds.featuresheet.loc[indi]['GeneName'].values,
            name='GeneName',
            )

    #print('Log counts')
    #ds.counts.log(inplace=True)

    print('Compare correlations in the two experiments')
    from scipy.stats import pearsonr
    corrg = pearsonr(cosg['10017003'].values, cosg['10017006'].values)
    corr = pearsonr(cos['10017003'].values, cos['10017006'].values)

    print('Scatter correlations in the two experiments')
    fig, ax = plt.subplots(figsize=(5, 4))
    cos.plot(
            kind='scatter',
            x='10017003',
            y='10017006',
            s=20,
            ax=ax,
            grid=True,
            alpha=0.05,
            color='grey',
            zorder=2,
            label='All genes: {:.2f}'.format(corr[0]),
            )
    cosg.plot(
            kind='scatter',
            x='10017003',
            y='10017006',
            s=20,
            ax=ax,
            grid=True,
            color='steelblue',
            alpha=0.6,
            zorder=3,
            label='$|\\rho| \geq 0.3$: {:.2f}'.format(corrg[0]),
            )
    xline = np.linspace(-0.7, 0.7, 100)
    ax.plot(xline, xline, lw=1.5, color='k', alpha=0.8, zorder=2)
    ax.legend(loc='best', title='Pearson\'s r:')
    ax.set_xlim(-0.62, 0.62)
    ax.set_ylim(-0.62, 0.62)
    ax.set_xlabel('Small scale experiment')
    ax.set_ylabel('Large scale experiment')

    ## Tag outliers
    #outliers = (np.abs(cosg) > 0.55).any(axis=1) | (np.abs(cosg.iloc[:, 1] - cosg.iloc[:, 0]) > 0.4)
    #outliers = outliers[outliers].index
    #for o in outliers:
    #    xo = cosg.loc[o, '10017003']
    #    yo = cosg.loc[o, '10017006']
    #    ax.text(xo + 0.02, yo - 0.03, o)

    plt.tight_layout()
    plt.ion()
    plt.show()

    ##print('Plot single correlation plots')
    ##ds.plot_single_correlations(co, figsize=(6.5, 4.6))
    #ds.plot_correlations_candidates(co, genes_validation, figsize=(10, 18))

    ##print('Get time correlations for somewhat expressed genes')
    ##dsn = ds.copy()
    ##dsn.counts = dsn.counts.loc[(dsn.counts >= 1).sum(axis=1) >= 10]
    ##dstimes = dsn.split(phenotypes=['time'])
    ##cotimes = []
    ##for t in sorted(dstimes.keys(), key=int):
    ##    dst = dstimes[t]
    ##    co = dst.correlation.correlate_features_phenotypes(
    ##            phenotypes='virus_reads_per_million',
    ##            fillna=0).fillna(0)
    ##    co.name = t
    ##    cotimes.append(co)
    ##cotimes = pd.concat(cotimes, axis=1)
    ##cotimes.columns.name = 'time'

    ##print('Get time switchers')
    ##coP = (cotimes > 0.3).any(axis=1)
    ##coN = (cotimes < -0.3).any(axis=1)
    ##coPN = pd.concat([coP, coN], axis=1)
    ##coPN.columns = ['positive', 'negative']
    ##coPN['all'] = 1
    ##switch = coPN.groupby(['positive', 'negative']).get_group((True, True)).index
    ##gnames = dsn.featuresheet.loc[switch, 'GeneName']

    ###print('Plot switchers')
    ###dsn.plot_time_correlation_switchers(cotimes.loc[switch], gnames, coPN)

    ##print('Plot time correlations for COPE')
    ##dsn.plot_time_correlation_COPE(cotimes, coPN)

    #plt.ion()
    #plt.show()
