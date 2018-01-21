'''
author:     Fabio Zanini
date:       08/09/17
content:    Script for figure 3.
'''
# Modules
import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

os.environ['SINGLET_CONFIG_FILENAME'] = 'singlet.yml'
sys.path.append('/home/fabio/university/postdoc/singlet')
from singlet.dataset import Dataset


# Globals
genes_validation = (
        'HSPA5',
        'NUCB2',
        'SEC61B',
        'SLC9A3R1',
        'CTSC',
        'RAB5A',
        'SSR3',
        'RPL31',
        'ISG15',
        'COPE',
        'SEC11C',
        'SEC61A1',
        'SELENOK',
        'ID2',
        'RSRC2',
        'PLPP5',
        'CALU',
        'SEC13',
        'TRAM1',
        'SLC38A2',
        'DDIT3',
        'ISG20',
        'RPN1',
        'FKBP11',
        'SCFD1',
        'WDFY1',
        'GORAB',
        'SPCS2',
        'TMED2',
        'CCND1',
        'RAB1B',
        'SSR1',
        'IRAK2',
        'CTNNB1')

switchers = (
        'SLC25A5',
        'MYL6',
        'DPYSL2',
        'HM13',
        'COPE',
        'PFN1',
        'GORASP2',
        'SNRPD2',
        'ANP32A',
        'UBC',
        'PSMB4',
        'SQSTM1',
        'RPN1',
        'TOMM5',
        'LAMTOR4',
        'NDUFA4',
        'PSMB3')

co_incg = (
        'HM13',
        'COPE',
        'GORASP2',
        'UBC',
        'SQSTM1',
        'RPN1')
co_decg = (
        'SLC25A5',
        'MYL6',
        'DPYSL2',
        'SNRPD2',
        'ANP32A',
        'PSMB3',
        'PSMB4',
        'TOMM5',
        'LAMTOR4',
        'NDUFA4',
        'PFN1')


# Classes / functions
class Fig3(Dataset):
    pass


# Script
if __name__ == '__main__':

    print('Load dataset')
    ds = Fig3(
            counts_table='dengueAndZika',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)
    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

    print('Normalize counts')
    ds.counts.normalize('counts_per_million', inplace=True)

    print('Add virus counts')
    ds.samplesheet['virus_reads_per_million'] = 0
    for virus in ('dengue', 'zika'):
        ind = ds.samplesheet['virus'] == virus
        n = ds.samplesheet.loc[ind, 'number'+virus.capitalize()+'Reads'].astype(int)
        ds.samplesheet.loc[ind, 'virus_reads_per_million'] = 1e6 * n / (cov.loc[ind] + n)

    print('Log counts')
    ds.counts.log(inplace=True)

    print('Select only some cells for comparison')
    dsc = ds.copy()
    dsc.samplesheet = dsc.samplesheet.query('500 < virus_reads_per_million')

    print('Get correlations')
    dsv = dsc.split(phenotypes='virus')
    vs = []
    cos = []
    for virus, dsvi in dsv.items():
        co = dsvi.correlation.correlate_features_phenotypes(
                phenotypes='virus_reads_per_million',
                fillna=0).fillna(0)
        cos.append(co)
        vs.append(virus)
    cos = pd.concat(cos, axis=1)
    cos.columns = pd.Index(vs, name='virus')

    print('Bootstrap examples')
    examples = ['ID2', 'ATF3', 'ACTG1', 'HSPA5']
    gids = [dsc.featuresheet.index[dsc.featuresheet['GeneName'] == gname][0] for gname in examples]
    dboot = dsc.copy()
    dboot.counts = ds.counts.loc[gids]
    n_bootstraps = 100
    co_mean = []
    co_std = []
    for virus, dbvi in dboot.split(phenotypes='virus').items():
        co_booti = []
        for iboot in range(n_bootstraps):
            dbi = dboot.bootstrap()
            coi = dbi.correlation.correlate_features_phenotypes(
                    phenotypes='virus_reads_per_million',
                    fillna=0).fillna(0)
            coi.name = iboot
            co_booti.append(coi)
        co_booti = pd.concat(co_booti, axis=1)
        co_booti.columns.name = 'bootstrap#'
        co_meani = co_booti.mean(axis=1)
        co_meani.name = virus
        co_stdi = co_booti.std(axis=1)
        co_stdi.name = virus
        co_mean.append(co_meani)
        co_std.append(co_stdi)
    co_mean = pd.concat(co_mean, axis=1)
    co_std = pd.concat(co_std, axis=1)

    print('Plot examples')
    fig, axs = plt.subplots(
            nrows=2, ncols=4,
            sharex=True,
            sharey=True,
            figsize=(10.2, 5.3))
    axs = axs.ravel()

    genes = ['ID2', 'ATF3', 'ACTG1', 'HSPA5']
    colors = ['steelblue', 'green', 'darkred', 'purple']
    markers = {'dengue': 'o', 'zika': 's'}
    for ig, gene in enumerate(genes):
        gid = dsc.featuresheet.loc[dsc.featuresheet['GeneName'] == gene].index[0]
        for iv, virus in enumerate(['dengue', 'zika']):
            dsvi = dsv[virus]
            ax = axs[ig * 2 + iv]

            x = dsvi.samplesheet['virus_reads_per_million']
            x = np.maximum(x, 0.5)
            y = dsvi.counts.loc[[gid]].copy().unlog().iloc[0]
            y = np.maximum(y, 0.5)
            ax.scatter(
                    x, y,
                    s=15,
                    alpha=0.4,
                    edgecolor='none',
                    facecolor=colors[ig],
                    marker=markers[virus],
                    zorder=5)

            loc = (0.985, 0.98)
            ha = 'right'
            va = 'top'
            ax.text(*loc, '$\\rho = {:1.2f}({:.0f})$'.format(
                cos.loc[gid, virus],
                100 * co_std.loc[gid, virus]),
                    ha=ha,
                    va=va,
                    transform=ax.transAxes)

            ax.text(0.02, 0.98, virus[0].upper(),
                    ha='left',
                    va=va,
                    transform=ax.transAxes,
                    fontweight='bold',
                    fontsize=12)

            ax.set_xlim(400, 1e6)
            ax.set_ylim(0.4, 1e5)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.grid(True)

        fig.text(
                0.29 + 0.46*(ig % 2), 0.98 - 0.45*(ig >= 2),
                gene,
                ha='center',
                va='top',
                fontweight='bold')

    fig.text(0.5, 0.022,
             'Virus per million transcripts',
             ha='center')
    fig.text(0.02, 0.5,
             'Counts per million transcripts',
             ha='center',
             va='center',
             rotation=90)

    plt.tight_layout(rect=(0.03, 0.03, 1, 0.98), h_pad=1.3, w_pad=0.5)

    plt.ion()
    plt.show()
