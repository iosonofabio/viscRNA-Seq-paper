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
class Fig4(Dataset):
    pass


# Script
if __name__ == '__main__':

    ds = Fig4(
            counts_table='dengueAndZika',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)

    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)
    ds.counts.normalize('counts_per_million', inplace=True)
    ds.samplesheet['virus_reads_per_million'] = 0
    for virus in ('dengue', 'zika'):
        ind = ds.samplesheet['virus'] == virus
        n = ds.samplesheet.loc[ind, 'number'+virus.capitalize()+'Reads'].astype(int)
        ds.samplesheet.loc[ind, 'virus_reads_per_million'] = 1e6 * n / (cov.loc[ind] + n)
    ds.counts.log(inplace=True)

    cos_fn = '../../tables/genes_fig4.tsv'
    if os.path.isfile(cos_fn):
        features = pd.read_csv(cos_fn, sep='\t', index_col=0).index
    else:
        # Select only some cells for comparison
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
        features = cos.index[np.abs(cos).max(axis=1) >= 0.4]
        with open(cos_fn, 'wt') as f:
            f.write('\n'.join([features.name] + features.tolist()))

    # Plot dimensionality reductions for all genes at the ourskirts of the correlation plot
    dsd = ds.copy()
    dsd.counts = dsd.counts.loc[features]
    dsd.samplesheet['virus_polarized'] = np.log10(1 + dsd.samplesheet['virus_reads_per_million'])
    ind_zika = dsd.samplesheet['virus'] == 'zika'
    dsd.samplesheet.loc[ind_zika, 'virus_polarized'] = -dsd.samplesheet.loc[ind_zika, 'virus_polarized']
    dsd.samplesheet['time_polarized'] = dsd.samplesheet['time'].astype(int)
    dsd.samplesheet.loc[ind_zika, 'time_polarized'] = -dsd.samplesheet.loc[ind_zika, 'time_polarized']
    # t-SNE is not deterministic, so save to file to make sure we can reproduce
    # the paper figures with variations if needed
    vs_fn = '../../tables/tsne.tsv'
    if os.path.isfile(vs_fn):
        vs = pd.read_csv(vs_fn, sep='\t', index_col=0)
    else:
        raise ValueError("t-SNE file not found?")
        #vs = dsd.dimensionality.tsne(perplexity=10)
        #vs.to_csv(vs_fn, sep='\t')
    # Original cmap
    #cmap = sns.diverging_palette(220, 20, as_cmap=True, center='light')
    from matplotlib import cm
    cmap  = cm.RdYlBu_r
    fig, axs = plt.subplots(
            nrows=2, ncols=2, figsize=(10, 6),
            gridspec_kw={'height_ratios': [1, 20]})
    axs = axs[::-1].T.ravel()
    dsd.plot.scatter_reduced_samples(
            vs,
            color_by='virus_polarized',
            cmap=cmap,
            alpha=0.8,
            ax=axs[0],
            s=10,
            lw=.25,
            edgecolor=[0.4, 0.4, 0.4],
            zorder=10)
    dsd.plot.scatter_reduced_samples(
            vs,
            color_by='time_polarized',
            cmap=cmap,
            alpha=0.8,
            ax=axs[2],
            s=10,
            lw=.25,
            edgecolor=[0.4, 0.4, 0.4],
            zorder=10)
    axs[2].set_ylabel('')
    axs[2].set_yticklabels([])
    norm1 = mpl.colors.Normalize(
            vmin=dsd.samplesheet['virus_polarized'].values.min(),
            vmax=dsd.samplesheet['virus_polarized'].values.max())
    cb1 = mpl.colorbar.ColorbarBase(
            ax=axs[1], cmap=cmap, norm=norm1, orientation='horizontal')
    cb1.set_ticks([-4, -2, 0, 2, 4])
    cb1.set_ticklabels(
            ['$10^4$ Zika cpmt', '$10^2$ Zika cpmt', 'no virus',
             '$10^2$ dengue cpmt', '$10^4$ dengue cpmt'])
    for tkl in cb1.ax.get_xticklabels():
        tkl.set_rotation(45)
        tkl.set_fontsize(8)
    cb2 = mpl.colorbar.ColorbarBase(
            ax=axs[3], cmap=cmap,
            boundaries=[-48, -24, -12, -4, 0, 4, 12, 24, 48],
            values=[-48, -24, -12, -4, 4, 12, 24, 48],
            ticks=[-36, -18, -8, -2, 2, 8, 18, 36],
            orientation='horizontal')
    cb2.set_ticklabels(
            ['Zika 48 hrs', 'Zika 24 hrs', 'Zika 12 hrs', 'Zika 4 hrs',
             'dengue 4 hrs', 'dengue 12 hrs', 'dengue 24 hrs', 'dengue 48 hrs'])
    for tkl in cb2.ax.get_xticklabels():
        tkl.set_rotation(45)
        tkl.set_fontsize(8)
    axs[0].grid(False)
    axs[2].grid(False)

    # Plot mean of samples at increasing viral load
    bins = [0, 30, 200, 1000, 10000, 70000, 500000]
    for virus in ['dengue', 'zika']:
        ph = dsd.samplesheet['virus_reads_per_million']
        vi = dsd.samplesheet['virus']
        vms = pd.concat([vs.loc[
            (ph >= bins[i]) & (ph < bins[i+2]) & (vi == virus)
            ].mean(axis=0) for i in range(len(bins) - 2)],
            axis=1).T
        x = vms.values[:, 0]
        y = vms.values[:, 1]
        axs[0].quiver(
                x[:-1], y[:-1], x[1:]-x[:-1], y[1:]-y[:-1],
                scale_units='xy', angles='xy', scale=1,
                zorder=20)

    plt.tight_layout(h_pad=0.5, w_pad=0.1)

    # Plot example genes
    gnames = ['ID2', 'ATF3', 'ACTG1', 'HSPA5']
    cmap = 'plasma'
    fig, axs = plt.subplots(
            sharex=True, sharey=True,
            nrows=1, ncols=4, figsize=(10, 3))
    for gname, ax in zip(gnames, axs):
        gid = dsd.featuresheet.index[dsd.featuresheet['GeneName'] == gname][0]
        dsd.plot.scatter_reduced_samples(
                vs,
                color_by=gid,
                cmap=cmap,
                alpha=0.7,
                ax=ax,
                s=10,
                lw=.25,
                edgecolor=[0.4, 0.4, 0.4],
                zorder=10)
        ax.set_title(gname)
        ax.set_xlabel('')
    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    plt.tight_layout(rect=(0, 0.03, 1, 1), w_pad=0.1)

    # Plot apoptotic genes
    dsd = ds.copy()
    #gnames = ['CASP2', 'CASP3', 'CASP4', 'CASP6', 'CASP7', 'CASP8', 'CASP9', 'CASP10']
    gnames = ['CASP3', 'CASP9', 'BCL2L11', 'BCL2', 'MAPK8', 'TP53']
    gids = [ds.featuresheet.index[ds.featuresheet['GeneName'] == gname][0] for gname in gnames]
    dsd.counts = dsd.counts.loc[gids]
    cmap = 'plasma'
    fig, axs = plt.subplots(
            sharex=True, sharey=True,
            nrows=2, ncols=4, figsize=(11, 7))
    axs = axs.ravel()
    for gname, ax in zip(gnames, axs):
        gid = dsd.featuresheet.index[dsd.featuresheet['GeneName'] == gname][0]
        dsd.plot.scatter_reduced_samples(
                vs,
                color_by=gid,
                cmap=cmap,
                alpha=0.7,
                ax=ax,
                s=10,
                lw=.25,
                edgecolor=[0.4, 0.4, 0.4],
                zorder=10)
        ax.set_title(gname)
        ax.set_xlabel('')
        ax.set_ylabel('')
    fig.text(0.5, 0.02, 'dimension 1', ha='center')
    fig.text(0.02, 0.5, 'dimension 2', ha='right', va='center', rotation=90)
    plt.tight_layout(rect=(0.02, 0.03, 1, 1), w_pad=0.1)

    dsh = dsd.copy()
    dsh.query_samples_by_metadata('virus_reads_per_million > 5e3', inplace=True)
    fig, axs = plt.subplots(
            sharex=True, sharey=True,
            nrows=2, ncols=4, figsize=(11, 7))
    axs = axs.ravel()
    coh = dsh.correlation.correlate_features_phenotypes(['virus_reads_per_million'], fillna=0)
    y = dsh.samplesheet.loc[:, 'virus_reads_per_million'].values
    for gname, ax in zip(gnames, axs):
        gid = dsh.featuresheet.index[dsd.featuresheet['GeneName'] == gname][0]
        x = dsh.counts.loc[gid].values
        ax.scatter(
                x, y,
                alpha=0.5,
                s=10,
                zorder=10)
        ax.set_title(gname)
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.grid(True)
        ax.set_xlim(-1.1, 5.5)
        ax.set_ylim(3e3, 10**(5.5))
        ax.set_yscale('log')
    fig.text(0.5, 0.02, 'Log10 Counts per million transcripts', ha='center')
    fig.text(0.02, 0.5, 'Virus amount per million transcripts', ha='right', va='center', rotation=90)
    plt.tight_layout(rect=(0.02, 0.03, 1, 1), w_pad=0.1)

    plt.ion()
    plt.show()
