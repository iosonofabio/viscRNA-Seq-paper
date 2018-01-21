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
        #'ISG15',
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
        #'ISG20',
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
class Fig2G(Dataset):
    def plot_time_correlation_switchers(self, coswitch, gnames, coPN, yerr=None):
        # Summary plot with all lines
        x = coswitch.columns.astype(int)

        co_inc = np.array([self.featuresheet.index[(self.featuresheet['GeneName'] == g)][0] for g in co_incg])
        co_dec = np.array([self.featuresheet.index[(self.featuresheet['GeneName'] == g)][0] for g in co_decg])
        co_inc = coswitch.loc[co_inc]
        co_dec = coswitch.loc[co_dec]

        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(13, 5))
        colors = sns.color_palette()
        colors.append((0.2, 0.2, 0.2, 1))
        for ax, cot in zip(axs, [co_inc, co_dec]):
            for ig, (gid, y) in enumerate(cot.iterrows()):
                gname = gnames.loc[gid]
                if yerr is None:
                    ax.plot(x, y,
                            lw=2,
                            alpha=0.8,
                            label=gname,
                            color=colors[ig],
                            zorder=5-0.01*ig,
                            )
                else:
                    dy = yerr.loc[gid]
                    ax.errorbar(
                            x, y,
                            yerr=dy,
                            lw=2,
                            alpha=0.8,
                            label=gname,
                            color=colors[ig],
                            zorder=5-0.01*ig,
                            )
            ax.set_xlim(0, 50)
            ax.set_ylim(-0.6, 0.6)
            ax.axhline(0, color=(0.3, 0.3, 0.3, 0.5), lw=2)
            ax.grid(True)

        fig.text(0.47, 0.02, 'Time post infection [hrs]', ha='center')
        axs[0].set_ylabel('Correlation with dengue virus')
        axs[0].legend(loc='lower right', bbox_to_anchor=(1.01, -0.01), title='Genes', fontsize=8)
        axs[1].legend(loc='upper right', bbox_to_anchor=(1.01, 1.015), title='Genes', fontsize=8)
        axs[0].set_title('Anticorrelated $\\rightarrow$ correlated')
        axs[1].set_title('Correlated $\\rightarrow$ anticorrelated')

        # Add a table about gene statistics
        ax = axs[0]
        table = coPN.groupby(['positive', 'negative']).count()['all'].unstack()
        cellText = [[str(cell) for cell in row] for key, row in table.iterrows()]
        ax.table(
                cellText=cellText,
                cellColours=[['w', 'w'], ['w', (0, 1, 0, 0.1)]],
                bbox=(0.0754, 0.732, 0.25, 0.2), zorder=10)
        ax.text(0.21, 0.955, '$\\rho < -0.3$', ha='center', transform=ax.transAxes, zorder=10)
        ax.text(0.04, 0.87, '$\\rho > 0.3$', ha='center', transform=ax.transAxes, zorder=11, rotation=90)
        from matplotlib.patches import Rectangle
        r1 = Rectangle((0.1, 0.201), 19.9, 0.395, edgecolor='none', facecolor='white', zorder=9)
        ax.add_patch(r1)

        plt.tight_layout(rect=(0, 0.04, 0.98, 1))

        return {'fig': fig, 'ax': ax}


class Fig4C(Dataset):
    def plot_time_correlation_switchers(self, coswitch, gnames, coPN, co_incg, co_decg, yerr=None):
        # Summary plot with all lines
        x = coswitch.columns.astype(int)

        co_inc = np.array([self.featuresheet.index[(self.featuresheet['GeneName'] == g)][0] for g in co_incg])
        co_dec = np.array([self.featuresheet.index[(self.featuresheet['GeneName'] == g)][0] for g in co_decg])
        co_inc = coswitch.loc[co_inc]
        co_dec = coswitch.loc[co_dec]

        fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 4))
        colors = sns.color_palette()
        colors.append((0.2, 0.2, 0.2, 1))
        for ax, cot in zip(axs, [co_inc, co_dec]):
            for ig, (gid, y) in enumerate(cot.iterrows()):
                gname = gnames.loc[gid]
                if yerr is None:
                    ax.plot(
                            x, y,
                            lw=2,
                            alpha=0.8,
                            label=gname,
                            color=colors[ig],
                            zorder=5-0.01*ig,
                            )
                else:
                    dy = yerr.loc[gid]
                    ax.errorbar(
                            x, y,
                            yerr=dy,
                            lw=2,
                            alpha=0.8,
                            label=gname,
                            color=colors[ig],
                            zorder=5-0.01*ig,
                            )
            ax.set_xlim(0, 50)
            ax.set_ylim(-0.8, 0.8)
            ax.axhline(0, color=(0.3, 0.3, 0.3, 0.5), lw=2)
            ax.grid(True)

        fig.text(0.47, 0.02, 'Time post infection [hrs]', ha='center')
        axs[0].set_ylabel('Correlation with Zika virus')
        axs[0].legend(loc='lower right', bbox_to_anchor=(1.01, -0.01), title='Genes', fontsize=8)
        axs[1].legend(loc='upper right', bbox_to_anchor=(1.01, 1.015), title='Genes', fontsize=7, ncol=2)
        axs[0].set_title('Anticorrelated $\\rightarrow$ correlated')
        axs[1].set_title('Correlated $\\rightarrow$ anticorrelated')

        # Add a table about gene statistics
        ax = axs[0]
        table = coPN.groupby(['positive', 'negative']).count()['all'].unstack()
        cellText = [[str(cell) for cell in row] for key, row in table.iterrows()]
        ax.table(
                cellText=cellText,
                cellColours=[['w', 'w'], ['w', (0, 1, 0, 0.1)]],
                bbox=(0.0754, 0.732, 0.25, 0.2), zorder=10)
        ax.text(0.21, 0.955, '$\\rho < -0.3$', ha='center', transform=ax.transAxes, zorder=10)
        ax.text(0.04, 0.87, '$\\rho > 0.3$', ha='center', transform=ax.transAxes, zorder=11, rotation=90)
        from matplotlib.patches import Rectangle
        r1 = Rectangle((0.1, 0.201), 19.9, 0.595, edgecolor='none', facecolor='white', zorder=9)
        ax.add_patch(r1)

        plt.tight_layout(rect=(0, 0.04, 0.98, 1))

        return {'fig': fig, 'ax': ax}


# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument(
            'panels', nargs='+', choices=['2G', '4C'],
            default=['2G', '4C'],
            help='Panels to show')
    args = parser.parse_args()

    if '2G' in args.panels:
        print('Load dengue dataset')
        ds = Fig2G(
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

        #print('Get correlations')
        #co = ds.correlation.correlate_features_phenotypes(
        #        phenotypes='virus_reads_per_million',
        #        fillna=0).fillna(0)

        print('Get time correlations for somewhat expressed genes')
        dsn = ds.copy()
        dsn.counts = dsn.counts.loc[(dsn.counts >= 1).sum(axis=1) >= 10]
        dstimes = dsn.split(phenotypes=['time'])
        cotimes = []
        for t in sorted(dstimes.keys(), key=int):
            dst = dstimes[t]
            co = dst.correlation.correlate_features_phenotypes(
                    phenotypes='virus_reads_per_million',
                    fillna=0).fillna(0)
            co.name = int(t)
            cotimes.append(co)
        cotimes = pd.concat(cotimes, axis=1)
        cotimes.columns.name = 'time'

        print('Get time switchers')
        coP = (cotimes > 0.3).any(axis=1)
        coN = (cotimes < -0.3).any(axis=1)
        coPN = pd.concat([coP, coN], axis=1)
        coPN.columns = ['positive', 'negative']
        coPN['all'] = 1
        switch = coPN.groupby(['positive', 'negative']).get_group((True, True)).index
        gnames = dsn.featuresheet.loc[switch, 'GeneName']

        print('Bootstrap over samples after selecting only switcher genes')
        import xarray as xr
        n_bootstraps = 100
        dboot = dsn.copy()
        dboot.counts = dboot.counts.loc[switch]
        cotimes_boot = []
        for iboot in range(n_bootstraps):
            cot = []
            dbi = dboot.bootstrap()
            dbitimes = dbi.split(phenotypes=['time'])
            for t in sorted(dbitimes.keys(), key=int):
                dst = dbitimes[t]
                co = dst.correlation.correlate_features_phenotypes(
                        phenotypes='virus_reads_per_million',
                        fillna=0).fillna(0)
                co.name = int(t)
                cot.append(co)
            cot = pd.concat(cot, axis=1)
            cot.columns.name = 'time'
            cotimes_boot.append(cot)
        cotimes_boot = xr.concat(
                [xr.DataArray(cot) for cot in cotimes_boot],
                dim=pd.Index(np.arange(n_bootstraps), name='bootstrap#'),
                )

        cotimes_mean = cotimes_boot.mean(dim='bootstrap#').to_dataframe('ciao')['ciao'].unstack()
        cotimes_std = cotimes_boot.std(dim='bootstrap#').to_dataframe('ciao')['ciao'].unstack()

        print('Plot switchers with cell boostrapping')
        d = dsn.plot_time_correlation_switchers(
                cotimes,
                gnames,
                coPN,
                yerr=cotimes_std,
                )

    if '4C' in args.panels:
        print('Load Zika dataset')
        ds = Fig4C(
                counts_table='zika',
                samplesheet='virus',
                featuresheet='humanGC38',
                )
        ds.query_samples_by_counts('total >= 50000', inplace=True)
        ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)

        cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)
        ds.counts.normalize('counts_per_million', inplace=True)

        n = ds.samplesheet['numberZikaReads'].astype(int)
        ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)
        ds.counts.log(inplace=True)

        print('Get time correlations for somewhat expressed genes')
        dsn = ds.copy()
        dsn.counts = dsn.counts.loc[(dsn.counts >= 1).sum(axis=1) >= 10]
        dstimes = dsn.split(phenotypes=['time'])
        cotimes = []
        for t in sorted(dstimes.keys(), key=int):
            dst = dstimes[t]
            co = dst.correlation.correlate_features_phenotypes(
                    phenotypes='virus_reads_per_million',
                    fillna=0).fillna(0)
            co.name = int(t)
            cotimes.append(co)
        cotimes = pd.concat(cotimes, axis=1)
        cotimes.columns.name = 'time'

        print('Get time switchers')
        coP = (cotimes > 0.3).any(axis=1)
        coN = (cotimes < -0.3).any(axis=1)
        coPN = pd.concat([coP, coN], axis=1)
        coPN.columns = ['positive', 'negative']
        coPN['all'] = 1
        switch = coPN.groupby(['positive', 'negative']).get_group((True, True)).index
        gnames = dsn.featuresheet.loc[switch, 'GeneName']
        co_inc = switch[cotimes.loc[switch, 48] > 0]
        co_dec = switch[cotimes.loc[switch, 48] < 0]
        co_incg = dsn.featuresheet.loc[co_inc, 'GeneName'].values
        co_decg = dsn.featuresheet.loc[co_dec, 'GeneName'].values

        #print('Plot switchers')
        #dsn.plot_time_correlation_switchers(cotimes.loc[switch], gnames, coPN, co_incg, co_decg)

        print('Bootstrap over samples after selecting only switcher genes')
        import xarray as xr
        n_bootstraps = 100
        dboot = dsn.copy()
        dboot.counts = dboot.counts.loc[switch]
        cotimes_boot = []
        for iboot in range(n_bootstraps):
            cot = []
            dbi = dboot.bootstrap()
            dbitimes = dbi.split(phenotypes=['time'])
            for t in sorted(dbitimes.keys(), key=int):
                dst = dbitimes[t]
                co = dst.correlation.correlate_features_phenotypes(
                        phenotypes='virus_reads_per_million',
                        fillna=0).fillna(0)
                co.name = int(t)
                cot.append(co)
            cot = pd.concat(cot, axis=1)
            cot.columns.name = 'time'
            cotimes_boot.append(cot)
        cotimes_boot = xr.concat(
                [xr.DataArray(cot) for cot in cotimes_boot],
                dim=pd.Index(np.arange(n_bootstraps), name='bootstrap#'),
                )

        cotimes_mean = cotimes_boot.mean(dim='bootstrap#').to_dataframe('ciao')['ciao'].unstack()
        cotimes_std = cotimes_boot.std(dim='bootstrap#').to_dataframe('ciao')['ciao'].unstack()

        print('Plot switchers with cell boostrapping')
        dsn.plot_time_correlation_switchers(
                cotimes,
                gnames,
                coPN,
                co_incg, co_decg,
                yerr=cotimes_std,
                )

    plt.ion()
    plt.show()
