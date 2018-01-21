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


# Classes / functions
class Fig2F(Dataset):
    def plot_time_correlation_COPE(self, cot, yerr=None):
        gid = self.featuresheet.loc[self.featuresheet['GeneName'] == 'COPE'].index[0]

        fig, axs = plt.subplots(
                nrows=1, ncols=4,
                sharex=True, sharey=True,
                figsize=(13, 2.7))

        st = self.split(phenotypes=['time'])
        for iax, (ax, t) in enumerate(zip(axs, sorted(st.keys(), key=int))):
            x = st[t].samplesheet['virus_reads_per_million']
            x = np.maximum(x, 0.5)

            y = st[t].counts.loc[[gid]].copy().unlog().iloc[0]
            y = np.maximum(y, 0.5)
            ax.scatter(
                    x, y,
                    s=15,
                    alpha=0.4,
                    edgecolor='none',
                    zorder=5)

            if t in ():
                loc = (0.02, 0.02)
                ha = 'left'
                va = 'bottom'
            else:
                loc = (0.98, 0.98)
                ha = 'right'
                va = 'top'

            if yerr is None:
                txt = '$\\rho = {:1.2f}$'.format(cot.loc[gid, int(t)])
            else:
                txt = '$\\rho = {:1.2f}({:.0f})$'.format(
                        cot.loc[gid, int(t)],
                        100 * yerr.loc[gid, int(t)],
                        )
            ax.text(*loc, txt,
                    ha=ha,
                    va=va,
                    transform=ax.transAxes,
                    zorder=6)

            from matplotlib.patches import Rectangle
            r = Rectangle(
                    (0.5, 0.85), 0.48, 0.13,
                    facecolor='white',
                    edgecolor='none',
                    transform=ax.transAxes,
                    zorder=4)
            ax.add_patch(r)

            ax.set_xlim(100, 1e6)
            ax.set_ylim(10, 1e4)
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True)
            ax.set_title(t+' hrs')

        fig.text(0.45, 0.025,
                 'Dengue per million transcripts',
                 ha='center')
        fig.text(0.02, 0.52,
                 'COPE counts per million transcripts',
                 ha='center',
                 va='center',
                 rotation=90)

        plt.tight_layout(rect=(0.03, 0.02, 1, 1), w_pad=1.5)

        return {'fig': fig, 'ax': ax}


# Script
if __name__ == '__main__':

    print('Load dataset')
    ds = Fig2F(
            counts_table='dengue',
            samplesheet='virus',
            featuresheet='humanGC38',
            )
    ds.query_samples_by_counts('total >= 50000', inplace=True)
    ds.samplesheet.rename(columns={'time [h]': 'time'}, inplace=True)
    cov = ds.samplesheet['coverage'] = ds.counts.sum(axis=0)

    print('Normalize')
    ds.counts.normalize('counts_per_million', inplace=True)

    print('Add virus reads')
    n = ds.samplesheet['numberDengueReads'].astype(int)
    ds.samplesheet['virus_reads_per_million'] = 1e6 * n / (cov + n)

    print('Log counts')
    ds.counts.log(inplace=True)

    print('Get correlations')
    co = ds.correlation.correlate_features_phenotypes(
            phenotypes='virus_reads_per_million',
            fillna=0).fillna(0)

    print('Get time correlations for COPE')
    dsn = ds.copy()
    gid_COPE = dsn.featuresheet.loc[dsn.featuresheet['GeneName'] == 'COPE'].index[0]
    dsn.counts = dsn.counts.loc[[gid_COPE]]
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

    print('Bootstrap over samples')
    import xarray as xr
    n_bootstraps = 100
    dboot = dsn.copy()
    dboot.counts = dboot.counts.loc[[gid_COPE]]
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

    print('Plot time correlations for COPE')
    d = dsn.plot_time_correlation_COPE(cotimes, yerr=cotimes_std)

    plt.ion()
    plt.show()
