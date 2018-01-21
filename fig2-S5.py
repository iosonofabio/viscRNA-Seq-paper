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
genes_by_complex = {
        'SEC61': ('SEC61A1', 'SEC61A2', 'SEC61B', 'SEC61G'),
        'TRAP': ('SSR1', 'SSR2', 'SSR3', 'SSR4'),
        'SRP': ('SRP9', 'SRP14'),
        'SPCS': ('SEC11A', 'SEC11C', 'SPCS1', 'SPCS2', 'SPCS3'),
        #'COPI': ('COPA', 'COPB1', 'COPB2', 'COPG1', 'COPG2',
        #         'ARCN1',  # = delta-COP
        #         'COPE', 'COPZ1', 'COPZ2'),
        'OST': ('STT3A', 'STT3B',
                'OSTC', 'OST4', 'DDOST',
                'RPN1', 'RPN2',
                'DAD1', 'MAGT1', 'KRTCAP2',
                'TUSC3', 'TMEM258'),
        'RPL': ('RPL31', 'RPL38'),
        'other': ('HM13', 'TRAM1'),
        }


# Classes / functions
class Fig2(Dataset):
    def plot_correlations_complexes(self, co, genes_dict, figsize=(10, 18)):
        fig, axs = plt.subplots(
                nrows=8, ncols=4,
                figsize=figsize)
        axs = axs.ravel()

        complexes = ['SEC61', 'TRAP', 'SRP', 'SPCS', 'OST', 'RPL', 'other']
        color_rims = sns.color_palette(n_colors=len(complexes))

        #TODO: separate by complex
        genes = sum((list(genes_dict[c]) for c in complexes), [])

        # Find IDs
        gids = [ds.featuresheet.loc[ds.featuresheet['GeneName'] == g].index[0]
                for g in genes]

        thresholds = {
            'SEC61A1': 3e3,
            'RPL31': 3e1,
            'TRAM1': 3e3,
            'SEC11C': 4e3,
            'SSR3': 4e3,
            }

        x = ds.samplesheet['virus_reads_per_million']
        x = np.maximum(x, 0.5)
        for iax, (gid, ax) in enumerate(zip(gids, axs)):
            gname = ds.featuresheet.loc[gid].iloc[0]
            y = ds.counts.loc[[gid]].copy().unlog().iloc[0]
            y = np.maximum(y, 0.5)
            ax.scatter(
                    x, y,
                    s=15,
                    alpha=0.4,
                    edgecolor='none',
                    zorder=5)

            # KDEplot is finnicky because of the log transforms, so I do it by
            # hand
            h2d, bx, by = np.histogram2d(
                    x, y,
                    bins=[0.4] + list(np.logspace(1, 6, 16)[:-1]),
                    normed=False)

            from matplotlib.patches import Rectangle
            from matplotlib import cm
            hmin = 9e-1
            h2d = np.maximum(h2d, hmin)
            hmax = h2d.max()
            cmap = cm.get_cmap('Greens')
            for ix in range(len(bx) - 1):
                for iy in range(len(by) - 1):
                    val = h2d[ix, iy]
                    lex, rex = bx[ix: ix+2]
                    ley, rey = by[iy: iy+2]
                    vcol = cmap(1.0 - 1.0 * np.log(val / hmax) / np.log(hmin / hmax))
                    r = Rectangle(
                            (lex, ley), rex - lex, rey - ley,
                            facecolor=vcol,
                            edgecolor='none',
                            alpha=0.5,
                            zorder=2)
                    ax.add_patch(r)

            loc = (0.985, 0.98)
            ha = 'right'
            va = 'top'
            ax.text(*loc, '$\\rho = {:1.2f}$'.format(co.loc[gid]),
                    ha=ha,
                    va=va,
                    transform=ax.transAxes)

            if gname in thresholds:
                th = thresholds[gname]
                if th is not None:
                    ax.axvline(th, lw=2, color='darkred')

            for colrim, comp in zip(color_rims, complexes):
                if gname in genes_dict[comp]:
                    break
            else:
                raise IndexError('Gene not in any complex!')
            rim = Rectangle(
                    (0.01, 0.01), 0.98, 0.98,
                    facecolor='none',
                    edgecolor=colrim,
                    lw=3,
                    zorder=70,
                    transform=ax.transAxes)
            ax.add_patch(rim)

            # FIXME: this does not quite work yet
            #sns.kdeplot(
            #        x, y,
            #        zorder=3,
            #        n_levels=40,
            #        shade=True,
            #        ax=ax)

            ax.set_xlim(100, 1e6)
            ax.set_ylim(0.4, 1e5)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.grid(True)
            ax.set_title(gname)

        if len(axs) > len(genes):
            for i in range(len(genes), len(axs)):
                axs[i].set_axis_off()

        # Legend in last axes
        from matplotlib.patches import Patch
        ax = axs[-1]
        patches = [Patch(color=cr, label=cp) for (cr, cp) in zip(color_rims, complexes)]
        ax.legend(handles=patches, loc='center',
                  title='Complex')

        fig.text(0.5, 0.022,
                 'Dengue per million transcripts',
                 ha='center')
        fig.text(0.02, 0.5,
                 'Counts per million transcripts',
                 ha='center',
                 rotation=90)

        plt.tight_layout(rect=(0.03, 0.03, 1, 1))



# Script
if __name__ == '__main__':

    ds = Fig2(
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

    print('Get correlations')
    co = ds.correlation.correlate_features_phenotypes(
            phenotypes='virus_reads_per_million',
            fillna=0).fillna(0)

    ds.plot_correlations_complexes(co, genes_by_complex, figsize=(10, 18))


    plt.ion()
    plt.show()
