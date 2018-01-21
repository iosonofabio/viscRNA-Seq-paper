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
class Fig2A_E(Dataset):
    def plot_correlation_distribution(
            self,
            co,
            figsize=None,
            yerr=None,
            cumulative=False):
        nL = co.nlargest(10)
        nS = co.nsmallest(10)
        fig, ax = plt.subplots(figsize=figsize)

        if not cumulative:
            sns.distplot(
                    co,
                    bins=(-0.7, -0.6, -0.5, -0.45, -0.4, -0.35, -0.3, -0.25, -0.2, -0.15, -0.1, -0.05, -0.02,
                          0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7),
                    kde_kws={'bw': 0.25},
                    hist=True,
                    kde=True,
                    ax=ax)
            sns.rugplot(
                    co.loc[np.concatenate([nL.index, nS.index])],
                    color='darkred',
                    lw=2,
                    ax=ax)
            ax.set_ylim(0.001, 100)
            ax.set_xlim(-0.7, 0.7)
            ax.set_yscale('log')
            ax.grid(True)
            ax.set_xlabel('Correlation with dengue virus')
            ax.set_ylabel('Density')

        else:
            x = np.sort(co.values)
            y = 1.0 - np.linspace(0, 1, len(x))
            ax.plot(x, y, lw=2, color='darkblue')

            for gn in nL.index:
                xi = co.loc[gn]
                ax.axvline(xi, ymin=0, ymax=0.05, color='darkred', lw=2)
            for gn in nS.index:
                xi = co.loc[gn]
                ax.axvline(xi, ymin=0.95, ymax=1, color='darkred', lw=2)

            ax.set_ylim(0.0001, 0.9999)
            ax.set_xlim(-0.7, 0.7)
            #ax.set_yticks([0.0001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.9999])
            #ax.set_yscale('logit')
            ax.grid(True)
            ax.set_xlabel('Correlation with dengue virus')
            ax.set_ylabel('Fraction genes with $\\rho \geq x$')

        plt.tight_layout()

        # Inset with highest correlations
        nLt = nL.copy()
        nSt = nS.copy()
        nLt.index = ds.featuresheet.loc[nL.index, 'GeneName']
        nSt.index = ds.featuresheet.loc[nS.index, 'GeneName']
        if yerr is None:
            ax.table(
                    cellText=[[key, '{:1.2f}'.format(val)] for key, val in nLt.items()],
                    cellColours=[['white'] * 2 for key in nLt],
                    rowLabels=None,
                    colLabels=None,
                    bbox=(0.7, 0.53, 0.28, 0.45),
                    zorder=5,
                    )
            ax.table(
                    cellText=[[key, '{:1.2f}'.format(val)] for key, val in nSt.items()],
                    cellColours=[['white'] * 2 for key in nSt],
                    rowLabels=None,
                    colLabels=None,
                    bbox=(0.02, 0.53 - 0.51 * cumulative, 0.28, 0.45),
                    zorder=5,
                    )
        else:
            stdLt = yerr.loc[nL.index].copy()
            stdLt.index = ds.featuresheet.loc[nL.index, 'GeneName']
            stdSt = yerr.loc[nS.index].copy()
            stdSt.index = ds.featuresheet.loc[nS.index, 'GeneName']
            ax.table(
                    cellText=[[key, '{:1.2f}({:.0f})'.format(val, 100 * stdLt[key])] for key, val in nLt.items()],
                    cellColours=[['white'] * 2 for key in nLt],
                    rowLabels=None,
                    colLabels=None,
                    bbox=(0.7, 0.53, 0.28, 0.45),
                    zorder=5,
                    )
            ax.table(
                    cellText=[[key, '{:1.2f}({:.0f})'.format(val, 100 * stdSt[key])] for key, val in nSt.items()],
                    cellColours=[['white'] * 2 for key in nSt],
                    rowLabels=None,
                    colLabels=None,
                    bbox=(0.02, 0.53 - 0.51 * cumulative, 0.28, 0.45),
                    zorder=5,
                    )

    def plot_single_correlations(
            self,
            co,
            figsize=(10, 8),
            yerr=None,
            fit_thresholds=True):
        from matplotlib.patches import Rectangle
        from matplotlib import cm
        import matplotlib.path as mpath
        import matplotlib.patches as mpatches
        Path = mpath.Path

        nL = co.nlargest(10)
        nS = co.nsmallest(10)

        fig, axs = plt.subplots(
                nrows=2, ncols=2,
                sharex=True, sharey=True,
                figsize=figsize)
        axs = axs.ravel()

        # Negative correlation, positive correlation, control, dynamic
        gids = [nS.index[0],
                nL.index[0],
                ds.featuresheet.loc[ds.featuresheet['GeneName'] == 'PSMB2'].index[0],
                ds.featuresheet.loc[ds.featuresheet['GeneName'] == 'COPE'].index[0],
                ]

        if fit_thresholds:
            mod = self.fit.fit_single(xs=['log_virus_reads_per_million'], ys=gids[:2], model='threshold-linear')
            thresholds = []
            for gid in gids[:2]:
                b, i, s, _ = mod.loc['log_virus_reads_per_million', gid].values
                t = (b - i) / s
                thresholds.append(10**t)
            thresholds.extend([None, None])
        else:
            thresholds = [1e4, 2e3, None, None]

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

            if gname in ('DDIT3', 'COPE'):
                loc = (0.98, 0.08)
                ha = 'right'
            else:
                loc = (0.02, 0.08)
                ha = 'left'
            if yerr is None:
                ax.text(*loc, '$\\rho = {:1.2f}$'.format(co.loc[gid]),
                        ha=ha,
                        transform=ax.transAxes)
            else:
                ax.text(*loc, '$\\rho = {:1.2f}({:.0f})$'.format(
                    co.loc[gid],
                    100 * yerr.loc[gid]),
                        ha=ha,
                        transform=ax.transAxes)

            th = thresholds[iax]
            if th is not None:
                ax.axvline(th, lw=2, color='darkred')

            # FIXME: this does not quite work yet
            #sns.kdeplot(
            #        x, y,
            #        zorder=3,
            #        n_levels=40,
            #        shade=True,
            #        ax=ax)

            ax.set_xlim(100, 3e5)
            ax.set_ylim(0.4, 1e5)
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.grid(True)
            ax.set_title(gname)

            # Add Bezier curve for some genes
            if gname == 'COPE':
                path = Path([(1.5e2, 1.2e3), (8e2, 6e2), (2e4, 3e3)],
                            [Path.MOVETO, Path.CURVE3, Path.CURVE3])
            elif gname == 'ACTB':
                path = Path([(thresholds[0], 6e4), (3e4, 6e4), (2e5, 1e3)],
                            [Path.MOVETO, Path.CURVE3, Path.CURVE3])
            elif gname == 'DDIT3':
                path = Path([(thresholds[1], 3e2), (3e3, 3e2), (3e4, 5e3)],
                            [Path.MOVETO, Path.CURVE3, Path.CURVE3])
            else:
                path = None

            if path is not None:
                pp1 = mpatches.FancyArrowPatch(
                    path=path,
                    arrowstyle='simple',
                    mutation_scale=20,
                    edgecolor='none',
                    facecolor='darkred',
                    zorder=10,
                    transform=ax.transData)
                ax.add_patch(pp1)

            # Add threshold-linear fits to first two
            if fit_thresholds:
                if gname in ['ACTB', 'DDIT3']:
                    b, i, s, _ = mod.loc['log_virus_reads_per_million', gid].values
                    xfit = np.linspace(-1, 5.5, 100)
                    yfit = i + s * xfit
                    t = (b - i) / s
                    yfit[xfit <= t] = b
                    xfit = 10**xfit - 0.1
                    yfit = 10**yfit - 0.1
                    ax.plot(
                            xfit,
                            yfit,
                            color='darkred',
                            lw=2,
                            ls='--',
                            alpha=0.5,
                            zorder=40)

        fig.text(0.5, 0.022,
                 'Dengue per million transcripts',
                 ha='center')
        fig.text(0.03, 0.5,
                 'Counts per million transcripts',
                 ha='center',
                 va='center',
                 rotation=90)

        plt.tight_layout(rect=(0.03, 0.03, 1, 1))

    def plot_correlations_candidates(self, co, genes, figsize=(10, 18)):
        fig, axs = plt.subplots(
                nrows=8, ncols=4,
                figsize=figsize)
        axs = axs.T.ravel()

        # Find IDs
        gids = [ds.featuresheet.loc[ds.featuresheet['GeneName'] == g].index[0]
                for g in genes]

        # Sort by correlation
        gids = co.loc[gids].sort_values().index

        thresholds = {
            'CTSC': 3e3,
            'CCND1': 5e3,
            'ID2': 4e3,
            'SEC61A1': 3e3,
            'RPL31': 3e1,
            'TRAM1': 3e3,
            'SEC11C': 4e3,
            'PLPP5': 5e3,
            'TMED2': 4e3,
            'SSR3': 4e3,
            'DDIT3': 3e3,
            'SELENOK': 2e3,
            'NUCB2': 2e3,
            'HSPA5': 2e3}

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

    print('Load dataset')
    ds = Fig2A_E(
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
    ds.samplesheet['log_virus_reads_per_million'] = np.log10(0.1 + ds.samplesheet['virus_reads_per_million'])

    print('Log counts')
    ds.counts.log(inplace=True)

    print('Get correlations')
    co = ds.correlation.correlate_features_phenotypes(
            phenotypes='virus_reads_per_million',
            fillna=0).fillna(0)

    print('Bootstrap top genes')
    top_pos = [
            'DDIT3', 'SELENOK', 'HERPUD1', 'CDK2AP2', 'CTH',
            'NUCB2', 'SDF2L1', 'SELENOS', 'HSPA5', 'DDIT4',
            ]
    top_neg = [
            'ACTB', 'TUBB', 'TUBB4B', 'ACTG1', 'HSPA8', 'LDHA',
            'SLC25A1', 'TUBA1B', 'HYAL2', 'PKM']
    examples = ['PSMB2', 'COPE']
    top_both = top_pos + top_neg + examples
    gids = [ds.featuresheet.index[ds.featuresheet['GeneName'] == gname][0] for gname in top_both]
    dboot = ds.copy()
    dboot.counts = ds.counts.loc[gids]
    n_bootstraps = 100
    co_boot = []
    for iboot in range(n_bootstraps):
        dbi = dboot.bootstrap()
        coi = dbi.correlation.correlate_features_phenotypes(
                phenotypes='virus_reads_per_million',
                fillna=0).fillna(0)
        coi.name = iboot
        co_boot.append(coi)
    co_boot = pd.concat(co_boot, axis=1)
    co_boot.columns.name = 'bootstrap#'
    co_mean = co_boot.mean(axis=1)
    co_std = co_boot.std(axis=1)

    print('Plot correlation distribution')
    ds.plot_correlation_distribution(
            co,
            figsize=(7, 4.6),
            yerr=co_std,
            #cumulative=False,
            cumulative=True,
            )

    #print('Plot single correlation plots')
    #ds.plot_single_correlations(co, figsize=(6, 4.6), yerr=co_std)

    #print('Plot correlations for supplementary')
    #ds.plot_correlations_candidates(co, genes_validation, figsize=(10, 18))

    plt.ion()
    plt.show()
