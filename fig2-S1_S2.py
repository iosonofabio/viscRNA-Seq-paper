# vim: fdm=indent
'''
author:     Fabio Zanini
date:       08/09/17
content:    Test parametric fitting of thresholds.
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

    print('Log counts')
    ds.counts.log(inplace=True)

    print('Get correlations')
    co = ds.correlation.correlate_features_phenotypes(
            phenotypes='virus_reads_per_million',
            fillna=0).fillna(0)

    print('Focus on decently correlated genes')
    top_pos = [
            'DDIT3', 'SELENOK', 'HERPUD1', 'CDK2AP2', 'CTH',
            'NUCB2', 'SDF2L1', 'SELENOS', 'HSPA5', 'DDIT4',
            ]
    top_neg = [
            'ACTB', 'TUBB', 'TUBB4B', 'ACTG1', 'HSPA8',
            'LDHA', 'SLC25A1', 'TUBA1B', 'HYAL2', 'PKM']
    examples = ['PSMB2', 'COPE']
    top_both = top_pos + top_neg + examples
    co_min = 0.3
    gids = co.index[np.abs(co) >= co_min]
    #gids = [ds.featuresheet.index[ds.featuresheet['GeneName'] == gname][0] for gname in top_both]
    dsn = ds.copy()
    dsn.counts = dsn.counts.loc[gids]
    dsn.rename(axis='features', column='GeneName', inplace=True)

    # Fit threshold-linear model
    mod = dsn.fit.fit_single(xs=['log_virus_reads_per_million'], ys='mapped', model='threshold-linear')

    # Plot a few (Fig. S8)
    gnames = top_pos + top_neg
    fig, axs = plt.subplots(
            4, 5, figsize=(14, 10),
            sharex=True,
            sharey=True,
            squeeze=False)
    axs = axs.ravel()
    for iax, (gname, ax) in enumerate(zip(gnames, axs)):
        x = dsn.samplesheet['log_virus_reads_per_million']
        y = dsn.counts.loc[gname]
        ax.scatter(x, y, alpha=0.5, label='')
        ax.set_xlim(-1.1, 5.5)
        ax.set_ylim(-1.1, 5.5)
        ax.grid(True)
        ax.set_title(gname)

        # Plot fit
        b, i, s, _ = mod.loc['log_virus_reads_per_million', gname].values
        xfit = np.linspace(-1, 5.5, 100)
        yfit = i + s * xfit
        t = (b - i) / s
        yfit[xfit <= t] = b
        ax.plot(xfit, yfit, color='darkred', lw=2, label='T=$10^{{{:.1f}}}$'.format(t))

        if iax < 10:
            leg_loc = 'lower right'
        else:
            leg_loc = 'lower left'
        ax.legend(loc=leg_loc)

    ax.set_yticks([0, 2, 4])
    ax.set_yticklabels(['$1$', '$10^2$', '$10^4$'])
    ax.set_xticks([-1, 0, 1, 2, 3, 4, 5])
    ax.set_xticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])

    fig.text(0.038, 0.75, 'Correlated genes', rotation=90, va='center', ha='right')
    fig.text(0.038, 0.25, 'Anticorrelated genes', rotation=90, va='center', ha='right')
    fig.text(0.017, 0.5,
             'Gene counts per million transcripts',
             rotation=90, va='center', ha='right')
    fig.text(0.5, 0.03,
             'Virus reads per million transcripts',
             ha='center', va='top')

    plt.tight_layout(rect=(0.04, 0.02, 1, 1), h_pad=0.5, w_pad=0)

    # Plot correlation versus threshold (Fig. S9)
    fig, axs = plt.subplots(1, 2, figsize=(9, 5), gridspec_kw={'width_ratios': [5, 1]}, sharey=True)
    ax = axs[0]
    tran = ds.featuresheet.loc[gids, 'GeneName']
    x = co.loc[gids].values
    mm = mod.loc['log_virus_reads_per_million', tran.values]
    y = ((mm.loc[:, 'baseline'] - mm.loc[:, 'intercept']) / mm.loc[:, 'slope']).values

    xmax = np.abs(x).max()
    alphas = 0.05 + 0.85 * (np.abs(x) - co_min) / (xmax - co_min)
    rgba_colors = np.zeros((len(x), 4))
    rgba_colors[:, :3] = [0.27, 0.51, 0.71]
    rgba_colors[:, 3] = alphas
    ax.scatter(x, y, color=rgba_colors, zorder=6)
    sns.kdeplot(x, y, cmap='viridis', zorder=5, ax=ax)

    ax.set_xlabel('Correlation with intracellular amount of virus')
    ax.set_ylabel('Minimal virus amount\nfor gene expression change')
    ax.set_xlim(-0.75, 0.75)
    ax.set_ylim(-1.3, 5.7)
    ax.set_yticks(np.arange(-1, 6))
    ax.set_yticklabels(['$0$', '$1$', '$10$', '$10^2$', '$10^3$', '$10^4$', '$10^5$'])
    ax.grid(True)

    # Stratify distributions
    ax = axs[1]
    #bins = [(-1, -0.4), (-0.4, -0.3), (0.3, 0.4), (0.4, 1)]
    bins = [(-1, -0.3), (0.3, 1)]
    colors = ['steelblue', 'darkred']
    ax.set_xlim(-0.5, len(bins) + 0)
    tmp = []
    for ib, bini in enumerate(bins):
        indi = (x >= bini[0]) & (x <= bini[1])
        xi = x[indi]
        yi = y[indi]
        cati = [str(bini) for i in range(sum(indi))]
        tmp.append(pd.DataFrame(data=[xi, yi, cati], index=['x', 'y', 'cat']).T)

        ax.hist(yi, bottom=ib, orientation='horizontal', density=True, color=colors[ib], alpha=0.5)
        ax.axvline(ib, lw=1.5, color='k')
        ax.axhline(y=np.median(yi), lw=1.5, alpha=1, color=colors[ib])

    tmp = pd.concat(tmp)

    ax.set_xticks(np.arange(len(bins)))
    if len(bins) != 2:
        ax.set_xticklabels([str(bini) for bini in bins], rotation=90)
    else:
        ax.set_xticklabels(['anticorrelated', 'correlated'], rotation=30)
    ax.set_ylim(-1.3, 5.7)
    ax.grid(True)

    plt.tight_layout()

    # Identify the order of events
    indp = x > 0.4
    tranp = tran.loc[indp]
    xp = x[indp]
    yp = y[indp]
    sp = mm.loc[:, 'slope'].values[indp]
    inds = np.argsort(yp)
    print('{:10s}\t{:}\t{:}\t{:}'.format('Gene', 'Th', 'rho', 'slope'))
    for i in inds:
        print('{:10s}\t{:.2f}\t{:.2f}\t{:.2f}'.format(tranp.iloc[i], yp[i], xp[i], sp[i]))

    plt.ion()
    plt.show()
