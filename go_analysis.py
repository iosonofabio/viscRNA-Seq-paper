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
class GOAnalysis(Dataset):
    pass


# Script
if __name__ == '__main__':

    ds = GOAnalysis(
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

    #print('Plot correlation distribution')
    #ds.plot_correlation_distribution(co, figsize=(6.5, 4.6))

    #print('Plot single correlation plots')
    #ds.plot_single_correlations(co, figsize=(6.5, 4.6))
    ##ds.plot_correlations_candidates(co, genes_validation, figsize=(18, 12))

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
        co.name = t
        cotimes.append(co)
    cotimes = pd.concat(cotimes, axis=1)
    cotimes.columns.name = 'time'


