#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import theano
import pandas as pd
import numpy as np
from lib.util import idxwhere, normalize_rows


def list_samples(pileup):
    return list(pileup.columns.get_level_values(level='sample_id').unique())


def list_bases(pileup):
    return list(pileup.columns.get_level_values(level='base').unique())


def load_pileup_data(pileup_path):
    flat_data = pd.read_table(
                    pileup_path,
                    # TODO: Drop skiprows if input pileups files evolve
                    # to not include a header.
                    skiprows=1,
                    names=['segment_id', 'position',
                           'sample_id', 'base', 'tally'],
                    index_col=['segment_id', 'position',
                               'base', 'sample_id'],
                    squeeze=True
                             )
    pileup = (flat_data
              .unstack('base')
              .unstack('sample_id')
              .reorder_levels(['sample_id', 'base'], axis='columns')
              .fillna(0)
              .astype(theano.config.floatX)
              .sort_index(axis='columns'))
    # pileup = (flat_data.unstack(['sample_id', 'base'], fill_value=0)
    #                    .astype(theano.config.floatX)
    #                    .sort_index(axis='columns'))
    sample_list = list_samples(pileup)
    base_list = list_bases(pileup)
    assert pileup.shape[1] == len(sample_list) * len(base_list)
    return pileup


def filter_sites(pileup, min_max_counts,
                 max_major_allele_freq,
                 max_z=2, min_z=-3,
                 max_third_allele_frac=0.05,
                 subset=None):
    if subset is None:
        subset = list_samples(pileup)

    # Max counts thresholding.
    pileup_ = pileup.loc[:, subset]
    total_counts = pileup_.groupby(level='sample_id', axis='columns').sum()
    max_counts = total_counts.max(1)
    pileup_ = pileup_[(max_counts >= min_max_counts)]

    base_coverage = pileup_.groupby(level='base', axis='columns').sum()
    ranked_allele_freq = (normalize_rows(base_coverage)
                          .apply(lambda x:
                                 pd.Series(x.sort_values(ascending=False)
                                            .values,
                                           # FIXME: Why is it range(1, 5)?
                                           index=range(1, 5)), axis=1))
    pileup_ = pileup_[(ranked_allele_freq[1] <= max_major_allele_freq)
                      & (ranked_allele_freq[3] <= max_third_allele_frac)]

    log_counts_p1 = np.log(total_counts.loc[pileup_.index] + 1)
    log_counts_z = ((log_counts_p1 - log_counts_p1.mean())
                    / log_counts_p1.std())
    pileup_ = pileup_[(log_counts_z.max(1) < max_z) &
                      (log_counts_z.min(1) > min_z)]

    return pileup.loc[pileup_.index]


def filter_samples(pileup, min_median_coverage):
    median_coverage = (pileup.groupby(level='sample_id', axis='columns')
                             .sum()
                             .median())
    return pileup.loc[:, idxwhere(median_coverage >= min_median_coverage)]


def ranked_alleles(pileup):
    a = len(list_bases(pileup))
    base_coverage = pileup.groupby(level='base', axis='columns').sum()

    rank_df_columns = ([f'f{i}' for i in range(1, a + 1)] +
                       [f'b{i}' for i in range(1, a + 1)])

    def ranked(x):
        x = x.sort_values(ascending=False)
        values = list(x.values)
        idx = list(x.index)
        return pd.Series(values + idx, index=rank_df_columns)

    return normalize_rows(base_coverage).apply(ranked, axis=1)


def convert_to_major_minor_allele(pileup):
    classified_base = ranked_alleles(pileup).rename({'b1': 'major',
                                                     'b2': 'minor'},
                                                    axis='columns')
    major = (pileup.stack('base')
             .loc[classified_base.set_index('major', append=True).index]
             .reset_index(level='major', drop=True)
             .stack()
             .rename('major'))
    minor = (pileup.stack('base')
             .loc[classified_base.set_index('minor', append=True).index]
             .reset_index(level='minor', drop=True)
             .stack()
             .rename('minor'))
    pileup = (pd.concat([major, minor],
                        axis='columns')
              .rename_axis(columns='base')
              .stack().unstack(['sample_id', 'base']))
    return pileup, classified_base[['major', 'minor']]


def sample_sites(pileup, n, random_state=0):
    index = pileup.index
    sample = pileup.sample(n, random_state=random_state)
    sample_index = index[index.isin(sample.index)]
    return sample.loc[sample_index]


def get_pileup_dims(pileup):
    n = len(list_samples(pileup))
    a = len(list_bases(pileup))
    g = pileup.shape[0]
    return g, n, a
