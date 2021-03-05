import patsy
import scipy as sp
import numpy as np
from skbio.stats.composition import clr as _clr
import pandas as pd
from skbio.stats.distance import anosim as skb_anosim
from skbio.stats.distance import DistanceMatrix
from dataclasses import dataclass
from warnings import warn


def clr(rabund, pad=1e-10):
    return pd.DataFrame(_clr(rabund + pad), index=rabund.index, columns=rabund.columns)


def normalize(x, pseudo=0):
    x = x + pseudo
    return x / np.sum(x)


class RaiseLowTransform:
    """Replace unmeasured values with a fraction of the smallest measured value.

    -   Values are considered unmeasured if they're below *thresh*.
    -   Values are replaced with the smallest measured value divided by *factor*.

    """

    def __init__(self, thresh=0, factor=2):
        self.thresh = thresh
        self.minimum_measured = None
        self.factor = factor

    def memorize_chunk(self, x, thresh=0, factor=2):
        chunk_min = x[x > thresh].min()
        if self.minimum_measured is None:
            self.minimum_measured = chunk_min
        else:
            self.minimum_measured = min(self.minimum_measured, chunk_min)

    def memorize_finish(self, thresh=0, factor=2):
        assert self.minimum_measured > thresh
        self.replace_value = self.minimum_measured / factor

    def transform(self, x, thresh=0, factor=2):
        return np.where(x > thresh, x, np.full(x.shape, self.replace_value))


raise_low = patsy.stateful_transform(RaiseLowTransform)


def mannwhitneyu(x, y, data, reference=None, alternative="two-sided"):
    levels = data[x].unique()
    assert len(levels) == 2
    if reference is not None:
        assert reference in levels
        if not reference == levels[1]:
            levels = levels[1], reference
    x0 = data[data[x] == levels[0]][y]
    x1 = data[data[x] == levels[1]][y]
    try:
        u, pvalue = sp.stats.mannwhitneyu(x0, x1, alternative=alternative)
    except ValueError as err:
        warn(f"sp.stats.mannwhitneyu raised a ValueError: {err}")
        u, pvalue = (np.nan, np.nan)
    
    # common language effect size
    cles = u / (len(x0) * len(x1))
    return cles, pvalue


def wilcoxon(x, y, data):
    levels = data[x].unique()
    assert len(levels) == 2
    return sp.stats.wilcoxon(
        data[data[x] == levels[0]][y], data[data[x] == levels[1]][y]
    )


def anosim(dmat, groups, subset=None, n=999):
    if subset is not None:
        dmat = dmat.loc[subset, subset]
    groups = groups.loc[dmat.index]
    return skb_anosim(DistanceMatrix(dmat), groups, permutations=n)
