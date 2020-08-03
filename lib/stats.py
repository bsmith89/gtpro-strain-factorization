import patsy
import scipy as sp
import numpy as np

class RaiseLowTransform():
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

def mannwhitneyu(x, y, data):
    levels = data[x].unique()
    assert len(levels) == 2
    return sp.stats.mannwhitneyu(data[data[x] == levels[0]][y], data[data[x] == levels[1]][y])
