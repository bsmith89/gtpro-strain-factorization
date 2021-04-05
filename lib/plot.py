import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from scipy.spatial.distance import pdist, squareform
import numpy as np
from itertools import cycle
from sklearn.decomposition import PCA
from collections import defaultdict
import seaborn as sns
import scipy as sp
from lib.pandas import align_indexes


DEFAULT_MARKER_LIST = ["o", "v", "s", ">", "D", "X", "h", "^"]
DEFAULT_COLOR_LIST = ["black", "blue", "green", "orange", "purple", "red", "pink"]
DEFAULT_LINESTYLE_LIST = ["-", "--", "-.", ":"]


def construct_ordered_pallete(x, cm="Spectral", other="grey"):
    labels = pd.Series(x).unique()
    cm = mpl.cm.get_cmap(cm)
    colormap = defaultdict(lambda: other)
    for i, s in enumerate(labels):
        colormap[s] = cm(i / len(labels))
    return colormap


def demo_pallete(pallete):
    df = pd.DataFrame(pallete)
    print(df)


def pca_ordination(data, **kwargs):
    pca = PCA().fit(data)
    d1 = pca.transform(data)

    d1 = pd.DataFrame(
        d1, index=data.index, columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1]
    )
    frac_explained = pd.Series(pca.explained_variance_ratio_, index=d1.columns)
    return d1, frac_explained, {}


def nmds_ordination(data, **kwargs):
    # PCoA  # TODO: Modularize out?
    dmat = pd.DataFrame(
        squareform(pdist(data, **kwargs)), index=data.index, columns=data.index
    )
    init = MDS(
        n_components=2,
        max_iter=3000,
        eps=1e-9,
        random_state=1,
        dissimilarity="precomputed",
        n_jobs=1,
    ).fit_transform(dmat)
    nmds = MDS(
        n_components=2,
        metric=False,
        max_iter=3000,
        eps=1e-12,
        dissimilarity="precomputed",
        random_state=1,
        n_jobs=1,
        n_init=1,
    )
    d1 = nmds.fit_transform(dmat, init=init)

    d1 = pd.DataFrame(
        d1, index=data.index, columns=[f"PC{i}" for i in np.arange(d1.shape[1]) + 1]
    )
    frac_explained = pd.Series(np.nan, index=d1.columns)
    return d1, frac_explained, {"dmat": dmat}


def scatterplot(
    x,
    y,
    data,
    subset=None,
    colorby="__none__",
    color_palette=None,
    colorby_order=None,
    markerby="__none__",
    marker_palette=None,
    markerby_order=None,
    markersizeby="__none__",
    markersize_palette=None,
    markersizeby_order=None,
    edgecolorby="__none__",
    edgecolor_palette=None,
    edgecolorby_order=None,
    edgestyleby="__none__",
    edgestyle_palette=None,
    edgestyleby_order=None,
    ax=None,
    scatter_kws={},
):

    data = data.copy()

    # Set up plotting defaults
    if not ax:
        fig, ax = plt.subplots()

    # Add constant column "__none__":
    data["__none__"] = "__none__"

    if colorby_order is None:
        colorby_order = data[colorby].sort_values().unique()
    if color_palette is None:
        color_palette = {k: v for k, v in zip(colorby_order, cycle(DEFAULT_COLOR_LIST))}

    if markerby_order is None:
        markerby_order = data[markerby].sort_values().unique()
    if marker_palette is None:
        marker_palette = {
            k: v for k, v in zip(markerby_order, cycle(DEFAULT_MARKER_LIST))
        }

    if markersizeby_order is None:
        markersizeby_order = data[markersizeby].sort_values().unique()
    if markersize_palette is None:
        markersize_palette = {k: v for k, v in zip(markersizeby_order, cycle([None]))}

    if edgecolorby_order is None:
        edgecolorby_order = data[edgecolorby].sort_values().unique()
    if edgecolor_palette is None:
        edgecolor_palette = {
            k: v for k, v in zip(edgecolorby_order, cycle(DEFAULT_COLOR_LIST))
        }

    if edgestyleby_order is None:
        edgestyleby_order = data[edgestyleby].sort_values().unique()
    if edgestyle_palette is None:
        edgestyle_palette = {
            k: v for k, v in zip(edgestyleby_order, cycle(DEFAULT_LINESTYLE_LIST))
        }

    # Default scatter options, replaced if also in scatter_kws.
    scatter_kws_ = dict(lw=2)
    scatter_kws_.update(scatter_kws)

    for (
        feat_color,
        feat_marker,
        feat_markersize,
        feat_edgecolor,
        feat_edgestyle,
    ), d in data.groupby([colorby, markerby, markersizeby, edgecolorby, edgestyleby]):
        if (
            (feat_color not in colorby_order)
            or (feat_marker not in markerby_order)
            or (feat_edgecolor not in edgecolorby_order)
            or (feat_edgestyle not in edgestyleby_order)
            or (feat_markersize not in markersizeby_order)
        ):
            continue
        ax.scatter(
            x,
            y,
            data=d,
            marker=marker_palette[feat_marker],
            c=[color_palette[feat_color]],
            label="__nolegend__",
            edgecolor=[edgecolor_palette[feat_edgecolor]],
            linestyle=edgestyle_palette[feat_edgestyle],
            s=markersize_palette[feat_markersize],
            **scatter_kws_,
        )

    for feat_color in colorby_order:
        ax.scatter(
            [],
            [],
            marker="o",
            c=[color_palette[feat_color]],
            label=feat_color,
            **scatter_kws_,
        )
    for feat_marker in markerby_order:
        ax.scatter(
            [],
            [],
            marker=marker_palette[feat_marker],
            c="grey",
            label=feat_marker,
            **scatter_kws_,
        )
    for feat_markersize in markersizeby_order:
        ax.scatter(
            [],
            [],
            marker="o",
            s=markersize_palette[feat_markersize],
            c="grey",
            label=feat_markersize,
            **scatter_kws_,
        )
    for feat_edgecolor in edgecolorby_order:
        ax.scatter(
            [],
            [],
            marker="o",
            c="lightgrey",
            edgecolor=edgecolor_palette[feat_edgecolor],
            label=feat_edgecolor,
            **scatter_kws_,
        )
    for feat_edgestyle in edgestyleby_order:
        ax.scatter(
            [],
            [],
            marker="o",
            c="none",
            edgecolor="black",
            linestyle=edgestyle_palette[feat_edgestyle],
            label=feat_edgestyle,
            **scatter_kws_,
        )

    return ax


def ordination_plot(
    data,
    meta,
    subset=None,
    ordin=nmds_ordination,
    ordin_kws=dict(metric="cosine"),
    xy=("PC1", "PC2"),
    frac_explained_label=True,
    **kwargs,
):
    """Plot nMDS ordination with markers colored/shaped by metadata features.

    TODO: Sane/minimal defaults for kwarg=None.
    TODO: Align this function to be similar to Seaborn plotting functions.
    """
    x, y = xy
    data, meta, subset = align_indexes(data, meta, subset)
    data = data.loc[subset].copy()
    meta = meta.loc[subset].copy()

    d1, frac_explained, ordin_out = ordin(data, **ordin_kws)
    assert not d1.columns.isin(
        meta.columns
    ).any(), "Overlap between ordination dimension names and some metadata column."
    d1 = d1.join(meta)

    ax = scatterplot(x, y, d1, **kwargs)

    if frac_explained_label:
        ax.set_xlabel(f"{x} ({frac_explained[x]:0.1%})")
        ax.set_ylabel(f"{y} ({frac_explained[y]:0.1%})")
    else:
        ax.set_xlabel(f"{x}")
        ax.set_ylabel(f"{y}")
    ax.legend(bbox_to_anchor=(1, 1))
    return ax, ordin_out


def rotate_xticklabels(ax=None, rotation=45, ha="right", **kwargs):
    if ax is None:
        ax = plt.gca()
    ax.set_xticklabels(
        [x.get_text() for x in ax.get_xticklabels()], rotation=rotation, ha=ha, **kwargs
    )


def boxplot_with_points(
    x,
    y,
    hue=None,
    data=None,
    legend=True,
    ax=None,
    dist_plotter=sns.boxplot,
    dist_kwargs={},
    points_plotter=sns.stripplot,
    points_kwargs={},
    **kwargs,
):
    if not ax:
        fig, ax = plt.subplots()

    if data is None:
        data = {}
        data["x"] = x
        data["y"] = y
        x = "x"
        y = "y"
        if hue:
            data["hue"] = hue
            hue = "hue"
        data = pd.DataFrame(data)

    data_kwargs = dict(data=data, x=x, y=y, hue=hue, **kwargs)

    distkw = data_kwargs.copy()
    distkw.update(dict(showfliers=False, ax=ax, saturation=0.35))
    distkw.update(dist_kwargs)
    dist_plotter(**distkw)
    handles, labels = ax.get_legend_handles_labels()

    pointskw = data_kwargs.copy()
    pointskw.update(dict(dodge=True, linewidth=1, ax=ax))
    pointskw.update(points_kwargs)
    points_plotter(**pointskw)
    if legend:
        ax.legend(handles=handles)

    return ax


def residuals_plot(fit, ax=None, data=None, **kwargs):
    if not ax:
        fig, ax = plt.subplots()

    x = fit.predict()
    y = fit.resid

    # Add trend line, but x must be strictly increasing for UnivariateSpline
    # so here's a serious hack. :-/
    delta = lambda x: x[-1:] - x[:-1]
    min_delta = delta(np.unique(x)).min()
    x = x + np.random.uniform(-min_delta / 1000, min_delta / 1000, len(x))
    x, y = zip(*sorted(zip(x, y)))
    spline = sp.interpolate.UnivariateSpline(x, y)

    _data = fit.model.data.row_labels.to_series().to_frame()[[]].assign(x=x, y=y)
    if data is not None:
        _data = _data.join(data, rsuffix="_")

    ax.scatter("x", "y", data=_data, **kwargs)
    ax.plot(x, spline(x))
    ax.set_ylabel("Standardized Residuals")
    ax.set_xlabel("Fitted Value")
