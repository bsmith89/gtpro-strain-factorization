#!/usr/bin/env python3

import xarray as xr
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
from lib.util import info, idxwhere
import sys

DEFAULT_COMPRESSION_LEVEL = 6


def xr_idxwhere(x):
    return idxwhere(x.to_series())


def select_informative_positions(data, incid_thresh):
    minor_allele_incid = (data > 0).mean("library_id").min("allele")
    informative_positions = xr_idxwhere(minor_allele_incid > incid_thresh)
    return informative_positions


def write_netcdf(path, data, compression_level=DEFAULT_COMPRESSION_LEVEL):
    (
        data.to_dataset(name="tally").to_netcdf(
            path,
            encoding=dict(tally=dict(zlib=True, complevel=compression_level)),
        )
    )


if __name__ == "__main__":
    inpath = sys.argv[1]  # 'data/ucfmt.sp-102506.gtpro-pileup.nc'
    minor_allele_thresh = float(sys.argv[2])  # 0.05
    position_thresh = float(sys.argv[3])  # 0.25
    npos_subsample = int(sys.argv[4])
    dist_thresh = float(sys.argv[5])  # 0.2
    clust_size_thresh = int(sys.argv[6])  # 3
    clust_pos_frac_thresh = float(sys.argv[7])  # 0.5
    pos_incid_thresh = float(sys.argv[8])  # 1.0
    outpath = sys.argv[9]  # 'data/ucfmt.sp-102506.gtpro-pileup.filt.nc'

    info(f"Loading data from `{inpath}`.")
    data = xr.open_dataset(inpath).tally.sum('read').squeeze()

    info(
        f"Identifying positions with "
        f">{minor_allele_thresh:.1%} minor allele incidence."
    )
    inf_positions = select_informative_positions(data, minor_allele_thresh)
    npos = len(inf_positions)
    info(f"Found {npos} informative positions.")

    info(
        f"Identifying libraries with observations at "
        f">{position_thresh:.1%} of informative positions."
    )
    frac_positions = (
        (data.sel(position=inf_positions) > 0).any("allele").mean("position")
    )
    high_cvrg_libs = xr_idxwhere(frac_positions > position_thresh)
    nlibs = len(high_cvrg_libs)
    info(f"Found {nlibs} high-coverage libraries.")

    npos_subsample = min(npos_subsample, len(inf_positions))
    info(f"Subsample {npos_subsample} from informative positions.")
    position_ss = np.random.choice(
        inf_positions,
        size=min(npos_subsample, len(inf_positions)),
        replace=False,
    )

    info(f"Transforming counts into allele scores.")
    trsnf_data = data.sel(library_id=high_cvrg_libs, position=position_ss)
    cvrg = trsnf_data.sum("allele")
    trsnf_data = (trsnf_data + 1) / (cvrg + 2)
    trsnf_data = trsnf_data.sel(allele="alt") * 2 - 1

    info(
        f"Clustering {nlibs} metagenotypes at a {dist_thresh} "
        f"maximum cosine distance threshold across {npos} positions."
    )
    clust = AgglomerativeClustering(
        n_clusters=None,
        affinity="cosine",
        linkage="complete",
        distance_threshold=dist_thresh,
    ).fit(trsnf_data)
    clust = pd.Series(clust.labels_, index=trsnf_data.library_id)

    clust_size = clust.value_counts()
    large_clusts = idxwhere(clust_size >= clust_size_thresh)
    nclusts = len(large_clusts)
    total_clust_libs = clust_size.loc[large_clusts].sum()
    info(
        f"Found {nclusts} ('large') clusters with at least "
        f"{clust_size_thresh} members, encompassing a total of "
        f"{total_clust_libs} libraries."
    )

    info(
        f"Constructing cluster-by-position matrix indicating "
        f"which large clusters have counts at each position in "
        f">{clust_pos_frac_thresh:.1%} of samples."
    )
    cvrg_subset = cvrg.to_pandas().loc[clust.index]
    clust_pos_frac = {}
    for c in large_clusts:
        clust_pos_frac[c] = (cvrg_subset.loc[clust == c] > 0).mean()
    clust_pos_frac = pd.DataFrame(clust_pos_frac, index=inf_positions)
    clust_pos_incid = (clust_pos_frac >= clust_pos_frac_thresh).astype(int)
    info(f"Calculating weighted cluster incidence.")
    pos_incid = (clust_pos_incid * clust_size.loc[large_clusts]).sum(
        1
    ) / clust_size.loc[large_clusts].sum()
    info(
        f"Identifying positions with weighted incidence "
        f"of {pos_incid_thresh:.1%}."
    )
    pos_passing = idxwhere(pos_incid >= pos_incid_thresh)
    npos_pass = len(pos_passing)
    info(f"Total positions passing all filters: {npos_pass}.")

    info(f"Writing output to `{outpath}`.")
    write_netcdf(outpath, data.sel(position=pos_passing))
