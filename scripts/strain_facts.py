#!/usr/bin/env python3

import sys
import argparse
import warnings

import matplotlib.pyplot as plt
from lib.util import info, idxwhere
import numpy as np
import pyro
import pyro.distributions as dist
import torch
from functools import partial
from tqdm import tqdm
import xarray as xr
import pandas as pd

from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering


def as_torch(x, dtype=None, device=None):
    # Cast inputs and set device
    return torch.tensor(x, dtype=dtype, device=device)


def as_torch_all(dtype=None, device=None, **kwargs):
    # Cast inputs and set device
    return {k: as_torch(kwargs[k], dtype=dtype, device=device) for k in kwargs}


def rss(x, y):
    return np.sqrt(np.sum((x - y) ** 2))


def binary_entropy(p):
    q = 1 - p
    return -p * np.log2(p) - q * np.log2(q)


def genotype_distance(x, y):
    x = x * 2 - 1
    y = y * 2 - 1
    weight = (x * y) ** 2
    dist = ((x - y) / 2) ** 2
    wmean_dist = (weight * dist).sum() / (weight.sum())
    return wmean_dist


def plot_loss_history(loss_history):
    min_loss = loss_history.min()
    plt.plot(loss_history - min_loss)
    plt.plot(
        np.linspace(0, len(loss_history), num=1000),
        np.linspace(len(loss_history), 0, num=1000),
        lw=1,
        linestyle="--",
        color="grey",
    )
    plt.title(f"+{min_loss:0.3e}")
    plt.yscale("log")
    return plt.gca()


def mean_residual_count(expect_frac, obs_count, m):
    frac_obs = obs_count / m
    out = np.abs(((frac_obs - expect_frac)))
    out[np.isnan(out)] = 0
    return (out * m).sum() / m.sum()


def _load_input_data(allpaths):
    data = []
    for path in allpaths:
        info(path)
        d = xr.open_dataarray(path).squeeze()
        info(f"Shape: {d.sizes}.")
        data.append(d)
    info("Concatenating data from {} files.".format(len(data)))
    data = xr.concat(data, "library_id", fill_value=0)
    return data


def select_informative_positions(data, incid_thresh):
    minor_allele_incid = (data > 0).mean("library_id").min("allele")
    informative_positions = idxwhere(
        minor_allele_incid.to_series() > incid_thresh
    )
    return informative_positions


def model(
    s,
    m,
    y=None,
    gamma_hyper=1.0,
    pi0=1.0,
    rho0=1.0,
    epsilon0=0.01,
    alpha0=1000.0,
    dtype=torch.float32,
    device="cpu",
):

    # Cast inputs and set device
    m, gamma_hyper, pi0, rho0, epsilon0, alpha0 = [
        torch.tensor(v, dtype=dtype, device=device)
        for v in [m, gamma_hyper, pi0, rho0, epsilon0, alpha0]
    ]
    if y is not None:
        y = torch.tensor(y)

    n, g = m.shape

    with pyro.plate("position", g, dim=-1):
        with pyro.plate("strain", s, dim=-2):
            gamma = pyro.sample("gamma", dist.Beta(gamma_hyper, gamma_hyper))
    # gamma.shape == (s, g)

    rho_hyper = pyro.sample("rho_hyper", dist.Gamma(rho0, 1.0))
    rho = pyro.sample(
        "rho",
        dist.Dirichlet(torch.ones(s, dtype=dtype, device=device) * rho_hyper),
    )

    epsilon_hyper = pyro.sample("epsilon_hyper", dist.Beta(1.0, 1 / epsilon0))
    alpha_hyper = pyro.sample("alpha_hyper", dist.Gamma(alpha0, 1.0))
    pi_hyper = pyro.sample("pi_hyper", dist.Gamma(pi0, 1.0))

    with pyro.plate("sample", n, dim=-1):
        pi = pyro.sample("pi", dist.Dirichlet(rho * s * pi_hyper))
        alpha = pyro.sample("alpha", dist.Gamma(alpha_hyper, 1.0)).unsqueeze(
            -1
        )
        epsilon = pyro.sample(
            "epsilon", dist.Beta(1.0, 1 / epsilon_hyper)
        ).unsqueeze(-1)
    # pi.shape == (n, s)
    # alpha.shape == epsilon.shape == (n,)

    p_noerr = pyro.deterministic("p_noerr", pi @ gamma)
    p = pyro.deterministic(
        "p", (1 - epsilon / 2) * (p_noerr) + (epsilon / 2) * (1 - p_noerr)
    )
    # p.shape == (n, g)

    y = pyro.sample(
        "y",
        dist.BetaBinomial(
            concentration1=alpha * p,
            concentration0=alpha * (1 - p),
            total_count=m,
        ),
        obs=y,
    )
    # y.shape == (n, g)
    return y


def conditioned_model(
    model,
    data={},
    dtype=torch.float32,
    device="cpu",
    **kwargs,
):
    data = {
        k: torch.tensor(v, dtype=dtype, device=device) for k, v in data.items()
    }
    return partial(
        pyro.condition(model, data=data),
        dtype=dtype,
        device=device,
        **kwargs,
    )


def find_map(
    model,
    init={},
    lag=10,
    stop_at=1.0,
    max_iter=int(1e5),
    learning_rate=1e-0,
    clip_norm=100.0,
):
    guide = pyro.infer.autoguide.AutoLaplaceApproximation(
        model,
        init_loc_fn=pyro.infer.autoguide.initialization.init_to_value(init),
    )

    svi = pyro.infer.SVI(
        model,
        guide,
        pyro.optim.Adamax(
            optim_args={"lr": learning_rate},
            clip_args={"clip_norm": clip_norm},
        ),
        loss=pyro.infer.JitTrace_ELBO(),
    )

    pyro.clear_param_store()
    pbar = tqdm(range(max_iter), position=0, leave=True)
    history = []
    try:
        for i in pbar:
            elbo = svi.step()

            if np.isnan(elbo):
                break

            # Fit tracking
            history.append(elbo)

            # Reporting/Breaking
            if i < 2:
                pbar.set_postfix({"ELBO": history[-1]})
            elif i < lag + 1:
                pbar.set_postfix(
                    {
                        "ELBO": history[-1],
                        "delta_1": history[-2] - history[-1],
                    }
                )
            else:
                delta_lag = (history[-lag] - history[-1]) / lag
                pbar.set_postfix(
                    {
                        "ELBO": history[-1],
                        "delta_1": history[-2] - history[-1],
                        f"delta_{lag}": delta_lag,
                    }
                )
                if delta_lag < stop_at:
                    info("Optimization converged")
                    break
    except KeyboardInterrupt:
        info("Optimization interrupted")
    pbar.refresh()
    assert delta_lag < stop_at, (
        f"Reached {args.max_iter} iterations with a per-step improvement of "
        f"{args.delta_lag}. Consider setting --max-iter "
        f"or --stop-at larger; increasing --learning-rate may also help, "
        f"although it could also lead to numerical issues."
    )
    # Gather MAP from parameter-store
    mapest = {
        k: v.detach().cpu().numpy().squeeze()
        for k, v in pyro.infer.Predictive(
            model, guide=guide, num_samples=1
        )().items()
    }
    return mapest, np.array(history)


def initialize_parameters_by_clustering_samples(
    y, m, thresh, addition_strains_factor=0.5
):
    n, g = y.shape

    sample_metagenotype = y + 1 / m + 2
    geno_dist = squareform(
        pdist(sample_metagenotype, metric=genotype_distance)
    )

    clust = pd.Series(
        AgglomerativeClustering(
            n_clusters=None,
            affinity="precomputed",
            linkage="complete",
            distance_threshold=args.collapse,
        )
        .fit(geno_dist)
        .labels_
    )

    y_total = (
        pd.DataFrame(pd.DataFrame(y))
        .groupby(clust, axis="columns")
        .sum()
        .values
    )
    m_total = (
        pd.DataFrame(pd.DataFrame(m))
        .groupby(clust, axis="columns")
        .sum()
        .values
    )
    clust_genotype = (y_total + 1) / (m_total + 2)
    additional_haplotypes = addition_strains_factor * clust_genotype.shape[0]

    init_genotype = pd.concat(
        [
            clust_genotype,
            pd.DataFrame(np.ones((additional_haplotypes, g)) * 0.5),
        ]
    ).reset_index(drop=True)

    s_init = init_genotype.shape[0]
    init_frac = np.ones((n, s_init))
    for i in range(n):
        init_frac[i, clust[i]] = s_init - 1
    init_frac /= init_frac.sum(1, keepdims=True)

    return init_genotype, init_frac


def parse_args(argv):
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Input
    p.add_argument(
        "pileup",
        nargs="+",
        help="""
Single, fully processed, pileup table in NetCDF format
with the following dimensions:
    * library_id
    * position
    * allele
                        """,
    )

    # Shape of the model
    p.add_argument("--nstrains", metavar="INT", type=int, default=1000)
    p.add_argument(
        "--npos",
        metavar="INT",
        default=2000,
        type=int,
        help=("Number of positions to sample for model fitting."),
    )

    # Data filtering
    p.add_argument(
        "--incid-thresh",
        metavar="FLOAT",
        type=float,
        default=0.02,
        help=(
            "Minimum fraction of samples that must have the minor allele "
            "for the position to be considered 'informative'."
        ),
    )
    p.add_argument(
        "--cvrg-thresh",
        metavar="FLOAT",
        type=float,
        default=0.5,
        help=(
            "Minimum fraction of 'informative' positions with counts "
            "necessary for sample to be included."
        ),
    )

    # Regularization
    p.add_argument(
        "--gamma-hyper",
        metavar="FLOAT",
        default=1e-2,
        type=float,
        help=("Ambiguity regularization parameter."),
    )
    p.add_argument(
        "--pi-hyper",
        metavar="FLOAT",
        default=1e-1,
        type=float,
        help=(
            "Heterogeneity regularization parameter (will be scaled by 1 / s)."
        ),
    )
    p.add_argument(
        "--rho-hyper",
        metavar="FLOAT",
        default=1e0,
        type=float,
        help=("Diversity regularization parameter."),
    )
    p.add_argument(
        "--epsilon-hyper", metavar="FLOAT", default=0.01, type=float
    )
    p.add_argument(
        "--epsilon",
        metavar="FLOAT",
        default=None,
        type=float,
        help=("Fixed error rate for all samples."),
    )
    p.add_argument("--alpha-hyper", metavar="FLOAT", default=100.0, type=float)
    p.add_argument(
        "--alpha",
        metavar="FLOAT",
        default=None,
        type=float,
        help=("Fixed concentration for all samples."),
    )
    p.add_argument(
        "--collapse",
        metavar="FLOAT",
        default=0.0,
        type=float,
        help=("Merge strains with a cosine distance of less than this value."),
    )

    # Fitting
    p.add_argument("--random-seed", default=0, type=int, help=("FIXME"))
    p.add_argument("--max-iter", default=10000, type=int, help=("FIXME"))
    p.add_argument("--lag", default=50, type=int, help=("FIXME"))
    p.add_argument("--stop-at", default=5.0, type=float, help=("FIXME"))
    p.add_argument("--learning-rate", default=1e-0, type=float, help=("FIXME"))
    p.add_argument("--clip-norm", default=100.0, type=float, help=("FIXME"))

    # Hardware
    p.add_argument("--device", default="cpu", help=("PyTorch device name."))

    # Output
    p.add_argument(
        "--outpath",
        metavar="PATH",
        help=("Path for genotype fraction output."),
    )

    args = p.parse_args(argv)

    return args


if __name__ == "__main__":
    warnings.filterwarnings(
        "error",
        message="Encountered NaN: loss",
        category=UserWarning,
        # module="trace_elbo",  # FIXME: What is the correct regex for module?
        lineno=217,
    )
    warnings.filterwarnings(
        "ignore",
        message="CUDA initialization: Found no NVIDIA",
        category=UserWarning,
        lineno=130,
    )
    warnings.filterwarnings(
        "ignore",
        message="torch.tensor results are registered as constants",
        category=torch.jit.TracerWarning,
        # module="trace_elbo",  # FIXME: What is the correct regex for module?
        lineno=95,
    )

    args = parse_args(sys.argv[1:])
    info(args)

    info(f"Setting random seed: {args.random_seed}")
    np.random.seed(args.random_seed)

    info("Loading input data.")
    data = _load_input_data(args.pileup)
    info(f"Full data shape: {data.sizes}.")

    info("Filtering positions.")
    informative_positions = select_informative_positions(
        data, args.incid_thresh
    )
    npos_available = len(informative_positions)
    info(
        f"Found {npos_available} informative positions with minor "
        f"allele incidence of >{args.incid_thresh}"
    )
    npos = min(args.npos, npos_available)
    info(f"Randomly sampling {npos} positions.")
    position_ss = np.random.choice(
        informative_positions,
        size=npos,
        replace=False,
    )

    info("Filtering libraries.")
    suff_cvrg_samples = idxwhere(
        (
            (
                data.sel(position=informative_positions).sum(["allele"]) > 0
            ).mean("position")
            > args.cvrg_thresh
        ).to_series()
    )
    nlibs = len(suff_cvrg_samples)
    info(
        f"Found {nlibs} libraries with >{args.cvrg_thresh:0.1%} "
        f"of informative positions covered."
    )

    info("Constructing input data.")
    data_fit = data.sel(library_id=suff_cvrg_samples, position=position_ss)
    m_ss = data_fit.sum("allele")
    n, g_ss = m_ss.shape
    y_obs_ss = data_fit.sel(allele="alt")

    info("Initializing parameter estimates based on clustering analysis.")
    info("TODO: Info about clustering parameters.")
    init_genotype, init_frac = initialize_parameters_by_clustering_samples(
        y_obs_ss,
        m_ss,
        thresh=thresh,
        additional_strains_factor=args.additional_nstrains,
    )

    s = init_genotype.shape[1]
    info(f"Model shape: n={n}, g={g_ss}, s={s}")

    info(
        "Fitting intial strain fractions and genotypes with ambiguity "
        "regularization."
    )
    conditioning_data = dict(
        alpha_hyper=args.alpha_hyper,
        epsilon_hyper=args.epsilon_hyper,
        pi_hyper=args.pi_hyper / s,
        rho_hyper=args.rho_hyper,
        y=y_obs_ss.values,
    )
    if args.alpha:
        conditioning_data["alpha"] = np.ones(n) * args.alpha
    if args.epsilon:
        conditioning_data["epsilon"] = np.ones(n) * args.epsilon
    info(f"Conditioning model on: {conditioning_data.keys()}")

    model_fit = conditioned_model(
        model,
        data=conditioning_data,
        s=s,
        m=m_ss.values,
        gamma_hyper=args.gamma_hyper,
        dtype=torch.float32,
        device=args.device,
    )

    info("Optimizing model parameters.")
    mapest1, history1 = find_map(
        model_fit,
        init=as_torch_all(
            gamma=init_genotype,
            pi=init_frac,
            dtype=torch.float32,
            device=args.device,
        ),
        lag=args.lag,
        stop_at=args.stop_at,
        learning_rate=args.learning_rate,
        max_iter=args.max_iter,
        clip_norm=args.clip_norm,
    )
    if args.device.startswith("cuda"):
        info(
            "CUDA available mem: {}".format(
                torch.cuda.get_device_properties(0).total_memory
            ),
        )
        info("CUDA reserved mem: {}".format(torch.cuda.memory_reserved(0)))
        info("CUDA allocated mem: {}".format(torch.cuda.memory_allocated(0)))
        info(
            "CUDA free mem: {}".format(
                torch.cuda.memory_reserved(0) - torch.cuda.memory_allocated(0)
            )
        )
        torch.cuda.empty_cache()

    info(
        "Re-fitting genotypes, this time without ambiguity regularization, "
        "conditioned on strain fractions from Round 1."
    )
    conditioning_data = dict(
        alpha_hyper=args.alpha_hyper,
        epsilon_hyper=args.epsilon_hyper,
        pi=mapest1["pi"],
        y=y_obs_ss.values,
    )
    if args.alpha:
        conditioning_data["alpha"] = np.ones(n) * args.alpha
    if args.epsilon:
        conditioning_data["epsilon"] = np.ones(n) * args.epsilon
    info(f"Conditioning model on: {conditioning_data.keys()}")

    model_fit = conditioned_model(
        model,
        data=conditioning_data,
        s=s,
        m=m_ss.values,
        gamma_hyper=1e0,
        dtype=torch.float32,
        device=args.device,
    )

    info("Optimizing model parameters.")
    mapest2, history2 = find_map(
        model_fit,
        lag=args.lag,
        stop_at=args.stop_at,
        learning_rate=args.learning_rate,
        max_iter=args.max_iter,
        clip_norm=args.clip_norm,
    )
    if args.device.startswith("cuda"):
        torch.cuda.empty_cache()

    if args.collapse > 0:
        # Collapse genotypes
        info(
            f"Collapsing genotype clusters with maximum-distances "
            f"<{args.collapse}."
        )
        info(f"See the definition of the function: 'genotype_distance'.")
        info(
            "Calculating strain distances. "
            "This can take a while for large numbers of strains."
        )
        geno_dist = squareform(
            pdist(mapest2["gamma"], metric=genotype_distance)
        )

        info("Clustering.")
        clust = pd.Series(
            AgglomerativeClustering(
                n_clusters=None,
                affinity="precomputed",
                linkage="complete",
                distance_threshold=args.collapse,
            )
            .fit(geno_dist)
            .labels_
        )
        pi_collapse = (
            pd.DataFrame(mapest2["pi"])
            .groupby(clust, axis="columns")
            .sum()
            .values
        )
        s_collapse = pi_collapse.shape[1]
        info(f"Collapsed {s} genotypes into {s_collapse}.")

        # Re-fit genotypes based on collapsed strain abundances in order
        # to fully integrate abundances details into genotypes.
        info("Fitting genotypes based on collapsed strain fractions.")
        conditioning_data = dict(
            alpha_hyper=args.alpha_hyper,
            epsilon_hyper=args.epsilon_hyper,
            pi=pi_collapse,
            y=y_obs_ss.values,
        )
        if args.alpha:
            conditioning_data["alpha"] = np.ones(n) * args.alpha
        if args.epsilon:
            conditioning_data["epsilon"] = np.ones(n) * args.epsilon
        info(f"Conditioning model on: {conditioning_data.keys()}")

        model_fit = conditioned_model(
            model,
            data=conditioning_data,
            s=s_collapse,
            m=m_ss.values,
            gamma_hyper=1e0,
            dtype=torch.float32,
            device=args.device,
        )

        info("Optimizing model parameters.")
        mapest3, history3 = find_map(
            model_fit,
            lag=args.lag,
            stop_at=args.stop_at,
            learning_rate=args.learning_rate,
            max_iter=args.max_iter,
            clip_norm=args.clip_norm,
        )
        if args.device.startswith("cuda"):
            torch.cuda.empty_cache()

    info("Finished fitting model.")
    result = xr.Dataset(
        {
            "gamma": (["strain", "position"], mapest3["gamma"]),
            "rho": (["strain"], mapest3["rho"]),
            "alpha_hyper": ([], mapest3["alpha_hyper"]),
            "pi": (["library_id", "strain"], mapest3["pi"]),
            "epsilon": (["library_id"], mapest3["epsilon"]),
            "rho_hyper": ([], mapest3["rho_hyper"]),
            "epsilon_hyper": ([], mapest3["epsilon_hyper"]),
            "pi_hyper": ([], mapest3["pi_hyper"]),
            "alpha": (["library_id"], mapest3["alpha"]),
            "p_noerr": (["library_id", "position"], mapest3["p_noerr"]),
            "p": (["library_id", "position"], mapest3["p"]),
            "y": (["library_id", "position"], y_obs_ss),
            "m": (["library_id", "position"], m_ss),
            "elbo_trace": (["iteration"], history1),
        },
        coords=dict(
            strain=np.arange(s_collapse),
            position=data_fit.position,
            library_id=data_fit.library_id,
        ),
    )

    if args.outpath:
        info("Saving results.")
        result.to_netcdf(
            args.outpath,
            encoding=dict(
                gamma=dict(dtype="float32", zlib=True, complevel=6),
                pi=dict(dtype="float32", zlib=True, complevel=6),
                p_noerr=dict(dtype="float32", zlib=True, complevel=6),
                p=dict(dtype="float32", zlib=True, complevel=6),
                y=dict(dtype="uint16", zlib=True, complevel=6),
                m=dict(dtype="uint16", zlib=True, complevel=6),
                elbo_trace=dict(dtype="float32", zlib=True, complevel=6),
            ),
        )
