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


def rss(x, y):
    return np.sqrt(np.sum((x - y) ** 2))


def binary_entropy(p):
    q = 1 - p
    return -p * np.log2(p) - q * np.log2(q)


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
    plt.title(f"+{min_loss}")
    plt.yscale("log")
    return plt.gca()


def mean_residual_count(expect_frac, obs_count, m):
    frac_obs = obs_count / m
    out = np.abs(((frac_obs - expect_frac)))
    out[np.isnan(out)] = 0
    return (out * m).sum() / m.sum()


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
    model, data={}, dtype=torch.float32, device="cpu", **kwargs,
):
    data = {
        k: torch.tensor(v, dtype=dtype, device=device) for k, v in data.items()
    }
    return partial(
        pyro.condition(model, data=data), dtype=dtype, device=device, **kwargs,
    )


def find_map(
    model,
    lag=10,
    stop_at=1.0,
    max_iter=int(1e5),
    learning_rate=1e-0,
    clip_norm=100.0,
):
    guide = pyro.infer.autoguide.AutoLaplaceApproximation(model)
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
            if i % 1 == 0:
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

    # Gather MAP from parameter-store
    mapest = {
        k: v.detach().cpu().numpy().squeeze()
        for k, v in pyro.infer.Predictive(
            model, guide=guide, num_samples=1
        )().items()
    }
    return mapest, np.array(history)


def parse_args(argv):
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # Input
    p.add_argument(
        "pileup",
        help="""
Single, fully processed, pileup table in NetCDF format
with the following dimensions:
    * library_id
    * position
    * read
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
        "--gamma-hyper",
        metavar="FLOAT",
        default=1e-2,
        type=float,
        help=("Ambiguity regularization parameter."),
    )
    p.add_argument(
        "--epsilon-hyper", metavar="FLOAT", default=0.01, type=float
    )
    p.add_argument(
        "--alpha",
        metavar="FLOAT",
        default=100.0,
        type=float,
        help=("Concentration parameter of BetaBinomial observation."),
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

    if args.outpath is None:
        args.outpath = args.pileup + ".sfacts.nc"

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
        message="torch.tensor results are registered as constants",
        category=torch.jit.TracerWarning,
        # module="trace_elbo",  # FIXME: What is the correct regex for module?
        lineno=65,
    )

    args = parse_args(sys.argv[1:])
    info(args)

    info(f"Setting random seed: {args.random_seed}")
    np.random.seed(args.random_seed)

    info("Loading input data.")
    data = xr.open_dataarray(args.pileup).squeeze()
    info(f"Input data shape: {data.sizes}.")
    data = data.sum("read")

    info("Filtering positions.")
    minor_allele_incid = (data > 0).mean("library_id").min("allele")
    informative_positions = idxwhere(
        minor_allele_incid.to_series() > args.incid_thresh
    )
    npos_available = len(informative_positions)
    info(
        f"Found {npos_available} informative positions with minor "
        f"allele incidence of >{args.incid_thresh}"
    )
    npos = min(args.npos, npos_available)
    info(f"Randomly sampling {npos} positions.")
    position_ss = np.random.choice(
        informative_positions, size=npos, replace=False,
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

    info("Building conditioned model.")
    data_fit = data.sel(library_id=suff_cvrg_samples, position=position_ss)
    m = data_fit.sum("allele")
    n, g = m.shape
    y_obs = data_fit.sel(allele="alt")
    s = args.nstrains
    info(f"Model shape: n={n}, g={g}, s={s}")
    model_fit = conditioned_model(
        model,
        data=dict(
            alpha=np.ones(n) * args.alpha,
            epsilon_hyper=args.epsilon_hyper,
            pi_hyper=args.pi_hyper / s,
            rho_hyper=args.rho_hyper,
            y=y_obs.values,
        ),
        s=s,
        m=m.values,
        gamma_hyper=args.gamma_hyper,
        dtype=torch.float32,
        device=args.device,
    )

    info("Fitting model.")
    mapest, history = find_map(
        model_fit,
        lag=args.lag,
        stop_at=args.stop_at,
        learning_rate=args.learning_rate,
        max_iter=args.max_iter,
        clip_norm=args.clip_norm,
    )
    if args.device.startswith('cuda'):
        torch.cuda.empty_cache()

    info("Finished fitting model.")
    result = xr.Dataset(
        {
            "gamma": (["strain", "position"], mapest["gamma"]),
            "rho": (["strain"], mapest["rho"]),
            "alpha_hyper": ([], mapest["alpha_hyper"]),
            "pi": (["library_id", "strain"], mapest["pi"]),
            "epsilon": (["library_id"], mapest["epsilon"]),
            "rho_hyper": ([], mapest["rho_hyper"]),
            "epsilon_hyper": ([], mapest["epsilon_hyper"]),
            "pi_hyper": ([], mapest["pi_hyper"]),
            "alpha": (["library_id"], mapest["alpha"]),
            "p_noerr": (["library_id", "position"], mapest["p_noerr"]),
            "p": (["library_id", "position"], mapest["p"]),
            "y": (["library_id", "position"], y_obs),
            "m": (["library_id", "position"], m),
            "elbo_trace": (["iteration"], history),
        },
        coords=dict(
            strain=np.arange(s),
            position=data_fit.position,
            library_id=data_fit.library_id,
        ),
    )

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
