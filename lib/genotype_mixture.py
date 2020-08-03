#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pymc3 as pm
from pymc3.distributions.transforms import t_stick_breaking, logodds
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from lib.util import info
from lib.pileup import list_bases, get_pileup_dims
import tqdm
import theano.tensor as tt
from itertools import product


BASE_TO_IUPAC = {frozenset('a'): 'a',
                 frozenset('c'): 'c',
                 frozenset('g'): 'g',
                 frozenset('t'): 't',
                 frozenset('ag'): 'r',
                 frozenset('ct'): 'y',
                 frozenset('cg'): 's',
                 frozenset('at'): 'w',
                 frozenset('gt'): 'k',
                 frozenset('ac'): 'm',
                 frozenset('cgt'): 'b',
                 frozenset('agt'): 'd',
                 frozenset('act'): 'h',
                 frozenset('acg'): 'v',
                 frozenset('acgt'): 'n'}


stick_breaking = t_stick_breaking(1e-10)


def build_polyallelic_model(n, g, s, a=4):
    with pm.Model() as model:
        # Fraction
        pi = pm.Dirichlet('pi', a=np.ones(s), shape=(n, s),
                          transform=stick_breaking
                          )
        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pm.Potential('heterogeneity_penalty',
                     # NOTE: we take the mean sqrt over the first axis so that
                     # it's invariant to the number of samples, but sum over
                     # the second axis so that it's affected by the actual
                     # number of latent strains, not the number allowed.
                     -(pm.math.sqrt(pi).mean(0).sum()**2)
                     * pi_hyper)

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        pm.Potential('diversity_penalty',
                     # NOTE: we take the mean over the first axis *before*
                     # taking the sqrt because we're interested in mean
                     # abundances of each strain.  This means that it *will* be
                     # affected by the choice of samples.
                     # We later sum over (instead of taking the mean over) the
                     # second axis so that it's affected by the actual number
                     # of latent strains, not the number allowed (s).
                     -(pm.math.sqrt(pi.mean(0)).sum()**2)
                     * rho_hyper)

        # Genotype
        gamma_ = pm.Dirichlet('gamma_', a=np.ones(a), shape=(g * s, a),
                              transform=stick_breaking
                              )
        gamma = pm.Deterministic('gamma', gamma_.reshape((g, s, a)))
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        pm.Potential('ambiguity_penalty',
                     # NOTE: we're taking the norm over the third axis.
                     # Then, after returning to the original scale, we take the
                     # mean over the first axis so that it is invariant to
                     # numbers of positions and finally we sum over the second
                     # axis so that each strain has an equal impact, regardless
                     # of the number of strains.
                     -(pm.math.sqrt(gamma).sum(2)**2).mean(0).sum(0)
                     * gamma_hyper)
        # NOTE: As a general rule, we sum over those dimensions where we want
        # to force all of the weight to one element (e.g. alleles or strains).
        # Conversely, we take the mean of dimensions that we want to be
        # invariant to decisions/facts that are independent of the goodness of
        # fit (e.g. number of positions or samples).

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        err_base_prob = tt.ones((n, g, a)) / a
        p_with_error = (true_p * (1 - epsilon_)) + (err_base_prob * epsilon_)

        # Observation
        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.Multinomial('data',
                       p=p_with_error.reshape((-1, a)),
                       n=observed.reshape((-1, a)).sum(1),
                       observed=observed)

    return model


def build_biallelic_model(g, n, s):
    a = 2

    with pm.Model() as model:
        # Fraction
        pi = pm.Dirichlet('pi', a=np.ones(s), shape=(n, s),
                          transform=stick_breaking,
                          )
        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pm.Potential('heterogeneity_penalty',
                     -(pm.math.sqrt(pi).sum(0).sum()**2) * pi_hyper)

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        pm.Potential('diversity_penalty',
                     -(pm.math.sqrt(pi.sum(0)).sum()**2)
                     * rho_hyper)

        # Genotype
        gamma_ = pm.Uniform('gamma_', 0, 1, shape=(g * s, 1))
        gamma = pm.Deterministic('gamma', (pm.math.concatenate([gamma_, 1 - gamma_], axis=1)
                                           .reshape((g, s, a))))
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        pm.Potential('ambiguity_penalty',
                     -(pm.math.sqrt(gamma).sum(2)**2).sum(0).sum(0)
                     * gamma_hyper)

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        err_base_prob = tt.ones((n, g, a)) / a
        p_with_error = (true_p * (1 - epsilon_)) + (err_base_prob * epsilon_)

        # Observation
        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.Binomial('data',
                    p=p_with_error.reshape((-1, a))[:,0],
                    n=observed.reshape((-1, a)).sum(1),
                    observed=observed[:,0])

    return model


def build_biallelic_model2(g, n, s):
    # EXPERIMENTAL: Dirichlet priors on latent variables.
    a = 2

    with pm.Model() as model:
        # Fraction

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        rho = pm.Dirichlet('rho', a=np.ones(s) * pm.math.exp(-rho_hyper), shape=(s,),
                           transform=stick_breaking,
                           )

        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pi = pm.Dirichlet('pi', a=rho * pm.math.exp(-pi_hyper), shape=(n, s),
                          transform=stick_breaking,
                          )

        # Genotype
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        gamma_ = pm.Beta('gamma_',
                         pm.math.exp(-gamma_hyper), pm.math.exp(-gamma_hyper),
                         shape=(g * s, 1))
        gamma = pm.Deterministic('gamma', (pm.math.concatenate([gamma_, 1 - gamma_], axis=1)
                                           .reshape((g, s, a))))

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        err_base_prob = tt.ones((n, g, a)) / a
        p_with_error = (true_p * (1 - epsilon_)) + (err_base_prob * epsilon_)

        # Observation
        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.Binomial('data',
                    p=p_with_error.reshape((-1, a))[:,0],
                    n=observed.reshape((-1, a)).sum(1),
                    observed=observed[:,0])

    return model


def build_biallelic_model3(g, n, s):
    # EXPERIMENTAL: Observations overdispersed as a BetaBinom w/ concentrations
    # 10.
    a = 2

    with pm.Model() as model:
        # Fraction
        pi = pm.Dirichlet('pi', a=np.ones(s), shape=(n, s),
                          transform=stick_breaking,
                          )
        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pm.Potential('heterogeneity_penalty',
                     -(pm.math.sqrt(pi).sum(0).sum()**2) * pi_hyper)

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        pm.Potential('diversity_penalty',
                     -(pm.math.sqrt(pi.sum(0)).sum()**2)
                     * rho_hyper)

        # Genotype
        gamma_ = pm.Uniform('gamma_', 0, 1, shape=(g * s, 1))
        gamma = pm.Deterministic('gamma', (pm.math.concatenate([gamma_, 1 - gamma_], axis=1)
                                           .reshape((g, s, a))))
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        pm.Potential('ambiguity_penalty',
                     -(pm.math.sqrt(gamma).sum(2)**2).sum(0).sum(0)
                     * gamma_hyper)

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        err_base_prob = tt.ones((n, g, a)) / a
        p_with_error = (true_p * (1 - epsilon_)) + (err_base_prob * epsilon_)

        # Observation
        _p = p_with_error.reshape((-1, a))[:, 0]
        # Overdispersion term
        # alpha = pm.Gamma('alpha', mu=100, sigma=5)
        # TODO: Figure out how to also fit this term.
        # FIXME: Do I want the default to be a valid value?
        #  Realistic or close to asymptotic?
        alpha = pm.Data('alpha', value=1000)

        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.BetaBinomial('data',
                        alpha=_p * alpha,
                        beta=(1 - _p) * alpha,
                        n=observed.reshape((-1, a)).sum(1),
                        observed=observed[:,0])

    return model


def build_biallelic_model4(g, n, s):
    # EXPERIMENTAL: Underlying allele freqs overdispersed before errors.
    a = 2

    with pm.Model() as model:
        # Fraction
        pi = pm.Dirichlet('pi', a=np.ones(s), shape=(n, s),
                          transform=stick_breaking,
                          )
        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pm.Potential('heterogeneity_penalty',
                     -(pm.math.sqrt(pi).sum(0).sum()**2) * pi_hyper)

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        pm.Potential('diversity_penalty',
                     -(pm.math.sqrt(pi.sum(0)).sum()**2)
                     * rho_hyper)

        # Genotype
        gamma_ = pm.Uniform('gamma_', 0, 1, shape=(g * s, 1))
        gamma = pm.Deterministic(
                'gamma',
                (pm.math.concatenate([gamma_, 1 - gamma_], axis=1)
                 .reshape((g, s, a)))
                                )
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        pm.Potential('ambiguity_penalty',
                     -(pm.math.sqrt(gamma).sum(2)**2).sum(0).sum(0)
                     * gamma_hyper)

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Overdispersion term
        # TODO: Consider making it different between samples.  How to shape?
        alpha = pm.ChiSquared('alpha', nu=50)
        _true_p = true_p.reshape((-1, a))[:, 0]
        _true_q = 1 - _true_p
        overdispersed_p_ = pm.Beta('overdispersed_p_',
                                   alpha=_true_p * alpha,
                                   beta=_true_q * alpha,
                                   shape=(n * g,))
        overdispersed_p = pm.Deterministic(
                'overdispersed_p',
                pm.math.concatenate([overdispersed_p_.reshape((-1, 1)),
                                    1 - overdispersed_p_.reshape((-1, 1))],
                                    axis=1).reshape((n, g, a))
                                          )

        # Sequencing error
        # epsilon_hyper = pm.Gamma('epsilon_hyper', alpha=100, beta=1)
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        p_with_error = (overdispersed_p * (1 - epsilon_) +
                        (1 - overdispersed_p) * (epsilon_ / (a - 1)))

        # Observation
        # _p = p_with_error.reshape((-1, a))[:,0]
        # _q = 1 - _p

        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.Binomial('data',
                    p=p_with_error.reshape((-1, a))[:, 0],
                    n=observed.reshape((-1, a)).sum(1),
                    observed=observed[:, 0])

    return model


def build_biallelic_model5(g, n, s):
    # EXPERIMENTAL: Dirichlet priors on latent variables.
    a = 2

    with pm.Model() as model:
        # Fraction

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        rho = pm.Dirichlet('rho', a=np.ones(s) * pm.math.exp(-rho_hyper),
                           shape=(s,),
                           transform=stick_breaking,
                           )

        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pi = pm.Dirichlet('pi', a=rho * pm.math.exp(-pi_hyper), shape=(n, s),
                          transform=stick_breaking,
                          )

        # Genotype
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        gamma_loc_ = pm.Uniform('gamma_loc_', 0, 1, shape=(g * s, 1))
        gamma_ = pm.Beta('gamma_',
                         alpha=gamma_loc_ * pm.math.exp(-gamma_hyper),
                         beta=(1 - gamma_loc_) * pm.math.exp(-gamma_hyper),
                         shape=(g * s, 1))
        gamma = pm.Deterministic('gamma',
                                 (pm.math.concatenate([gamma_, 1 - gamma_],
                                                      axis=1)
                                  .reshape((g, s, a))))

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        p_with_error = (true_p * (1 - epsilon_) +
                        (1 - true_p) * (epsilon_ / (a - 1)))

        # Observation
        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.Binomial('data',
                    p=p_with_error.reshape((-1, a))[:, 0],
                    n=observed.reshape((-1, a)).sum(1),
                    observed=observed[:, 0])

    return model

def build_biallelic_model6(g, n, s):
    # EXPERIMENTAL: Observations overdispersed as a BetaBinom w/ concentrations
    # 10.
    a = 2

    with pm.Model() as model:
        # Fraction
        pi = pm.Dirichlet('pi', a=np.ones(s), shape=(n, s),
                          transform=stick_breaking,
                          )
        pi_hyper = pm.Data('pi_hyper', value=0.0)
        pm.Potential('heterogeneity_penalty',
                     -(pm.math.sqrt(pi).sum(0).sum()**2) * pi_hyper)

        rho_hyper = pm.Data('rho_hyper', value=0.0)
        pm.Potential('diversity_penalty',
                     -(pm.math.sqrt(pi.sum(0)).sum()**2)
                     * rho_hyper)

        # Genotype
        gamma_ = pm.Uniform('gamma_', 0, 1, shape=(g * s, 1))
        gamma = pm.Deterministic('gamma', (pm.math.concatenate([gamma_, 1 - gamma_], axis=1)
                                           .reshape((g, s, a))))
        gamma_hyper = pm.Data('gamma_hyper', value=0.0)
        pm.Potential('ambiguity_penalty',
                     -((pm.math.sqrt(gamma).sum(2)**2).sum(0) * pi.sum(0)).sum(0)
                     * gamma_hyper)

        # Product of fraction and genotype
        true_p = pm.Deterministic('true_p', pm.math.dot(pi, gamma))

        # Sequencing error
        epsilon_hyper = pm.Data('epsilon_hyper', value=100)
        epsilon = pm.Beta('epsilon', alpha=2, beta=epsilon_hyper,
                          shape=n)
        epsilon_ = epsilon.reshape((n, 1, 1))
        err_base_prob = tt.ones((n, g, a)) / a
        p_with_error = (true_p * (1 - epsilon_)) + (err_base_prob * epsilon_)

        # Observation
        _p = p_with_error.reshape((-1, a))[:, 0]
        # Overdispersion term
        # alpha = pm.Gamma('alpha', mu=100, sigma=5)
        # TODO: Figure out how to also fit this term.
        # FIXME: Do I want the default to be a valid value?
        #  Realistic or close to asymptotic?
        alpha = pm.Data('alpha', value=1000)
        # alpha = pm.Gamma('alpha', mu=100, sigma=5)

        observed = pm.Data('observed', value=np.empty((g * n, a)))
        pm.BetaBinomial('data',
                        alpha=_p * alpha,
                        beta=(1 - _p) * alpha,
                        n=observed.reshape((-1, a)).sum(1),
                        observed=observed[:,0])

    return model

def pileup_to_model_input(pileup):
    g, n, a = get_pileup_dims(pileup)
    return pileup.values.reshape((g, n, a)).swapaxes(1, 0)


def compile_start_point(model, gamma=None, pi=None, epsilon=None):
    point = model.test_point
    if gamma is not None:
        point['gamma__stickbreaking__'] = \
                stick_breaking.forward(gamma).eval()
    if pi is not None:
        point['pi_stickbreaking__'] = \
                stick_breaking.forward(pi).eval()
    if epsilon is not None:
        point['epsilon_logodds__'] = \
                logodds.forward(epsilon).eval()
    return point


def optimize_parameters(model, observed,
                        gamma_hyper=0, pi_hyper=0, rho_hyper=0,
                        start=None, verbose=1, **kwargs):
    n, g, a = observed.shape
    model.gamma_hyper.set_value(gamma_hyper)
    model.pi_hyper.set_value(pi_hyper)
    model.rho_hyper.set_value(rho_hyper)
    model.observed.set_value(observed.reshape((-1, a)))
    i = 1
    while True:
        mapest, optim = pm.find_MAP(model=model,
                                    return_raw=True, maxeval=10000,
                                    start=start, progressbar=(verbose >= 2),
                                    **kwargs)
        # optim returns as None if it reaches maxeval or it gets a SIGKILL.
        if optim is not None:
            logp = model.logp(mapest)
            assert optim.success and np.isfinite(logp), \
                f'Optimization failed:\n{optim}'
            break

        start = mapest
        i += 1
        if verbose > 1:
            info(f'MAP estimate has not yet converged. '
                 f'Starting round {i} of gradient descent.')

    return mapest


def gamma_to_genotype_table(gamma, pileup):
    genotype = []
    for i in range(gamma.shape[1]):
        genotype_i = pd.DataFrame(gamma[:, i, :],
                                  index=pileup.index,
                                  columns=list_bases(pileup))
        genotype_i.columns.name = 'base'
        genotype_i = genotype_i.stack()
        genotype_i.name = i
        genotype.append(genotype_i)
    genotype = (pd.concat(genotype, axis='columns')
                  .drop_duplicates()
                  .rename_axis(columns='strain')
                  .unstack('base'))
    return genotype


def gamma_plot(gamma, savepath=None, width_per_s=0.4, height_per_g=0.1):
    g, s, a = gamma.shape
    max_height = 100
    fig, axs = plt.subplots(1, s, sharex=True, sharey=True,
                            figsize=(s * width_per_s,
                                     min(g * height_per_g, max_height)))
    for i, ax in enumerate(axs):
        ax.imshow(gamma[:, i, :], vmin=0, vmax=1, aspect='auto')
        ax.set_title(i)
    ax.set_xticks([])
    ax.set_yticks([])
    if savepath is not None:
        fig.savefig(savepath)


def genotype_plot(geno, savepath=None, width_per_s=0.4, height_per_g=0.1):
    strain_list = geno.columns.get_level_values(level='strain').unique()
    base_list = geno.columns.get_level_values(level='base').unique()
    g = geno.shape[0]
    s = len(strain_list)
    a = len(base_list)
    assert geno.shape == (g, s * a)
    max_height = 100
    fig, axs = plt.subplots(1, s, sharex=True, sharey=True,
                            figsize=(s * width_per_s,
                                     min(g * height_per_g, max_height)))
    for i, ax in zip(strain_list, axs):
        ax.imshow(geno[i], vmin=0, vmax=1, aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(i)
    if savepath is not None:
        fig.savefig(savepath)


def pi_plot(pi, savepath=None, width_per_s=0.4, height_per_n=0.05, pwr=1):
    n, s = pi.shape
    fig, ax = plt.subplots(figsize=(s * width_per_s, n * height_per_n))
    ax.imshow(pi, vmin=0, vmax=1, aspect='auto',
              norm=mpl.colors.PowerNorm(pwr))
    ax.set_xticks(range(s))
    ax.set_yticks([])
    if savepath is not None:
        fig.savefig(savepath)


def calculate_deviation(pi, gamma, observed):
    true_p = np.dot(pi, gamma)
    total_counts = observed.sum(2, keepdims=True)
    expect_counts = true_p * total_counts
    deviation = np.abs(observed - expect_counts).sum(-1)
    return deviation


def ambiguity_index(gamma, pi=None):
    if pi is None:
        pi = np.ones(gamma.shape[1])
    return (((np.sqrt(gamma).sum(2)**2) - 1).mean(0)
            * pi.mean(0)).sum()


def deviation_plot(pi, gamma, observed,
                   savepath=None,
                   width_per_n=0.1, height_per_g=0.1,
                   pwr=1/3, vmax=10):
    true_p = np.dot(pi, gamma)
    n, g, a = true_p.shape
    mean_sample_coverage = (observed.mean(1).sum(1, keepdims=True))
    deviation = calculate_deviation(pi, gamma, observed)
    scaled_deviation = deviation / mean_sample_coverage

    fig, ax = plt.subplots(figsize=(n * width_per_n, g * height_per_g))
    sns.heatmap(scaled_deviation.T, norm=mpl.colors.PowerNorm(pwr),
                vmin=0, vmax=vmax,
                cbar_kws={"ticks": [0, 1, 10, vmax]},
                ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    if savepath is not None:
        fig.savefig(savepath)


def genotype_to_fasta(genotype, path):
    with open(path, 'w') as f:
        for strain, geno in (genotype.stack('strain')
                                     .idxmax(1)
                                     .unstack('strain')
                                     .iteritems()):
            algn = ''.join(geno)
            print(f'>{strain}\n{algn}', file=f)


def fuzzy_allele_to_iupac(acgt, thresh=0.5):
    priority = acgt.sort_values(ascending=False).cumsum()
    possible = priority.loc[:(priority > thresh).idxmax()]
    return BASE_TO_IUPAC[frozenset(possible.index)]


def genotype_to_iupac(geno, thresh=0.5, progress=True):
    if progress:
        tqdm.pandas()
        out = (geno.stack('strain')
               .progress_apply(fuzzy_allele_to_iupac, axis='columns')
               .unstack('strain'))
    else:
        out = (geno.stack('strain')
               .apply(fuzzy_allele_to_iupac, axis='columns')
               .unstack('strain'))

    return out


def simulate_genotype(ancestral, mutation_rate):
    g, a = ancestral.shape
    ancestral = ancestral.astype(bool)
    assert (ancestral.sum(1) == 1).all()

    prob_genotype = (ancestral * (1 - mutation_rate)
                     + (1 - ancestral) * (mutation_rate / 3))
    assert np.allclose(prob_genotype.sum(1), 1)

    genotype = np.empty_like(prob_genotype)
    for i in range(g):
        genotype[i, :] = np.random.multinomial(1, prob_genotype[i, :])
    genotype = genotype.astype(bool)
    assert (genotype.sum(1) == 1).all()

    return genotype.astype(int)


def simulate_pileup(genotype, frac, cvrg, err_rate):
    g, s, a = genotype.shape
    n, s2 = frac.shape
    assert s == s2
    g, n2 = cvrg.shape
    assert n == n2

    frac_allele = frac @ genotype
    assert frac_allele.shape == (g, n, a)
    prob_allele = (frac_allele * (1 - err_rate)
                   + (1 - frac_allele) * (err_rate / (a - 1)))
    assert np.allclose(prob_allele.sum(2), 1)

    tally = np.empty_like(prob_allele)
    for i in range(g):
        for j in range(n):
            tally[i, j, :] = np.random.multinomial(cvrg[i, j],
                                                   prob_allele[i, j, :])
    assert (tally.sum(2) == cvrg).all()

    pileup = []
    for i in range(n):
        pileup.append(pd.DataFrame(tally[:, i, :]).stack())
    pileup = (pd.concat(pileup, axis=1)
              .unstack(1)
              .rename_axis(columns=('sample_id', 'base'), index='position'))
    return pileup
