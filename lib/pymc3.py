import pymc3 as pm
import numpy as np
from typing import Dict
import matplotlib.pyplot as plt


def _debug(*args, **kwargs):
    pass
    # print(*args, file=sys.stderr, **kwargs)


class RingBuff(pm.backends.ndarray.NDArray):
    """NDArray trace object

    Parameters
    ----------
    name : str
        Name of backend. This has no meaning for the NDArray backend.
    model : Model
        If None, the model is taken from the `with` context.
    keep_n : int
        Number of samples to keep.
    vars : list of variables
        Sampling values will be stored for these variables. If None,
        `model.unobserved_RVs` is used.
    """

    supports_sampler_stats = True

    def __init__(self, name=None, model=None, vars=None, keep_n=10,
                 test_point=None):
        super().__init__(name, model, vars, test_point)
        self.keep_n = keep_n

    def setup(self, draws, chain, sampler_vars=None) -> None:
        """Perform chain-specific setup.

        Parameters
        ----------
        draws : int
            Expected number of draws
        chain : int
            Chain number
        sampler_vars : list of dicts
            Names and dtypes of the variables that are
            exported by the samplers.
        """
        _debug(f"Setting up ring buffer backend of size {self.keep_n}.")
        super(pm.backends.ndarray.NDArray, self).setup(draws, chain, sampler_vars)
        self.chain = chain
        _debug(f"I am chain {chain}.")
        if self.samples:  # Concatenate new array if chain is already present.
            _debug("Concatenating old samples.")
            old_draws = len(self)
            self.draws = old_draws + draws
            self.draw_idx = old_draws
            for varname, shape in self.var_shapes.items():
                old_var_samples = self.samples[varname]
                _debug(f"Initializing container for {varname} of shape {shape} which has old samples.")
                new_var_samples = np.empty((max(0, self.keep_n - old_draws),)
                                           + shape,
                                           self.var_dtypes[varname])
                _debug(f"Concatenating old samples to {varname}.")
                self.samples[varname] = np.concatenate((old_var_samples,
                                                        new_var_samples),
                                                       axis=0)
                _debug(f"Finished concatenating old samples for {varname}.")
        else:  # Otherwise, make empty arrays for each variable.
            self.draws = draws
            for varname, shape in self.var_shapes.items():
                _debug(f"Initializing container for {varname} of shape {shape}")
                self.samples[varname] = \
                        np.empty((self.keep_n, ) + shape,
                                 dtype=self.var_dtypes[varname])

        if sampler_vars is None:
            return

        if self._stats is None:
            self._stats = []
            for sampler in sampler_vars:
                data = dict()  # type: Dict[str, np.ndarray]
                self._stats.append(data)
                for varname, dtype in sampler.items():
                    data[varname] = np.empty(draws, dtype=dtype)
        else:
            for data, vars in zip(self._stats, sampler_vars):
                if vars.keys() != data.keys():
                    raise ValueError("Sampler vars can't change")
                old_draws = len(self)
                for varname, dtype in vars.items():
                    old = data[varname]
                    new = np.empty(draws, dtype=dtype)
                    data[varname] = np.concatenate([old, new])

    def record(self, point, sampler_stats=None) -> None:
        """Record results of a sampling iteration.

        Parameters
        ----------
        point : dict
            Values mapped to variable names
        """
        for varname, value in zip(self.varnames, self.fn(point)):
            self.samples[varname][self.draw_idx % self.keep_n] = value

        if self._stats is not None and sampler_stats is None:
            raise ValueError("Expected sampler_stats")
        if self._stats is None and sampler_stats is not None:
            raise ValueError("Unknown sampler_stats")
        if sampler_stats is not None:
            for data, vars in zip(self._stats, sampler_stats):
                for key, val in vars.items():
                    data[key][self.draw_idx % self.keep_n] = val
        self.draw_idx += 1

    def close(self):
        if self.draw_idx == self.draws:
            return
        elif self.draw_idx < self.keep_n:
            # Remove trailing zeros if interrupted before completing enough
            # draws.
            self.samples = {var: vtrace[:self.draw_idx]
                            for var, vtrace in self.samples.items()}
        else:
            # Rearrange the trace given the pointer location.
            self.samples = \
                    {var: np.concatenate([vtrace[self.draw_idx % self.keep_n:],
                                          vtrace[:self.draw_idx % self.keep_n]
                                          ],
                                         axis=0)
                     for var, vtrace in self.samples.items()}
        if self._stats is not None:
            self._stats = [
                {var: trace[:self.draw_idx] for var, trace in stats.items()}
                for stats in self._stats]

    def __len__(self):
        if not self.samples:  # `setup` has not been called.
            return 0
        return min(self.draw_idx, self.keep_n)


def split_sampler_traces(trace, statname):
    nsteps = len(trace)
    nchains = len(trace.chains)
    varshape = trace[statname][0].shape
    out = np.empty((nchains, nsteps, *varshape))
    for i in trace.chains:
        chain_start = i * nsteps
        chain_stop = (i + 1) * nsteps
        out[i] = trace[statname][chain_start:chain_stop]
    return out


def trace_stat_plot(trace, statname,
                    savepath=None, exclude=[], skip_frac=0.2):
    nsteps = len(trace)
    skip = int(np.floor(nsteps * skip_frac))
    fig, ax = plt.subplots()
    stat = getattr(trace, statname)
    for i in trace.chains:
        if i in exclude:
            continue
        chain_start = i*nsteps
        chain_stop = i*nsteps + nsteps
        ax.plot(range(skip, nsteps),
                stat[chain_start:chain_stop][skip:], label=i)
    ax.legend(bbox_to_anchor=(1, 1))
    ax.set_title(statname)
    if savepath is not None:
        fig.savefig(savepath)
