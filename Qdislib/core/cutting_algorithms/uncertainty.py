#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# -*- coding: utf-8 -*-

"""Error bars on a reconstructed value, by bootstrap over the shots.

Each fragment circuit's counts are resampled (a multinomial draw of the same
number of shots from what was measured), with the ``num_samples`` draws resampled
the same way, and the cheap reconstruction is done again; the spread of those
values is the standard error. It costs reconstruction time, not device time, and
holds for every cut type, Pauli sums and drawn configurations.
"""

import math

import numpy as np

__all__ = ["Estimate", "shots_for"]

#: Resamples when ``error=True``.
DEFAULT_RESAMPLES = 200


class Estimate(float):
    """A reconstructed value with its standard error: a ``float`` (it is used as
    one everywhere) that also carries ``error`` and the ``shots`` it took.

    ``error`` is the bootstrap standard error over the shots (and the drawn
    configurations with ``num_samples``); 0 for an exact run.
    """

    def __new__(cls, value, error=0.0, shots=0):
        self = float.__new__(cls, value)
        self.error = float(error)
        self.shots = int(shots)
        return self

    def __reduce__(self):
        return (Estimate, (float(self), self.error, self.shots))

    def __repr__(self):
        return f"{float(self):.6g} ± {self.error:.2g}"

    __str__ = __repr__


def resamples(error):
    """How many bootstrap resamples ``error=`` asks for (0: none)."""
    if error is None or error is False:
        return 0
    if error is True:
        return DEFAULT_RESAMPLES
    if isinstance(error, int) and error >= 2:
        return error
    raise ValueError(
        f"error= takes True, False or a number of resamples >= 2, got {error!r}"
    )


def _resample_counts(corr, rng):
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

    if not isinstance(corr, _JointCounts):
        return corr  # exact: nothing to resample
    return corr.resampled(rng)


def _resample_draws(draws, rng):
    """The draws resampled within their strata (the allocated ones keep their
    counts; a group of draws taken at random is redrawn with its total)."""
    if draws is None:
        return None
    strata = getattr(draws, "strata", None) or [0] * len(draws)
    counts = np.array([d[4] for d in draws], dtype=float)
    new = np.zeros(len(draws), dtype=np.int64)
    for group in set(strata):
        where = np.array([g == group for g in strata])
        total = counts[where].sum()
        new[where] = rng.multinomial(int(total), counts[where] / total)
    return [d[:4] + (int(c),) for d, c in zip(draws, new) if c]


def total_shots(results):
    """The shots behind a reconstruction's results (0 for exact ones)."""
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

    return int(
        sum(
            sum(corr.joint.values())
            for corr, _ in results
            if isinstance(corr, _JointCounts)
        )
    )


def bootstrap(evaluate, results, draws, count, seed=None):
    """``(standard error, shots)`` of ``evaluate(results, draws)`` by ``count``
    resamples of the counts (and the draws); 0 when nothing is random."""
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

    shots = total_shots(results)
    random = shots > 0 or (draws is not None and len(draws) > 0)
    if not count or not random:
        return 0.0, shots
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(count):
        sample = [
            (
                (_resample_counts(corr, rng), original)
                if isinstance(corr, _JointCounts)
                else (corr, original)
            )
            for corr, original in results
        ]
        values.append(evaluate(sample, _resample_draws(draws, rng)))
    return float(np.std(values, ddof=1)), shots


def combine(coefficients, values):
    """``sum(c * v)``, an :class:`Estimate` (errors in quadrature: the groups of a
    Pauli sum run separate circuits) when any value is one."""
    total = float(sum(c * float(v) for c, v in zip(coefficients, values)))
    if not any(isinstance(v, Estimate) for v in values):
        return total
    error = math.sqrt(
        sum((c * getattr(v, "error", 0.0)) ** 2 for c, v in zip(coefficients, values))
    )
    shots = sum(getattr(v, "shots", 0) for v in values)
    return Estimate(total, error, shots)


def shots_for(circuit, cuts, epsilon, delta=0.05, **options):
    """Total shots that guarantee ``|estimate - <O>| <= epsilon`` with probability
    ``1 - delta`` for **any** circuit with this cut, by Hoeffding's bound on the
    quasiprobability Monte-Carlo estimate: ``2 gamma**2 / epsilon**2 * ln(2/delta)``.

    A worst-case number: with drawn configurations the actual factor over the uncut
    circuit is at least ``gamma**2`` (``(gamma**2 - <O>**2) / (1 - <O>**2)``), while
    running every term with its share of the shots can need fewer. ``options`` go to
    :func:`~Qdislib.cutting_overhead` (``hadamard``, ``optimized``, ``joint``).
    """
    from Qdislib.core.cutting_algorithms.mixed_cutting import cutting_overhead

    if not 0 < epsilon:
        raise ValueError("epsilon must be positive")
    if not 0 < delta < 1:
        raise ValueError("delta must be in (0, 1)")
    gamma = cutting_overhead(circuit, cuts, **options)["gamma"]
    return int(math.ceil(2 * gamma**2 / epsilon**2 * math.log(2 / delta)))
