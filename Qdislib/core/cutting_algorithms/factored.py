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

"""A cut circuit's output distribution, kept by fragment.

Over ``n`` qubits a distribution has ``2**n`` outcomes, more than fit in memory soon
after 30 qubits (and in a Python dict soon after 25). A cut circuit's does not need
them: it is a contraction of one small tensor per fragment, over the cut's bonds and
that fragment's own outcomes (``2**(n/k)`` of them for ``k`` even fragments). A
:class:`FactoredDistribution` keeps those tensors and contracts them when asked.
"""

import collections.abc
import itertools

import numpy as np

__all__ = ["FactoredDistribution"]

#: Widest distribution :meth:`FactoredDistribution.to_dict` builds as a dict, and the
#: widest output :func:`circuit_cutting_sampling` returns as one.
MAX_DENSE_WIDTH = 24


def _transform_outcomes(tensor):
    """Walsh-Hadamard transform along the last axis, divided by its size: from a
    fragment's Z-string correlators to its (quasi-)probabilities."""
    a = np.array(tensor, dtype=float)
    shape, n = a.shape, a.shape[-1]
    a = a.reshape(-1, n)
    h = 1
    while h < n:
        for start in range(0, n, 2 * h):
            left = a[:, start : start + h].copy()
            right = a[:, start + h : start + 2 * h].copy()
            a[:, start : start + h] = left + right
            a[:, start + h : start + 2 * h] = left - right
        h *= 2
    return (a / n).reshape(shape)


class FactoredDistribution(collections.abc.Mapping):
    """The output distribution of a cut circuit, by fragment: what
    :func:`~Qdislib.circuit_cutting_sampling` returns over more than 24 qubits.

    A mapping from bitstrings (qubit 0 first, or ``qubits[j]`` at position ``j`` with
    ``qubits=``) to probabilities, each computed when read: ``d["0101..."]``,
    ``d.get(bits, 0.0)``. Nothing here holds ``2**n`` values, unless :meth:`to_dict`
    is asked to::

        d = qd.circuit_cutting_sampling(qc, cuts, shots=100000)   # 40 qubits
        d["0" * 40]              # one outcome
        d.marginal([3, 17])      # any marginal, without running anything again
        d.most_likely(10)        # the ten likeliest outcomes (beam search)
        d.sample(1000, seed=1)   # {bitstring: count}

    With shots these are the estimator's quasi-probabilities: a small one can be
    negative, and they are not rescaled. :meth:`marginal` and :meth:`to_dict` give what
    a dict output gives: absolute values, normalised over the outcomes they return.

    The cutting functions build it, not you. ``width`` is the number of bits of a key.
    """

    def __init__(self, parts, outs, weights, width):
        #: ``(tensor, labels)`` per fragment: an axis per bond, then its outcomes.
        self._parts = [
            (_transform_outcomes(tensor), list(labels)) for tensor, labels in parts
        ]
        #: Per fragment, the position in a key of each bit of its outcomes.
        self._outs = [list(out) for out in outs]
        self._weights = weights
        self.width = int(width)
        # a position in no fragment's output belongs to an idle qubit: always 0
        self._owner = {q: f for f, out in enumerate(self._outs) for q in out}

    def bond_weights(self, label):
        """The weights of a bond's terms (as a plan's, for the contraction)."""
        return self._weights[label]

    def __repr__(self):
        sizes = [part[0].size for part in self._parts]
        return (
            f"FactoredDistribution({self.width} qubits, {len(self._parts)} fragments, "
            f"{sum(sizes)} stored values instead of 2**{self.width})"
        )

    def __len__(self):
        return 2**self.width

    def __bool__(self):  # len() has a limit, and past 62 qubits it raises
        return True

    def __iter__(self):
        for bits in itertools.product("01", repeat=self.width):
            yield "".join(bits)

    def __getitem__(self, bits):
        if len(bits) != self.width or set(bits) - {"0", "1"}:
            raise KeyError(bits)
        return self._partial(dict(enumerate(int(b) for b in bits)))

    def _partial(self, fixed):
        """The total probability of the outcomes with ``fixed`` ``{position: bit}``."""
        from Qdislib.core.cutting_algorithms.mixed_cutting import _contract

        if any(bit and pos not in self._owner for pos, bit in fixed.items()):
            return 0.0  # an idle qubit is never 1
        parts = []
        for (tensor, labels), out in zip(self._parts, self._outs):
            reduced = tensor
            if out:
                reduced = reduced.reshape(reduced.shape[:-1] + (2,) * len(out))
                # the last bit of the local index is the first in axis order: its
                # bit j is at axis -1 - j
                select = [slice(None)] * reduced.ndim
                summed = []
                for j, pos in enumerate(out):
                    axis = reduced.ndim - 1 - j
                    if pos in fixed:
                        select[axis] = fixed[pos]
                    else:
                        summed.append(axis)
                sums = tuple(summed)
                reduced = reduced.sum(axis=sums, keepdims=True) if sums else reduced
                reduced = reduced[tuple(select)].reshape(tensor.shape[:-1] + (-1,))
            parts.append((reduced, labels, None))
        return float(np.sum(_contract(self, parts)))

    def marginal(self, qubits):
        """The distribution over ``qubits`` (positions in this distribution's keys), as
        a dict over them, bit ``j`` of a key being ``qubits[j]``: nothing runs again.

        :raises ValueError: For repeated or unknown positions, or more than 24.
        """
        from Qdislib.core.cutting_algorithms.mixed_cutting import (
            _contract,
            _global_correlators,
        )

        qubits = [int(q) for q in qubits]
        if not qubits or len(set(qubits)) != len(qubits):
            raise ValueError(f"marginal takes distinct positions, got {qubits}")
        if not all(0 <= q < self.width for q in qubits):
            raise ValueError(f"positions go from 0 to {self.width - 1}, got {qubits}")
        if len(qubits) > MAX_DENSE_WIDTH:
            raise ValueError(
                f"a marginal over {len(qubits)} qubits has 2**{len(qubits)} outcomes: "
                f"at most {MAX_DENSE_WIDTH}"
            )
        where = {q: j for j, q in enumerate(qubits)}
        parts, outs = [], []
        for (tensor, labels), out in zip(self._parts, self._outs):
            kept = [j for j, q in enumerate(out) if q in where]
            reduced = tensor
            if out:
                reduced = reduced.reshape(tensor.shape[:-1] + (2,) * len(out))
                gone = tuple(
                    reduced.ndim - 1 - j for j in range(len(out)) if j not in kept
                )
                reduced = reduced.sum(axis=gone) if gone else reduced
                reduced = reduced.reshape(tensor.shape[:-1] + (-1,))
            parts.append((reduced, labels, None))
            outs.append([where[out[j]] for j in kept])
        tensor = _contract(self, parts)
        probabilities = _global_correlators(tensor, outs, len(qubits))
        return _as_dict(probabilities, self._idle_mask(qubits), len(qubits))

    def _idle_mask(self, qubits):
        """The bits, over ``qubits``, of the qubits in no fragment's output."""
        return sum(1 << j for j, q in enumerate(qubits) if q not in self._owner)

    def to_dict(self):
        """Every outcome as a dict, as a smaller circuit's output (2**width entries).

        :raises ValueError: Over more than 24 qubits.
        """
        if self.width > MAX_DENSE_WIDTH:
            raise ValueError(
                f"{self.width} qubits have 2**{self.width} outcomes: read outcomes by "
                "key, or use marginal, most_likely or sample"
            )
        return self.marginal(range(self.width))

    def most_likely(self, k=10, beam=None):
        """The ``k`` likeliest outcomes as ``{bitstring: probability}``, likeliest first.

        A beam search over the qubits, keeping the ``beam`` likeliest prefixes (default
        ``4 * k``, at least 16). It is a search: with enough probability spread thinly
        across prefixes it can miss one.

        :raises ValueError: For ``k`` below 1.
        """
        if k < 1:
            raise ValueError(f"k must be at least 1, got {k}")
        beam = max(4 * k, 16) if beam is None else int(beam)
        prefixes = [({}, 1.0)]
        for position in range(self.width):
            grown = []
            for fixed, _ in prefixes:
                for bit in (0, 1):
                    longer = {**fixed, position: bit}
                    grown.append((longer, self._partial(longer)))
            grown.sort(key=lambda item: -item[1])
            prefixes = grown[:beam]
        best = sorted(prefixes, key=lambda item: -item[1])[:k]
        return {
            "".join(str(fixed[q]) for q in range(self.width)): float(p)
            for fixed, p in best
        }

    def sample(self, shots, seed=None):
        """``shots`` bitstrings drawn from the distribution, as ``{bitstring: count}``:
        qubit by qubit, each conditioned on the ones before (a negative quasi-probability
        counts as zero). About ``width`` contractions per distinct prefix.

        :param shots: How many (not negative).
        :param seed: Seed for the draws.
        :raises ValueError: For negative ``shots``.
        """
        if shots < 0:
            raise ValueError(f"shots must not be negative, got {shots}")
        rng = np.random.default_rng(seed)
        counts = {}

        def descend(fixed, position, left):
            if not left:
                return
            if position == self.width:
                key = "".join(str(fixed[q]) for q in range(self.width))
                counts[key] = counts.get(key, 0) + int(left)
                return
            zero = max(self._partial({**fixed, position: 0}), 0.0)
            one = max(self._partial({**fixed, position: 1}), 0.0)
            total = zero + one
            ones = rng.binomial(left, one / total) if total > 0 else 0
            descend({**fixed, position: 0}, position + 1, left - ones)
            descend({**fixed, position: 1}, position + 1, ones)

        descend({}, 0, int(shots))
        return counts


def _as_dict(probabilities, idle, width):
    """``{bits: p}`` over ``width`` positions (bit ``q`` of the index is position ``q``),
    the bits of ``idle`` never set: absolute values normalised, as a dict output is."""
    values = np.abs(np.asarray(probabilities, dtype=float))
    if idle:
        values = values * ((np.arange(2**width) & idle) == 0)
    norm = float(values.sum()) or 1.0
    return {
        "".join(str((x >> q) & 1) for q in range(width)): float(v) / norm
        for x, v in enumerate(values)
    }


def factored_distribution(plan, parts, marginal=None):
    """The :class:`FactoredDistribution` of contracted-to-be fragment tensors ``parts``
    (``(tensor, labels, out)``), over every qubit or a marginal's."""
    where = None
    width = plan.width
    if marginal is not None:
        where = {q: j for j, q in enumerate(marginal.qubits)}
        width = len(marginal.qubits)
    outs = [out if where is None else [where[q] for q in out] for _, _, out in parts]
    labels = {label for _, part, _ in parts for label in part[:-1]}
    weights = {label: plan.bond_weights(label) for label in labels}
    return FactoredDistribution(
        [(tensor, part) for tensor, part, _ in parts], outs, weights, width
    )
