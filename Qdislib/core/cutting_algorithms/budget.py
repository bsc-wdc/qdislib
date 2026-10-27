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

"""``total_shots=``: one shot budget split across the fragment circuits of a cut.

In two rounds:

1. a pilot (``PILOT`` of the budget) split in proportion to each circuit's weight
   in the reconstruction -- the summed ``|weight|`` of the terms it serves, times
   the 1-norm of the cuts it does not touch;
2. the rest in proportion to weight x the spread the pilot measured
   (``sqrt(1 - <Z>**2)`` of what the reconstruction reads from it): Neyman
   allocation, which minimises the variance for the budget.

The second round's task receives the pilot's counts and returns them merged with
its own. The master waits for the pilot to decide the second round.
"""

import math

import numpy as np

__all__ = ["unit_weight", "run_with_budget"]

#: The share of the budget the pilot round spends.
PILOT = 0.2
#: The fewest shots a circuit gets in the pilot round, when the budget allows.
MIN_PILOT_SHOTS = 32


def unit_weight(plan, circuit):
    """How much of the quasiprobability weight flows through one fragment circuit:
    the summed ``|weight|`` of its cuts' terms it serves, times the 1-norm of every
    cut it does not touch."""
    from Qdislib.core.cutting_algorithms.mixed_cutting import (
        _WIRE_DIGIT_BASIS,
        _WIRE_DIGIT_STATE,
    )
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _SETTING_PAULI,
        _WIRE_PREP,
    )

    fragment = circuit.template.fragment
    (digits, config), states, setting = circuit.spec
    weight, touched = 1.0, set()
    for k, side in fragment["gate_sides"]:
        ops, weights = plan.gate_tables[k - 1][side - 1], plan.gate_tables[k - 1][2]
        op = ops[digits[k - 1]]
        weight *= sum(abs(weights[d]) for d in range(6) if ops[d] == op)
        touched.add(("gate", k))
    for g in fragment["had_groups"]:
        settings, weights = plan._hadamard_settings[g]
        weight *= abs(weights[list(settings).index(tuple(config[g]))])
        touched.add(("hadamard", g))
    for i, state in zip(fragment["preps"], states):
        terms = plan.wire_digits[i - 1]
        if plan.optimized:
            weight *= sum(
                0.5 * abs(c) for d in terms for c, s in _WIRE_PREP[d] if s == state
            )
        else:
            weight *= sum(0.5 for d in terms if _WIRE_DIGIT_STATE[d] == state)
        touched.add(("wire", i))
    measures = fragment["measures"]
    if measures:
        if setting and setting[0] == "joint":  # a joint setting reads a share of all
            weight *= math.prod(
                0.5 * len(plan.wire_digits[i - 1]) for i in measures
            ) / (2 ** len(measures) + 1)
        else:
            for i, basis in zip(measures, setting):
                terms = plan.wire_digits[i - 1]
                if plan.optimized:
                    letters = ("I", "Z") if basis == "Z" else (basis,)
                    weight *= sum(0.5 for d in terms if _SETTING_PAULI[d] in letters)
                else:
                    weight *= sum(0.5 for d in terms if _WIRE_DIGIT_BASIS[d] == basis)
        touched |= {("wire", i) for i in measures}
    for f in range(len(plan.fragments)):
        for bond in plan.fragment_bonds(f):
            if bond not in touched:
                touched.add(bond)
                weight *= float(np.abs(plan.bond_weights(bond)).sum())
    return weight


def _split(budget, scores, floor=0):
    """Integer shots in proportion to ``scores``, at least ``floor`` each, summing
    to ``budget`` (or to ``floor * len`` when that is more)."""
    scores = np.maximum(np.asarray(scores, dtype=float), 1e-12)
    raw = np.maximum(budget * scores / scores.sum(), floor)
    shots = np.floor(raw).astype(np.int64)
    rest = int(budget - shots.sum())
    if rest > 0:  # hand out what flooring left, largest remainders first
        for i in np.argsort(-(raw - shots))[:rest]:
            shots[i] += 1
    return shots


def _spread(circuit, result):
    """``sqrt(mean(1 - v**2))`` over what the reconstruction reads from a circuit,
    from its pilot counts (a small floor so no circuit is starved)."""
    from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts
    from Qdislib.core.cutting_algorithms.mixed_cutting import _needed_masks

    corr, original = result
    if not isinstance(corr, _JointCounts):
        return 1.0
    reader = corr.correlators()
    masks = _needed_masks(circuit.template, original, circuit.spec[2], circuit.mask)
    if not masks:
        return 1.0
    shots = max(1, sum(corr.joint.values()))
    values = [reader[m] for m in masks]
    variance = float(np.mean([1 - v * v for v in values])) + 1.0 / shots
    return math.sqrt(max(variance, 1e-6))


def run_with_budget(groups, total_shots, seed):
    """Run every fragment circuit of every group within ``total_shots`` (see the
    module): the results of each group, counts kept, in its circuits' order.

    :param groups: ``[(plan, circuits, keys, scale), ...]``: ``keys`` the
        reconstruction's key of each circuit, ``scale`` the group's summed
        ``|coefficient|``.
    """
    from Qdislib.core.cutting_algorithms._reconstruction import next_seed
    from Qdislib.core.cutting_algorithms.mixed_cutting import _run_fragment
    from Qdislib.utils.runtime import compss_wait_on

    units = [
        (g, n, circuit)
        for g, (_, circuits, _, _) in enumerate(groups)
        for n, circuit in enumerate(circuits)
    ]
    # identical circuits (the same key and setting) are averaged by the
    # reconstruction: each carries its share of their weight
    copies = {}
    for g, (_, circuits, keys, _) in enumerate(groups):
        for circuit, key in zip(circuits, keys):
            slot = (g, key, repr(circuit.spec[2]))
            copies[slot] = copies.get(slot, 0) + 1
    if total_shots < len(units):
        raise ValueError(
            f"total_shots={total_shots} is fewer than the {len(units)} fragment "
            "circuits to run (at least one shot each)"
        )
    weights = [
        groups[g][3]
        * unit_weight(groups[g][0], c)
        / copies[(g, groups[g][2][n], repr(c.spec[2]))]
        for g, n, c in units
    ]
    pilot_budget = max(len(units), int(PILOT * total_shots))
    floor = min(MIN_PILOT_SHOTS, pilot_budget // len(units))
    pilot_shots = _split(pilot_budget, weights, floor)

    def run(g, n, circuit, shots, key, merge=None):
        return _run_fragment(
            circuit.template,
            circuit.spec,
            circuit.software,
            int(shots),
            circuit.qpu,
            circuit.qpu_dict,
            None if seed is None else next_seed(seed, key, g, n),
            None,
            None,
            True,
            merge,
        )

    pilot = compss_wait_on(
        [run(g, n, c, s, "pilot") for (g, n, c), s in zip(units, pilot_shots)]
    )
    rest = int(total_shots - pilot_shots.sum())
    second = (
        _split(
            rest,
            [w * _spread(c, r) for w, (_, _, c), r in zip(weights, units, pilot)],
        )
        if rest > 0
        else np.zeros(len(units), dtype=np.int64)
    )
    results = [[None] * len(circuits) for _, circuits, _, _ in groups]
    for (g, n, circuit), result, shots in zip(units, pilot, second):
        results[g][n] = (
            run(g, n, circuit, shots, "main", merge=result) if shots > 0 else result
        )
    return results
