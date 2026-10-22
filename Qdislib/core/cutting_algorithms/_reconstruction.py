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

"""Shared reconstruction helpers for the cutting algorithms.

Gate and wire cutting both reconstruct an expectation value as a *signed product*
of subcircuit results: each quasiprobability term contributes the product of its
subcircuit values, and the sign of that term is read off the base-``N`` digits of
the term index. The two schemes differ only in the base and which digit values
flip the sign -- gate cutting uses base 6 with ``{3, 5}``, wire cutting base 8
with ``{3, 5, 7}`` -- so the shared logic lives here and each module keeps a thin
PyCOMPSs ``@task`` wrapper around it.
"""

import math


def signed_product(values, index, base, flip_digits):
    """Signed product of subcircuit values for one quasiprobability term.

    :param values: Subcircuit expectation values for the term. Semantic-cache
        placeholders are resolved here (on the COMPSs worker, so the subcircuit
        futures feed straight in without a barrier); ``None`` entries are dropped.
    :param index: The term index; its digits in ``base`` determine the sign.
    :param base: Numeral base of the measure-and-prepare decomposition (6 for
        gate cutting, 8 for wire cutting).
    :param flip_digits: The digit values that flip the reconstruction sign.
    :return: The (possibly sign-flipped) product of the resolved values.
    """
    from Qdislib.utils.circuit_cache import _resolve_value

    resolved = [_resolve_value(x) for x in values if x is not None]
    product = math.prod(resolved)

    flip = False
    number = index
    while number != 0:
        if number % base in flip_digits:
            flip = not flip
        number //= base

    return -product if flip else product


def check_cuts_separate(cut_graph, sides):
    """Raise unless every cut separates: its two sides must land in different
    fragments.

    A cut whose two sides end up in the same fragment (another, uncut path still
    connects them) separates nothing: summing its terms inside that fragment just
    rebuilds what was cut, so it only multiplies the circuit count and the shots.
    It is not a valid cut, whatever the other cuts do -- even when they do split
    the circuit.

    :param cut_graph: The circuit after the cut surgery of every cut, red
        (ordering) edges removed, so its connected components are the fragments.
    :param sides: ``{cut: (node, node)}``, the two marker nodes of each cut.
    :raises InvalidCutError: Naming every cut that does not separate.
    """
    import networkx as nx

    from Qdislib.utils.exceptions import InvalidCutError

    fragment = {}
    for index, nodes in enumerate(nx.connected_components(cut_graph.to_undirected())):
        for node in nodes:
            fragment[node] = index
    invalid = [
        cut
        for cut, (a, b) in sides.items()
        if a in fragment and fragment[a] == fragment.get(b)
    ]
    if invalid:
        raise InvalidCutError(
            f"Not a valid cut: {invalid}. Each of these leaves its two sides in the "
            "same fragment -- another, uncut path still connects them -- so it "
            "separates nothing and only multiplies the cost. Every cut has to "
            "separate its two sides; drop these (or also cut the path that "
            "connects them)."
        )
