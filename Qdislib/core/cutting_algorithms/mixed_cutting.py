#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Mixed wire + gate cutting on a single circuit (circuit knitting).

Cut some two-qubit gates (the measure-and-prepare gate cut) *and* some wires (the
wire cut) on the same circuit, in one integrated, fully parallel workflow. The two
quasiprobability decompositions are independent, so they **tensor**: the combined QPD
runs over the product of the two configuration spaces and the reconstruction is

    <O> = 1 / 2^(n_gate + n_wire)
          * sum_{gate cfg} sum_{wire cfg}
                gate_sign(cfg) * wire_sign(cfg) * prod_fragment <fragment>

Every fragment across the whole product space is independent, so they all fan out as
PyCOMPSs tasks and fold back through one reduction (single fan-in, or a balanced tree
when ``reduce_chunk`` is set) -- the same parallel shape as the standalone cutters.

The heavy lifting is reused unchanged: the gate-cut DAG surgery + base-6
measure-and-prepare substitution (:func:`_generate_cut`, :func:`_apply_gate_substitution`),
the wire-cut DAG surgery + base-8 basis substitution (:func:`_insert_wire_cut_markers`,
:func:`_apply_wire_substitution`), and the fragment builder/executor
(:func:`_build_and_run_gate_subcircuit`, which also handles the gate cut's
collapse-and-continue measures a mixed fragment can carry).
"""

from __future__ import annotations

import functools
import itertools
import math

import numpy as np

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN
from Qdislib.utils.backends import Subcircuit
from Qdislib.core.cutting_algorithms._reconstruction import next_seed
from Qdislib.core.cutting_algorithms._api import (
    SDKS,
    standard_arguments,
    resolve_backend,
    resolve_software,
    default_observable,
    is_simulator,
)

__all__ = [
    "gate_wire_cutting",
    "circuit_cutting",
    "circuit_cutting_subcircuits",
    "circuit_cutting_sampling",
    "circuit_cutting_sampling_subcircuits",
    "cutting_overhead",
    "CutCircuits",
    "result_from_counts",
]


def _config_sign(index, base, flip_digits):
    """The +/-1 quasiprobability sign of configuration ``index``: flip once per digit
    (in base ``base``) that lands in ``flip_digits``. Gate cut uses base 6 / {3, 5};
    wire cut uses base 8 / {3, 5, 7} -- matching the standalone ``_change_sign`` rules.
    """
    flip = False
    n = index
    while n:
        if n % base in flip_digits:
            flip = not flip
        n //= base
    return -1 if flip else 1


@task(returns=1, values=COLLECTION_IN)
def _signed_product(values, sign):
    """One configuration's contribution: ``sign * prod(fragment values)`` (the fan-in
    task; ``values`` is the ``COLLECTION_IN`` of that config's fragment futures)."""
    product = 1.0
    for value in values:
        if value is not None:
            product *= value
    return sign * product


def gate_wire_cutting(
    circuit,
    gates_cut,
    wire_cut,
    observable=None,
    software="qibo",
    shots=1024,
    method="automatic",
    sync=True,
    qpu=False,
    qpu_dict=None,
    hadamard=False,
    optimized=False,
    joint=None,
):
    """Deprecated: use :func:`circuit_cutting` with
    ``cuts={"gate": gates_cut, "wire": wire_cut}``.

    Runs as it did in 1.0 (measure-and-prepare gate cuts with the standard wire
    cut); with ``hadamard``, ``optimized`` or ``joint`` it is forwarded to
    :func:`circuit_cutting`.
    """
    import warnings

    warnings.warn(
        'gate_wire_cutting() is deprecated: use circuit_cutting(circuit, {"gate": '
        '[...], "wire": [...]}, observable, ...); it is removed in the next release.',
        DeprecationWarning,
        stacklevel=2,
    )
    if hadamard or optimized or joint is not None:
        if qpu or qpu_dict:
            raise ValueError("circuit_cutting takes backend=, not qpu/qpu_dict")
        return circuit_cutting(
            circuit,
            {"gate": gates_cut or [], "wire": wire_cut or []},
            observable,
            backend=software,
            shots=shots,
            hadamard=hadamard,
            optimized=optimized,
            joint=True if joint is None else joint,
            sync=sync,
        )
    import networkx as nx

    from Qdislib.utils.backends import convert_to_dag, change_basis_dag
    from Qdislib.utils.graph_qibo import _max_qubits_graph
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _generate_cut,
        _apply_gate_substitution,
        _build_and_run_gate_subcircuit,
        _expand_general_gate_cuts,
    )
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _expand_wire_cuts,
        _insert_wire_cut_markers,
        _apply_wire_substitution,
        _remove_red_edges,
    )

    dag = convert_to_dag(circuit, observable)
    if observable:
        dag = change_basis_dag(dag, observable)
    num_qubits = _max_qubits_graph(dag)  # original qubit count (before wire additions)

    gates_cut = list(gates_cut or [])
    # A non-CZ/CX cut gate is decomposed into CZ + single-qubit gates and its CZs cut.
    if gates_cut:
        dag, gates_cut = _expand_general_gate_cuts(dag, gates_cut)

    # Wire surgery FIRST, then gate surgery. Doing wire cuts while the gate-cut gates
    # are still present lets a wire cut that sits right after a gate cut attach its
    # O_i/PS_i markers between the gate and its successor; the subsequent gate surgery
    # then wires the gate's SUBS output straight into that marker. (Gate-first would
    # have already replaced the gate, so an adjacent wire cut could not reference it.)
    wire_cuts = _expand_wire_cuts(dag, wire_cut) if wire_cut else []
    if wire_cuts:
        _insert_wire_cut_markers(dag, wire_cuts, num_qubits)
    # Gate surgery: insert the SUBS1/SUBS2 markers (also strips red edges internally).
    dag = _generate_cut(dag, gates_cut) if gates_cut else _remove_red_edges(dag)

    # Every cut, gate or wire, has to separate its two sides.
    from Qdislib.core.cutting_algorithms._reconstruction import check_cuts_separate

    sides = {g: (f"SUBS1_{k}", f"SUBS2_{k}") for k, g in enumerate(gates_cut, 1)}
    sides.update({w: (f"O_{i}", f"PS_{i}") for i, w in enumerate(wire_cuts, 1)})
    check_cuts_separate(_remove_red_edges(dag), sides)

    n_gate, n_wire = len(gates_cut), len(wire_cuts)
    total_qubits = num_qubits + n_wire  # fresh wire lines are numbered above the rest

    # Enumerate the PRODUCT of the two QPD configuration spaces; each configuration
    # applies both substitutions to a fresh DAG copy, splits into independent fragments,
    # fans them out, and contributes sign * prod(fragments).
    results = []
    for gate_idx in range(6**n_gate):
        for wire_idx in range(8**n_wire):
            config_dag = _remove_red_edges(dag)
            if n_gate:
                _apply_gate_substitution(config_dag, gates_cut, gate_idx)
            if n_wire:
                _apply_wire_substitution(config_dag, total_qubits, wire_idx, wire_cuts)
            config_dag = _remove_red_edges(config_dag)

            fragment_values = [
                _build_and_run_gate_subcircuit(
                    config_dag.subgraph(component).copy(),
                    total_qubits,
                    software,
                    shots,
                    method,
                    False,  # gpu
                    0,  # gpu_min_qubits
                    qpu,
                    qpu_dict,
                    None,  # batch
                    None,  # ibm_backend
                )
                for component in nx.connected_components(
                    config_dag.to_undirected(as_view=True)
                )
            ]
            sign = _config_sign(gate_idx, 6, {3, 5}) * _config_sign(
                wire_idx, 8, {3, 5, 7}
            )
            results.append(_signed_product(fragment_values, sign))

    if not sync:
        return results
    results = compss_wait_on(results)
    return 1 / (2 ** (n_gate + n_wire)) * sum(results)


# --------------------------------------------------------------------------- #
# circuit_cutting: any combination of cuts
# --------------------------------------------------------------------------- #
#
# Every cut is a quasiprobability decomposition, and independent decompositions
# tensor: the value is a sum over the product of their terms ("configurations") of
# weight * prod_fragment <fragment>. A fragment only sees the cuts that touch it, so
# its circuit depends on its own part of the configuration -- its "legs":
#
# - a measure-and-prepare gate cut: what each side marker becomes (5 circuits per
#   side, out of the 6 terms);
# - a Hadamard gate cut: the (s, t, g) of each group that has an ancilla here;
# - a wire it sends: the measurement basis (3 per wire, or 2**k + 1 jointly);
# - a wire it receives: the prepared state (4 base states; |-> and |-i> are linear
#   combinations of them).
#
# optimized=False runs every fragment of every configuration, as each standalone
# cutter does with its standard decomposition. optimized=True runs each fragment
# once per distinct set of legs and lets the reconstruction reuse it, so the cost is
# a sum over fragments instead of a product over cuts.
#
# Every fragment is read in the Z basis, which gives its correlators for every
# subset of its qubits in one run; the same runs give an expectation value or the
# whole output distribution (inverse Walsh-Hadamard transform).

_WIRE_DIGIT_BASIS = {0: "Z", 1: "Z", 2: "X", 3: "X", 4: "Y", 5: "Y", 6: "Z", 7: "Z"}
_WIRE_DIGIT_STATE = {0: "0", 1: "1", 2: "+", 3: "-", 4: "+i", 5: "-i", 6: "0", 7: "1"}
#: A wire's eight terms, and the six left without its two Y terms (digits 4, 5).
_ALL_WIRE_DIGITS = tuple(range(8))
_NO_Y_DIGITS = (0, 1, 2, 3, 6, 7)
_MARKER_PREFIXES = ("O_", "PS_", "GATE", "SUBS")


def _cz_table():
    from Qdislib.core.cutting_algorithms.gate_cutting import _GATE_SIDE_OPS

    weights = tuple((-1) ** (d in (3, 5)) / 2 for d in range(6))
    return (_GATE_SIDE_OPS[1], _GATE_SIDE_OPS[2], weights)


def _rzz_table(phi):
    """The six terms of the angle-aware measure-and-prepare cut of ``rzz(phi)``."""
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _RZZ_SIDE_OPS,
        _rzz_cut_weights,
    )

    weights = tuple(float(w) for w in _rzz_cut_weights(float(phi)))
    return (_RZZ_SIDE_OPS[1], _RZZ_SIDE_OPS[2], weights)


#: Gates real whatever their angle (a free parameter keeps them golden-eligible).
_REAL_AT_ANY_ANGLE = {"ry", "cry"}


def _is_real_gate(data):
    """True if a DAG node's gate is a real matrix (so it keeps a real state real)."""
    from Qdislib.utils.statevector import gate_unitary

    from Qdislib.core.cutting_algorithms.parameters import is_symbolic

    name = str(data.get("gate") or "")
    if name.upper() == "OBSERVABLE I":
        return True
    params = data.get("parameters") or ()
    if any(is_symbolic(p) for p in params):
        # a free angle: real only for the gates that are real at every angle
        return name.lower() in _REAL_AT_ANY_ANGLE
    matrix = gate_unitary(name, params)
    return matrix is not None and not np.any(np.abs(np.imag(matrix)) > 1e-12)


def _golden_wires(graph, fragments):
    """The wire cuts whose two Y terms are zero ("golden" cut points), so they are
    left out: 6 terms instead of 8, ``gamma = 3`` instead of 4. Exact.

    A wire qualifies when the fragment that sends it has a real state and reads
    only real operators: every gate a real matrix (an X/Y observable's basis change
    is a gate of the fragment, so a Y in the observable is caught), no gate-cut or
    Hadamard leg, and every wire it receives golden too (no ``|+i>`` prepared).
    Then ``<O (x) Y> = 0`` for every Z/I string ``O``: ``O (x) Y`` is imaginary.
    The fragment must send this wire only: with two, ``<O (x) Y (x) Y>`` is real
    and need not vanish.
    """
    real_gates = []
    for fragment in fragments:
        real_gates.append(
            not fragment["gate_sides"]
            and not fragment["had_groups"]
            and len(fragment["measures"]) == 1
            and all(
                _is_real_gate(graph.nodes[node])
                for node in fragment["nodes"]
                if not str(node).startswith(_MARKER_PREFIXES)
            )
        )
    golden = set()
    changed = True
    while changed:  # a fragment is real once every wire it receives is golden
        changed = False
        for fragment, real in zip(fragments, real_gates):
            if real and set(fragment["preps"]) <= golden:
                wire = fragment["measures"][0]
                if wire not in golden:
                    golden.add(wire)
                    changed = True
    return golden


from Qdislib.utils.exceptions import InvalidCutError  # noqa: E402


class _Draws(list):
    """The drawn configurations, ``[(digits, config, wire, weight, count)]``, and
    the stratum of each (``strata``): the error bars resample draws within it."""

    def __init__(self, draws, strata):
        super().__init__(draws)
        self.strata = list(strata)

    def __reduce__(self):
        return (_Draws, (list(self), self.strata))


def _normalise_cuts(cuts, hadamard):
    """``(gate, hadamard, wire)`` cut lists from a dict or a flat list.

    A dict takes the keys ``"gate"`` (measure-and-prepare), ``"hadamard"`` and
    ``"wire"`` -- :func:`~Qdislib.find_cut` with ``mixed=True`` returns one. In a list,
    a gate name (``str``) is a gate cut and a ``(source, target)`` pair a wire cut.
    ``hadamard=True`` cuts the ``"gate"`` / ``str`` entries with the Hadamard test.
    """
    if isinstance(cuts, dict):
        unknown = set(cuts) - {"gate", "hadamard", "wire"}
        if unknown:
            raise ValueError(
                f"Unknown cut kind(s) {sorted(unknown)}: use 'gate', 'hadamard', 'wire'."
            )
        for kind in ("gate", "hadamard", "wire"):
            if isinstance(cuts.get(kind), str):
                raise InvalidCutError(
                    f"cuts[{kind!r}] must be a list of cuts, got the string "
                    f"{cuts[kind]!r} (write [{cuts[kind]!r}])."
                )
        gate = list(cuts.get("gate") or [])
        had = list(cuts.get("hadamard") or [])
        wire = [tuple(edge) for edge in cuts.get("wire") or []]
    else:
        if isinstance(cuts, str):
            raise InvalidCutError(
                f"cuts must be a list of cuts, got the string {cuts!r} "
                f"(write [{cuts!r}])."
            )
        gate, had, wire = [], [], []
        for cut in cuts or []:
            if isinstance(cut, str):
                gate.append(cut)
            else:
                wire.append(tuple(cut))
    if hadamard:
        gate, had = [], had + gate
    return gate, had, wire


class _CutPlan:
    """A circuit after the cut surgery of every cut, split into its fragments.

    Picklable (a DAG and plain data), so the reconstruction can run as a task.
    """

    def __init__(
        self,
        circuit,
        cuts,
        observable,
        hadamard,
        optimized,
        joint,
        validate=True,
        golden=True,
    ):
        import networkx as nx

        from Qdislib.utils.backends import convert_to_dag, change_basis_dag
        from Qdislib.utils.graph_qibo import _remove_red_edges
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _generate_cut,
            _expand_rotation_cut,
        )
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            _generate_cut_hadamard,
            _cut_groups,
            _cut_settings,
            gamma_factor,
            local_unitary_decomposition,
            _num_circuit_qubits,
        )
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _expand_wire_cuts,
            _place_wire_cuts,
            _wire_cut_line_maps,
        )

        gate_cuts, had_cuts, wire_edges = _normalise_cuts(cuts, hadamard)
        #: The measure-and-prepare gate cuts as given (before a general gate is
        #: replaced by the CZs it decomposes into).
        self.requested_gate_cuts = list(gate_cuts)
        overlap = set(gate_cuts) & set(had_cuts)
        if overlap:
            raise ValueError(f"Gate(s) {sorted(overlap)} cut twice.")
        self.optimized = bool(optimized)
        self.joint_wires = joint_gates = bool(joint)

        dag = convert_to_dag(circuit)
        from Qdislib.core.cutting_algorithms.parameters import (
            check_numeric,
            circuit_parameters,
        )

        #: The free parameters (names), in the order a list of values follows.
        self.parameters = circuit_parameters(circuit, dag)
        if self.parameters:
            check_numeric(dag, list(gate_cuts) + list(had_cuts))
        if observable and set(observable) & {"X", "Y"}:
            dag = change_basis_dag(dag, observable)
        self.num_qubits = 1 + max(
            (q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()),
            default=0,
        )
        self.width = max(self.num_qubits, _num_circuit_qubits(circuit))

        # Read the Hadamard decompositions before any surgery touches the gates.
        decomps = [
            local_unitary_decomposition(
                (dag.nodes[name].get("gate") or "").upper(),
                dag.nodes[name].get("parameters") or (),
            )
            for name in had_cuts
        ]
        #: A measure-and-prepare cut of a gate other than CZ/CX -> ``(its gamma, the
        #: Hadamard cut's gamma for it)``, to point out a cheaper cut.
        self.general_gate_cuts = {}
        #: Per measure-and-prepare cut: ``(side 1 ops, side 2 ops, weights)`` of its
        #: six terms -- the CZ table (+-1/2), or an ``rzz`` rotation's, at its angle.
        self.gate_tables = []
        label = {}  # a cut node -> how a message names it (the gate the user cut)
        if gate_cuts:
            expanded = []
            for name in gate_cuts:
                data = dag.nodes[name]
                kind = (data.get("gate") or "").lower()
                if kind in ("cz", "cx", "cnot"):
                    expanded.append(name)
                    self.gate_tables.append(_cz_table())
                    continue
                hadamard_gamma = None
                try:
                    hadamard_gamma = gamma_factor(
                        local_unitary_decomposition(
                            kind.upper(), data.get("parameters") or ()
                        )
                    )
                except Exception:  # noqa: BLE001 -- no Hadamard form known
                    pass
                # Cut as its rzz rotations, each at its own angle (not as CZs).
                rotations = _expand_rotation_cut(dag, name)
                for i, r in enumerate(rotations, 1):
                    label[r] = name if len(rotations) == 1 else f"{name} (rotation {i})"
                tables = [_rzz_table(dag.nodes[r]["parameters"][0]) for r in rotations]
                expanded.extend(rotations)
                self.gate_tables.extend(tables)
                self.general_gate_cuts[name] = (
                    math.prod(float(np.abs(t[2]).sum()) for t in tables),
                    hadamard_gamma,
                )
            gate_cuts = expanded
        # Wires first, as gate_wire_cutting does: a wire cut right after a cut gate
        # still finds the gate in place.
        wire_cuts = _expand_wire_cuts(dag, wire_edges) if wire_edges else []
        placed = _place_wire_cuts(dag, wire_cuts, self.num_qubits)
        self.cut_source_index, self.prepared_to_orig = _wire_cut_line_maps(
            placed, self.num_qubits
        )
        #: How the measure-and-prepare cuts are named to the user (the gate cut).
        self.gate_labels = [label.get(g, g) for g in gate_cuts]
        if gate_cuts:
            dag = _generate_cut(dag, gate_cuts, rotations=True)
            # Both gate cuts name their markers SUBS{side}_{k}; keep them apart.
            dag = nx.relabel_nodes(
                dag,
                {
                    f"SUBS{side}_{k}": f"GATE{side}_{k}"
                    for k in range(1, len(gate_cuts) + 1)
                    for side in (1, 2)
                    if f"SUBS{side}_{k}" in dag
                },
            )
        if had_cuts:
            dag = _generate_cut_hadamard(dag, had_cuts)
        graph = _remove_red_edges(dag)

        from Qdislib.core.cutting_algorithms._reconstruction import (
            check_cuts_separate,
        )

        sides = {
            label.get(g, g): (f"GATE1_{k}", f"GATE2_{k}")
            for k, g in enumerate(gate_cuts, 1)
        }
        sides.update(
            {g: (f"SUBS1_{k}", f"SUBS2_{k}") for k, g in enumerate(had_cuts, 1)}
        )
        sides.update({w: (f"O_{i}", f"PS_{i}") for i, w in enumerate(wire_cuts, 1)})
        #: Each cut -> its two marker nodes in ``graph``.
        self.sides = sides
        if validate:
            check_cuts_separate(graph, sides)
        self.graph = graph
        self.gate_cuts, self.had_cuts, self.wire_cuts = gate_cuts, had_cuts, wire_cuts
        self.ancilla_base = self.num_qubits + len(wire_cuts) + 1

        components = list(nx.connected_components(graph.to_undirected(as_view=True)))
        self.groups = (
            _cut_groups(graph, components, decomps, joint_gates) if had_cuts else []
        )
        group_of = {
            k: index
            for index, (cut_indices, _) in enumerate(self.groups)
            for k in cut_indices
        }
        self.fragments = []
        for nodes in components:
            w = range(1, len(wire_cuts) + 1)
            measures = [i for i in w if f"O_{i}" in nodes]
            had_nodes = sorted(
                (node for node in nodes if str(node).startswith("SUBS")),
                key=lambda node: (int(node.split("_")[1]), node),
            )
            had_groups = sorted({group_of[int(n.split("_")[1])] for n in had_nodes})
            keys = []  # the fragment's ancillas, one per key
            for node in had_nodes:
                group = group_of[int(node.split("_")[1])]
                key = (
                    ("group", group)
                    if len(self.groups[group][0]) > 1
                    else ("node", node)
                )
                if key not in keys:
                    keys.append(key)
            self.fragments.append(
                {
                    "nodes": frozenset(nodes),
                    "measures": measures,
                    "cut_qubits": {
                        i: graph.nodes[f"O_{i}"]["qubits"][0] for i in measures
                    },
                    "preps": [i for i in w if f"PS_{i}" in nodes],
                    "gate_sides": [
                        (k, side)
                        for k in range(1, len(gate_cuts) + 1)
                        for side in (1, 2)
                        if f"GATE{side}_{k}" in nodes
                    ],
                    "had_groups": had_groups,
                    "ancilla_of": {
                        key: self.ancilla_base + j for j, key in enumerate(keys)
                    },
                }
            )
        golden_wires = _golden_wires(graph, self.fragments) if golden else set()
        #: Per wire cut, the terms (base-8 digits) it keeps: all eight, or six
        #: when its Y terms are provably zero (see :func:`_golden_wires`).
        self.wire_digits = [
            _NO_Y_DIGITS if i in golden_wires else _ALL_WIRE_DIGITS
            for i in range(1, len(wire_cuts) + 1)
        ]

    # The configurations are built on first use: a plan made only to inspect or
    # draw a cut never enumerates them (there are 6**g * 8**w of them).

    @functools.cached_property
    def gate_configs(self):
        """``[(mp digits, Hadamard config, weight), ...]``: measure-and-prepare
        terms x Hadamard settings."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import _cut_settings

        mp = [
            (
                digits,
                math.prod(
                    float(table[2][d]) for table, d in zip(self.gate_tables, digits)
                ),
            )
            for digits in itertools.product(range(6), repeat=len(self.gate_cuts))
        ]
        had = []
        for config in itertools.product(
            *(_cut_settings(len(terms)) for _, terms in self.groups)
        ):
            weight = 1.0
            for (_, terms), (s, t, g) in zip(self.groups, config):
                weight *= terms[s][0] * terms[t][0] * (-1) ** g
            if weight != 0.0:
                had.append((config, weight))
        return [
            (digits, config, w_mp * w_had)
            for digits, w_mp in mp
            for config, w_had in had
        ]

    @functools.cached_property
    def wire_configs(self):
        """``[(base-8 digits, weight), ...]`` over the wire cuts."""
        return [
            (
                digits,
                (-1) ** sum(d in (3, 5, 7) for d in digits) / 2 ** len(self.wire_cuts),
            )
            for digits in itertools.product(*self.wire_digits)
        ]

    def invalid_cuts(self):
        """The cuts whose two sides share a fragment (see ``check_cuts_separate``)."""
        import networkx as nx

        fragment = {}
        for index, nodes in enumerate(
            nx.connected_components(self.graph.to_undirected(as_view=True))
        ):
            for node in nodes:
                fragment[node] = index
        # A gate cut as several rotations is named once, as the user named it.
        return list(
            dict.fromkeys(
                cut.split(" (rotation")[0] if isinstance(cut, str) else cut
                for cut, (a, b) in self.sides.items()
                if a in fragment and fragment[a] == fragment.get(b)
            )
        )

    @property
    def num_configurations(self):
        """``6**g * (Hadamard settings) * 8**w`` (6 per golden wire), without
        listing them."""
        return (
            6 ** len(self.gate_cuts)
            * math.prod(self._group_settings())
            * math.prod(len(digits) for digits in self.wire_digits)
        )

    def _group_settings(self):
        """How many settings of each Hadamard group carry a non-zero weight, in
        closed form (a joint group's settings are never listed)."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import _num_settings

        return [_num_settings(terms) for _, terms in self.groups]

    def num_units(self):
        """How many fragment circuits :meth:`units` returns, without listing them.

        Standard: every fragment of every configuration. Optimized: each fragment
        runs once per setting of its own cuts -- 5 per side of a measure-and-prepare
        gate (6 when both sides are in it), the settings of each Hadamard group it
        holds an ancilla of, 4 per wire it receives, and 3 per wire it sends
        (``2**k + 1`` for ``k`` sent jointly).
        """
        if not self.optimized:
            return self.num_configurations * len(self.fragments)
        group_settings = self._group_settings()
        total = 0
        for fragment in self.fragments:
            sides = {}
            for k, _ in fragment["gate_sides"]:
                sides[k] = sides.get(k, 0) + 1
            count = math.prod(6 if n == 2 else 5 for n in sides.values())
            count *= math.prod(group_settings[g] for g in fragment["had_groups"])
            count *= math.prod(len(self._prep_states(i)) for i in fragment["preps"])
            count *= len(self._measure_settings(fragment))
            total += count
        return total

    def gamma(self):
        """The sampling overhead: the product of every cut's 1-norm."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import phi_factor

        total = math.prod(float(np.abs(t[2]).sum()) for t in self.gate_tables)
        total *= math.prod(len(digits) / 2 for digits in self.wire_digits)
        for _, terms in self.groups:
            total *= phi_factor(terms)
        return total

    # The reconstruction is a tensor network: each measure-and-prepare cut, Hadamard
    # group and wire is a bond (an index over its terms, with its weights), each
    # fragment a tensor with one axis per bond it touches. Contracting it costs about
    # a sum over fragments, where listing every configuration costs the product.

    @functools.cached_property
    def _hadamard_settings(self):
        """Per Hadamard group: its settings ``(s, t, g)`` with a non-zero weight, and
        those weights."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import _cut_settings

        out = []
        for _, terms in self.groups:
            settings, weights = [], []
            for s, t, g in _cut_settings(len(terms)):
                weight = terms[s][0] * terms[t][0] * (-1) ** g
                if weight != 0.0:
                    settings.append((s, t, g))
                    weights.append(weight)
            out.append((settings, np.array(weights)))
        return out

    def fragment_bonds(self, f):
        """The bonds fragment ``f`` touches: ``("gate", k)``, ``("hadamard", group)``
        and ``("wire", i)``."""
        fragment = self.fragments[f]
        bonds = [("gate", k) for k, _ in fragment["gate_sides"]]
        bonds += [("hadamard", g) for g in fragment["had_groups"]]
        bonds += [("wire", i) for i in fragment["measures"] + fragment["preps"]]
        return list(dict.fromkeys(bonds))

    def bond_weights(self, bond):
        """A bond's term weights (measure-and-prepare: 6 terms of +-1/2; wire: 8;
        Hadamard: ``c_s c_t (-1)**g`` per setting)."""
        kind, index = bond
        if kind == "gate":
            return np.asarray(self.gate_tables[index - 1][2], dtype=float)
        if kind == "wire":
            return np.array(
                [(-1) ** (d in (3, 5, 7)) / 2 for d in self.wire_digits[index - 1]]
            )
        return self._hadamard_settings[index][1]

    def gate_legs(self, f, digits, config):
        """Fragment ``f``'s part of a gate configuration."""
        fragment = self.fragments[f]
        return (
            tuple(
                self.gate_tables[k - 1][s - 1][digits[k - 1]]
                for k, s in fragment["gate_sides"]
            ),
            tuple(config[g] for g in fragment["had_groups"]),
        )

    def draw(self, num_samples, seed=None, stratified=True):
        """Draw ``num_samples`` configurations with probability ``|weight| / gamma``
        (Monte-Carlo quasiprobability sampling), without listing them all.

        Measure-and-prepare and wire terms all weigh the same, so their digits are
        uniform; a Hadamard group's setting ``(s, t, g)`` is drawn with probability
        ``c_s c_t / phi`` (Harrow & Lowe), ``s`` and ``t`` gate by gate.

        ``stratified`` (the default) allocates the draws instead of drawing them
        independently (Dai & Koczor: allocating shots to the terms in proportion to
        their weight is never worse than drawing them at random). Every cut is
        stratified when all configurations fit in ``num_samples``; otherwise the
        cuts are taken smallest first while there are about 8 draws or more per
        joint configuration. Each of those configurations gets ``floor(N p)`` draws, and
        the remaining draws are spread by systematic sampling, so the total stays
        ``N`` and each configuration is drawn ``N p`` times on average (unbiased).
        Within it, the other cuts are drawn independently. With equal weights and
        ``N`` a multiple of the configuration count it is exactly every
        configuration, ``N / count`` times each.

        :return: ``[(digits, config, wire, weight, count), ...]``, one per distinct
            draw: ``weight`` per draw (``gamma * sign / num_samples``, or
            ``gamma * p * sign / n_h`` in a stratum of probability ``p`` drawn
            ``n_h`` times), ``count`` the times it was drawn.
        """
        from Qdislib.core.cutting_algorithms.hadamard_cut import _draw_settings

        rng = np.random.default_rng(seed)
        n = int(num_samples)
        gamma = self.gamma()
        digits = rng.integers(0, 6, size=(n, len(self.gate_cuts)))
        gate_probs = [np.abs(t[2]) / np.abs(t[2]).sum() for t in self.gate_tables]
        for k, q in enumerate(gate_probs):
            if not np.allclose(q, q[0]):  # a rotation: its terms by weight
                digits[:, k] = rng.choice(6, size=n, p=q)
        wires = rng.integers(0, 8, size=(n, len(self.wire_cuts)))
        for i, choices in enumerate(self.wire_digits):
            if len(choices) < 8:  # a golden wire: its six terms
                wires[:, i] = np.array(choices)[rng.integers(0, len(choices), n)]
        group_counts = self._group_settings()
        fixed = {}  # group -> (s, t, g) of the groups the strata set
        per_draw = np.full(n, gamma / n)  # the weight of each draw (before its sign)
        # each draw's stratum, for the error bars: draws resample within it
        group = np.zeros(n, dtype=np.int64)
        if stratified:
            # The cuts as dimensions with their term probabilities, smallest first.
            dims = [("gate", k, 6) for k in range(len(self.gate_cuts))]
            dims += [
                ("wire", i, len(digits)) for i, digits in enumerate(self.wire_digits)
            ]
            dims += [("group", g, c) for g, c in enumerate(group_counts)]
            dims.sort(key=lambda d: d[2])
            # Every cut when all configurations fit in N (that is running them all,
            # in proportion); otherwise only while each stratum gets about 8 draws
            # or more, where splitting N unevenly (n_h and n_h + 1) costs < 2%.
            limit = n if math.prod(d[2] for d in dims) <= n else n / 8
            chosen, size = [], 1
            for dim in dims:
                if size * dim[2] > limit:
                    break
                chosen.append(dim)
                size *= dim[2]
            probs = []
            for kind, index, count in chosen:
                if kind == "group":
                    weights = np.abs(self._hadamard_settings[index][1])
                    probs.append(weights / weights.sum())
                elif kind == "gate":
                    probs.append(gate_probs[index])
                else:
                    probs.append(np.full(count, 1.0 / count))
            # Every joint configuration of the chosen cuts, its probability, and
            # its draws: floor(N p) each, the rest by systematic sampling.
            p = np.ones(1)
            for q in probs:
                p = np.outer(p, q).reshape(-1)
            share = n * p
            counts = np.floor(share).astype(np.int64)
            remainder = share - counts
            extra = n - int(counts.sum())
            if extra > 0:
                points = rng.random() + np.arange(extra)
                cuts_at = np.cumsum(remainder)
                where = np.searchsorted(cuts_at, points, side="right")
                counts += np.bincount(
                    np.minimum(where, len(counts) - 1), minlength=len(counts)
                )
            stratum = np.repeat(np.arange(len(counts)), counts)
            rng.shuffle(stratum)
            # A stratum always drawn (N p >= 1) is weighted by what it actually got,
            # gamma p / n_h per draw: the stratified estimator, whose variance does
            # not take on the rounding of n_h. A rarer one keeps gamma / N
            # (Horvitz-Thompson); both are unbiased.
            sure = share >= 1.0
            per_stratum = np.where(sure, gamma * p / np.maximum(counts, 1), gamma / n)
            per_draw = per_stratum[stratum]
            # a sure stratum resamples on its own; the rare ones (Horvitz-Thompson)
            # together, as one group
            group = np.where(sure[stratum], stratum + 1, 0)
            # Decode each draw's stratum into its chosen cuts' terms.
            radix = [len(q) for q in probs]
            place = stratum.copy()
            values = []
            for r in reversed(radix):
                values.append(place % r)
                place //= r
            values.reverse()
            for (kind, index, _), value in zip(chosen, values):
                if kind == "gate":
                    digits[:, index] = value
                elif kind == "wire":
                    wires[:, index] = np.array(self.wire_digits[index])[value]
                else:
                    picked = np.array(self._hadamard_settings[index][0])[value]
                    fixed[index] = (picked[:, 0], picked[:, 1], picked[:, 2])
        flips = np.isin(wires, (3, 5, 7)).sum(axis=1)
        for k, table in enumerate(self.gate_tables):
            flips = flips + (np.asarray(table[2])[digits[:, k]] < 0)
        settings = []
        for index, (_, terms) in enumerate(self.groups):
            if index in fixed:
                s, t, g = fixed[index]
            else:
                s, t, g = _draw_settings(terms, n, rng)
            flips = flips + g
            settings.append((s, t, g))
        signs = np.where(flips % 2, -1.0, 1.0)
        distinct = {}
        group_of = {}
        for row in range(n):
            config = tuple(
                (int(s[row]), int(t[row]), int(g[row])) for s, t, g in settings
            )
            key = (
                tuple(int(d) for d in digits[row]),
                config,
                tuple(int(w) for w in wires[row]),
            )
            entry = distinct.setdefault(key, [signs[row], per_draw[row], 0])
            group_of.setdefault(key, int(group[row]))
            entry[2] += 1
        return _Draws(
            [
                (d, c, w, sign * weight, count)
                for (d, c, w), (sign, weight, count) in distinct.items()
            ],
            [group_of[key] for key in distinct],
        )

    def _setting_for(self, f, letters):
        """The measurement setting of fragment ``f`` that reads the Pauli
        ``letters`` on the wires it sends."""
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _joint_measurement_settings,
        )

        m = len(letters)
        if self.joint_wires and m > 1:
            for t, (_, readout) in enumerate(_joint_measurement_settings(m)):
                if any(tuple(p) == letters for p, _, _ in readout):
                    return ("joint", t)
        return tuple("Z" if letter in "IZ" else letter for letter in letters)

    def units(self, draws=None):
        """``[(key, f, gate_config, prep_states, setting, count), ...]``: the fragment
        circuits to run. ``key`` (``f`` first) is how the reconstruction finds them;
        ``count`` how many draws use the circuit (1 when every configuration is run)."""
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _SETTING_PAULI,
            _WIRE_PREP,
        )

        units = []
        if draws is not None:
            if not self.optimized:
                # The fragment circuits the draws use, one per distinct set of the
                # fragment's own terms, with the shots of every draw that uses it
                # (pooled, as running every configuration pools them).
                needed = {}
                for digits, config, wire, _, count in draws:
                    for f, fragment in enumerate(self.fragments):
                        setting = tuple(
                            _WIRE_DIGIT_BASIS[wire[i - 1]] for i in fragment["measures"]
                        )
                        states = tuple(
                            _WIRE_DIGIT_STATE[wire[i - 1]] for i in fragment["preps"]
                        )
                        key = (f, self.gate_legs(f, digits, config), states)
                        entry = needed.setdefault(
                            (key, setting), [(digits, config), states, 0]
                        )
                        entry[2] += count
                for (key, setting), (gate_config, states, count) in needed.items():
                    units.append((key, key[0], gate_config, states, setting, count))
                return units
            needed = {}  # the circuits the draws read, and how many draws read each
            for digits, config, wire, _, count in draws:
                for f, fragment in enumerate(self.fragments):
                    leg = self.gate_legs(f, digits, config)
                    letters = tuple(
                        _SETTING_PAULI[wire[i - 1]] for i in fragment["measures"]
                    )
                    setting = self._setting_for(f, letters)
                    for choice in itertools.product(
                        *(_WIRE_PREP[wire[i - 1]] for i in fragment["preps"])
                    ):
                        states = tuple(state for _, state in choice)
                        entry = needed.setdefault(
                            ((f, leg, states), setting), [(digits, config), 0]
                        )
                        entry[1] += count
            for (key, setting), (gate_config, count) in needed.items():
                f, _, states = key
                units.append((key, f, gate_config, states, setting, count))
            return units
        if not self.optimized:
            # Every fragment of every configuration; the circuits with the same legs
            # share a key and are averaged in the reconstruction.
            for (digits, config, _), (wire, _) in itertools.product(
                self.gate_configs, self.wire_configs
            ):
                for f, fragment in enumerate(self.fragments):
                    setting = tuple(
                        _WIRE_DIGIT_BASIS[wire[i - 1]] for i in fragment["measures"]
                    )
                    states = tuple(
                        _WIRE_DIGIT_STATE[wire[i - 1]] for i in fragment["preps"]
                    )
                    key = (f, self.gate_legs(f, digits, config), states)
                    units.append((key, f, (digits, config), states, setting, 1))
            return units
        for f, fragment in enumerate(self.fragments):
            for gate_config in self._local_gate_configs(f):
                for states in itertools.product(
                    *(self._prep_states(i) for i in fragment["preps"])
                ):
                    for setting in self._measure_settings(fragment):
                        leg = self.gate_legs(f, *gate_config)
                        units.append(
                            ((f, leg, states), f, gate_config, states, setting, 1)
                        )
        return units

    def _prep_states(self, wire):
        """The states the optimized decomposition prepares for a wire it receives
        (no ``|+i>`` for a golden wire)."""
        from Qdislib.core.cutting_algorithms.wire_cutting import _PREP_STATES

        if len(self.wire_digits[wire - 1]) < 8:
            return tuple(state for state in _PREP_STATES if state != "+i")
        return _PREP_STATES

    def _measure_settings(self, fragment):
        """The optimized measurement settings of the wires a fragment sends (no Y
        basis for a golden wire, which is always sent alone)."""
        from Qdislib.core.cutting_algorithms.wire_cutting import _measurement_settings

        measures = fragment["measures"]
        if len(measures) == 1 and len(self.wire_digits[measures[0] - 1]) < 8:
            return [("X",), ("Z",)]
        return _measurement_settings(len(measures), self.joint_wires)

    def _local_gate_configs(self, f):
        """One gate configuration per distinct set of fragment ``f``'s gate legs,
        from its own cuts only (the other entries are 0 / ``None``)."""
        fragment = self.fragments[f]
        local = list(dict.fromkeys(k for k, _ in fragment["gate_sides"]))
        settings = self._hadamard_settings
        seen = {}
        for local_digits in itertools.product(range(6), repeat=len(local)):
            digits = [0] * len(self.gate_cuts)
            for k, d in zip(local, local_digits):
                digits[k - 1] = d
            for chosen in itertools.product(
                *(settings[g][0] for g in fragment["had_groups"])
            ):
                config = [None] * len(self.groups)
                for g, setting in zip(fragment["had_groups"], chosen):
                    config[g] = setting
                gate_config = (tuple(digits), tuple(config))
                seen.setdefault(self.gate_legs(f, *gate_config), gate_config)
        return list(seen.values())

    def build(self, f, gate_config, states, setting):
        """The DAG of fragment ``f`` for one set of legs, and its readout."""
        return _FragmentTemplate(self, f).build(gate_config, states, setting)


class _FragmentTemplate:
    """Fragment ``f`` of a :class:`_CutPlan`, as its circuits are built from: its
    part of the cut DAG and its cut markers.

    All the fragment's circuits share one, so PyCOMPSs ships it to a worker once and
    each circuit is built there from a small spec (:func:`_run_fragment`).
    """

    def __init__(self, plan, f):
        self.f = f
        self.fragment = plan.fragments[f]
        self.graph = plan.graph.subgraph(self.fragment["nodes"]).copy()
        self.groups = plan.groups
        self.gate_tables = plan.gate_tables
        # names of the cuts, for the readable views (not part of the circuit)
        self.gate_labels = plan.gate_labels
        self.had_cuts = plan.had_cuts
        self.wire_cuts = plan.wire_cuts
        self.ancilla_base = plan.ancilla_base
        self.cut_source_index = plan.cut_source_index
        self.prepared_to_orig = plan.prepared_to_orig
        self.parameters = plan.parameters
        qubits = {
            q for _, d in self.graph.nodes(data=True) for q in d.get("qubits") or ()
        }
        #: The width of its circuits: its qubits and its Hadamard ancillas.
        self.num_qubits = len(qubits) + len(self.fragment["ancilla_of"])

    #: ``{name: float}`` bound into its circuits (:meth:`bound`); ``None``: as cut.
    values = None

    def bound(self, values):
        """This template with ``values`` bound into the circuits it builds (the
        graph is shared, not copied)."""
        import copy

        template = copy.copy(self)
        template.__dict__.pop("_fingerprint", None)
        template.__dict__.pop("_bound_graph", None)
        template.values = values
        #: the unbound template, what PyCOMPSs ships (once, being the same object
        #: for every set of values); the values travel apart (:func:`_ship`)
        template._unbound = getattr(self, "_unbound", self)
        return template

    def _graph(self):
        """Its graph with :attr:`values` bound, bound once per process (every
        circuit of the fragment starts from it)."""
        if not self.values:
            return self.graph
        graph = self.__dict__.get("_bound_graph")
        if graph is None:
            from Qdislib.core.cutting_algorithms.parameters import bind_dag

            graph = self.graph.copy()
            bind_dag(graph, self.values)
            self._bound_graph = graph
        return graph

    def __getstate__(self):
        state = dict(self.__dict__)
        state.pop("_bound_graph", None)  # rebuilt where it is needed, not shipped
        return state

    def leg(self, gate_config):
        """This fragment's own part of a gate configuration (as
        :meth:`_CutPlan.gate_legs`): what its circuit depends on."""
        digits, config = gate_config
        return (
            tuple(
                self.gate_tables[k - 1][s - 1][digits[k - 1]]
                for k, s in self.fragment["gate_sides"]
            ),
            tuple(config[g] for g in self.fragment["had_groups"]),
        )

    def build(self, gate_config, states, setting):
        """The fragment's DAG for one set of legs, and its readout."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _apply_gate_side,
        )
        from Qdislib.core.cutting_algorithms.hadamard_cut import _config_assignments
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _apply_measurement_setting,
        )

        fragment = self.fragment
        digits, config = gate_config
        sub = self._graph().copy()
        for k, side in fragment["gate_sides"]:
            node = f"GATE{side}_{k}"
            _apply_gate_side(
                sub,
                node,
                self.gate_tables[k - 1][side - 1][digits[k - 1]],
                f"{node}_rz",
            )
        for cut, state in zip(fragment["preps"], states):
            _prepare_wire_state(sub, cut, state)
        readout = _apply_measurement_setting(
            sub, fragment["measures"], fragment["cut_qubits"], setting
        )
        if fragment["had_groups"]:
            assignments = {}
            for g in fragment["had_groups"]:
                # one group at a time, under its own index (the ancilla key)
                groups = [((), [(1.0, ())])] * g + [self.groups[g]]
                configs = [(0, 0, 0)] * g + [config[g]]
                assignments.update(_config_assignments(groups, configs))
            _insert_hadamard_gadgets(
                sub,
                {n: a for n, a in assignments.items() if n in sub},
                fragment["ancilla_of"],
            )
        return sub, readout


def _needed_masks(template, original, setting, mask):
    if isinstance(mask, (tuple, list)):  # the terms of a commuting group
        return sorted(
            {m for one in mask for m in _needed_masks(template, original, setting, one)}
        )
    return _one_term_masks(template, original, setting, mask)


def _one_term_masks(template, original, setting, mask):
    """The correlators an expectation value reads from one circuit: the
    observable's output bits, every ancilla bit, and each readout's cut-wire bits
    -- a few masks instead of ``2**width``, so a big fragment's task computes and
    sends back almost nothing."""
    fragment = template.fragment
    base = 0
    for i, qubit in enumerate(original):
        if qubit >= template.ancilla_base:
            base |= 1 << i
        elif template.cut_source_index.get(qubit) not in fragment["measures"]:
            output = template.prepared_to_orig.get(qubit, qubit)
            base |= ((mask >> output) & 1) << i
    positions = [
        original.index(fragment["cut_qubits"][c]) for c in fragment["measures"]
    ]
    return sorted(
        {
            base | sum(1 << positions[b] for b in bits)
            for _, _, bits in _setting_readout(len(fragment["measures"]), setting)
        }
    )


def _setting_readout(num_measures, setting):
    """The readout of a measurement setting (what :func:`_apply_measurement_setting`
    returns), without building the circuit."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _joint_measurement_settings,
        _per_wire_readout,
    )

    if setting and setting[0] == "joint":
        return _joint_measurement_settings(num_measures)[setting[1]][1]
    return _per_wire_readout(setting)


def _prepare_wire_state(graph, cut, state):
    """Set the ``PS_cut`` marker to prepare ``state``: a base state, or ``-`` / ``-i``
    (an X in front of ``+`` / ``+i``), for the standard decomposition."""
    from Qdislib.core.cutting_algorithms.wire_cutting import _apply_wire_prep_state

    if state not in ("-", "-i"):
        _apply_wire_prep_state(graph, cut, state)
        return
    _apply_wire_prep_state(graph, cut, "+i" if state == "-i" else "+")
    first = f"PS_{cut}_Hi" if state == "-i" else f"PS_{cut}"
    flip = f"PS_{cut}_X"
    graph.add_node(flip, gate="x", qubits=graph.nodes[first]["qubits"], parameters=())
    graph.add_edge(flip, first, color="blue")


def _insert_hadamard_gadgets(dag, assignments, ancilla_of):
    """Replace each Hadamard-cut marker in a fragment DAG by its select gate, in place.

    ``assignments`` maps a marker to ``(u_i, u_j, g, key)``; markers with the same
    ``key`` share the ancilla ``ancilla_of[key]``, opened (H, then S if ``g``) before
    the first and closed (H) after the last -- the DAG form of
    :func:`~Qdislib.core.cutting_algorithms.hadamard_cut._component_hadamard_circuit`,
    so the fragment can also hold the other cuts' markers and mid-circuit measurements.
    """
    import networkx as nx

    from Qdislib.core.cutting_algorithms.hadamard_cut import _u3_angles

    order = [node for node in nx.topological_sort(dag) if node in assignments]
    last_marker = {assignments[node][3]: node for node in order}
    tail = {}  # ancilla -> its last node

    def append(ancilla, name, gate, params=()):
        dag.add_node(name, gate=gate, qubits=(ancilla,), parameters=params)
        if tail.get(ancilla) is not None:
            dag.add_edge(tail[ancilla], name, color="blue")
        tail[ancilla] = name

    for node in order:
        u_i, u_j, g, key = assignments[node]
        ancilla = ancilla_of[key]
        (qubit,) = dag.nodes[node]["qubits"]
        if ancilla not in tail:
            tail[ancilla] = None
            append(ancilla, f"{node}_open", "h")
            if g == 1:
                append(ancilla, f"{node}_open_S", "s")
        theta_i, phi_i, lam_i, _ = _u3_angles(u_i)
        theta_m, phi_m, lam_m, gphase = _u3_angles(
            np.asarray(u_j) @ np.asarray(u_i).conj().T
        )
        preds, succs = list(dag.pred[node]), list(dag.succ[node])
        dag.remove_node(node)
        u3, cu3 = f"{node}_U3", f"{node}_CU3"
        dag.add_node(u3, gate="u", qubits=(qubit,), parameters=(theta_i, phi_i, lam_i))
        dag.add_node(
            cu3,
            gate="cu",  # OpenQASM's cu, phased to control qibo's SU(2) U3 exactly
            qubits=(ancilla, qubit),
            parameters=(theta_m, phi_m, lam_m, -(phi_m + lam_m) / 2),
        )
        for pred in preds:
            dag.add_edge(pred, u3, color="blue")
        dag.add_edge(u3, cu3, color="blue")
        for succ in succs:
            dag.add_edge(cu3, succ, color="blue")
        if tail[ancilla] is not None:
            dag.add_edge(tail[ancilla], cu3, color="blue")
        tail[ancilla] = cu3
        if abs(gphase) > 1e-12:
            append(ancilla, f"{node}_U1", "p", (gphase,))
        if last_marker[key] == node:
            append(ancilla, f"{node}_close", "h")


def _exact_runner(runner):
    """The SDK an exact (statevector) run uses for ``runner``: the SDK itself when
    it is one, Qiskit for a registered backend or backend object."""
    return runner if isinstance(runner, str) else "qiskit"


def _sdk_probabilities(software, dag, num_qubits):
    """Basis-state probabilities (index bit ``q`` = qubit ``q``) of a DAG from
    ``software``'s own statevector, or ``None`` if it has none."""
    from Qdislib.utils.sdks import sdk

    entry = sdk(software)
    if entry is None or not entry.statevector or not entry.installed():
        return None
    built = entry.call("from_dag", dag, num_qubits)
    circuit = built[0] if isinstance(built, (list, tuple)) else built
    state = np.asarray(entry.call("statevector", circuit)).reshape(-1)
    n = int(round(math.log2(len(state))))
    probs = np.abs(state) ** 2  # qubit 0 the most significant bit
    return probs.reshape([2] * n).transpose(list(range(n))[::-1]).reshape(-1)


def _exact_fragment(dag, software="qibo"):
    """``(corr, original)`` of a fragment from its statevector (qibo or Qiskit).

    A mid-circuit measurement becomes a CNOT onto a fresh qubit (deferred
    measurement), whose Z then stands for the cut-point outcome: the gate cut's
    ``-(-1)**m`` factor is ``-Z`` of that qubit, folded in here.
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _fragment_correlators
    from Qdislib.utils.graph_qibo import _update_qubits

    dag = dag.copy()
    used = {q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()}
    if not used:
        return np.ones(1), []
    original = sorted(used)
    fresh = max(used) + 1
    deferred = 0
    for node in list(dag.nodes):
        if (dag.nodes[node].get("gate") or "").lower() == "measure":
            (qubit,) = dag.nodes[node]["qubits"]
            dag.nodes[node].update(gate="cx", qubits=(qubit, fresh + deferred))
            deferred += 1
    renumbered, n = _update_qubits(dag)
    probs = None
    if software in ("qilisdk",):  # SDKs whose own statevector serves exact mode
        probs = _sdk_probabilities(software, renumbered, n)
    if probs is None:  # Qdislib's own simulator: no SDK circuit to build
        from Qdislib.utils.statevector import dag_probabilities

        probs = dag_probabilities(renumbered, n)
    if probs is not None:
        from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard

        corr = _fast_walsh_hadamard(probs)
    elif software != "qibo":
        from qiskit.quantum_info import Statevector

        from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard
        from Qdislib.utils.graph_qiskit import dag_to_circuit_qiskit

        built = dag_to_circuit_qiskit(renumbered, n)
        qc = built[0] if isinstance(built, (list, tuple)) else built
        probs = Statevector(qc.remove_final_measurements(inplace=False)).probabilities()
        corr = _fast_walsh_hadamard(probs)  # qiskit index bit q = qubit q
    else:
        corr = _fragment_correlators(renumbered, n, n, "qibo", None, None, False, None)
    kept = len(original)
    folded = corr[np.arange(2**kept) | (((1 << deferred) - 1) << kept)]
    return (-1) ** deferred * folded, original


class _FragmentCircuit(Subcircuit):
    """One fragment circuit of a :func:`circuit_cutting` cut (three-step workflow).

    Running it -- with :func:`~Qdislib.expectation_value` or
    :func:`~Qdislib.sample_subcircuit`, either works -- returns the fragment's
    ``(corr, original)`` Z correlators, which serve both the expectation value and
    the distribution. ``shots=None`` is exact (statevector, on qibo or Qiskit).

    It holds its fragment's :class:`_FragmentTemplate` and its legs, not a built
    circuit: the circuit is built where it runs (a task), or on first access to
    ``dag`` / ``circuit``. ``dag`` keeps the circuit's qubit numbering; ``circuit``
    is renumbered.
    """

    #: How many drawn configurations use this circuit (``num_samples``); 1 when
    #: every configuration is run.
    draws = 1

    def __init__(self, template, spec, software, qpu=False, qpu_dict=None, mask=None):
        self.template = template
        self.spec = spec  # (gate_config, prep_states, setting)
        #: The observable's global Z mask (expectation lists): running it returns
        #: only the correlators the reconstruction reads. ``None``: all of them.
        self.mask = mask
        self.num_qubits = template.num_qubits
        self.software = software
        self.qpu = qpu
        self.qpu_dict = qpu_dict
        self._dag = None

    @property
    def dag(self):
        if self._dag is None:
            self._dag = self.template.build(*self.spec)[0]
        return self._dag

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_dag"] = None  # rebuilt where needed, not shipped
        return state

    def to_circuit(self, software=None, defer_measurements=False):
        """``[circuit, observable]`` in ``software`` (default: the list's), on qubits
        ``0..n-1``, like every subcircuit's ``to_circuit`` (the 1.x pairs; for just
        the native circuits use the list's :meth:`CutCircuits.circuits`). ``defer_measurements`` moves a gate cut's mid-circuit
        measurement onto an extra qubit read at the end (a CNOT onto it: the same
        channel), for hardware without mid-circuit measurement; it is done anyway
        when ``software`` cannot express the measurement (QiliSDK)."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _defer_all_cut_measurements,
        )
        from Qdislib.utils.backends import _materialize

        if software is not None and software != self.software:
            from Qdislib.core.cutting_algorithms._api import resolve_software

            software = resolve_software(software, None)
        software = software or self.software
        dag, width = self._renumbered()

        def deferred():  # every cut measurement, a qubit cut twice included
            moved, cut_qubits = _defer_all_cut_measurements(dag, width)
            return moved, width + len(cut_qubits)

        if defer_measurements:
            return _materialize(*deferred(), software)
        try:
            return _materialize(dag, width, software)
        except NotImplementedError:  # no mid-circuit measurement in this SDK
            moved, wider = deferred()
            if moved is dag:
                raise
            return _materialize(moved, wider, software)

    @property
    def observable(self):
        return None

    def __repr__(self):
        from Qdislib.core.cutting_algorithms.describe import item_repr

        return item_repr(self)

    def _sample(self, shots, backend, seed, keep_counts=False):
        """``keep_counts`` (with shots): the result keeps the counts, so any
        observable in the list's basis reads back from it, not only the list's."""
        runner = self.software if backend is None else backend
        if shots is None and not is_simulator(runner):
            raise ValueError("shots=None (exact) needs a simulator backend")
        template, values = _ship(self.template)
        if keep_counts and shots is not None:
            return _run_fragment(
                template,
                self.spec,
                runner,
                shots,
                self.qpu,
                self.qpu_dict,
                seed,
                None,
                None,
                True,
                values=values,
            )
        return _run_fragment(
            template,
            self.spec,
            runner,
            shots,
            self.qpu,
            self.qpu_dict,
            seed,
            self.mask,
            values=values,
        )

    def _sample_counts(self, shots, backend, seed):
        """:func:`~Qdislib.sample_subcircuit`'s run: the counts kept."""
        return self._sample(shots, backend, seed, keep_counts=True)

    def _expectation(self, shots, backend, seed, cache, cache_mode):
        if cache is None or cache is False:
            return self._sample(shots, backend, seed)
        # The engine's cache (cache_mode does not apply: the whole result is kept).
        from Qdislib.core.cutting_algorithms.engine_cache import EngineCache, _store

        runner = self.software if backend is None else backend
        if shots is None and not is_simulator(runner):
            raise ValueError("shots=None (exact) needs a simulator backend")
        return EngineCache(_store(cache)).run(self, shots, seed, runner=runner)

    def sample(self, shots=1024, seed=None):
        """Run it (:func:`~Qdislib.sample_subcircuit`). With shots the result keeps
        the counts, so ``reconstruct(results, observable)`` reads any observable
        measured in the list's basis from the same runs; an exact run
        (``shots=None``) keeps only what the list's own observable needs."""
        return self._sample(shots, None, seed, keep_counts=True)


@task(returns=1)
def _run_fragment(
    template,
    spec,
    runner,
    shots,
    qpu,
    qpu_dict,
    seed,
    mask=None,
    cache=None,
    keep_counts=False,
    merge=None,
    values=None,
):
    """Build one fragment circuit from its template and run it: a
    :class:`~Qdislib.core.cutting_algorithms.describe.FragmentResult`
    ``(corr, original)`` (see :func:`_fragment_run`). ``merge``: an earlier run's
    counts of the same circuit, returned added to this run's (``total_shots``).
    ``values``: the free parameters' values (:func:`_ship`), bound here."""
    return _run_fragment_here(
        _unship(template, values),
        spec,
        runner,
        shots,
        qpu,
        qpu_dict,
        seed,
        mask,
        cache,
        keep_counts,
        merge,
    )


@task(returns=1)
def _run_fragments(template, values, units):
    """``circuits_per_task=``: several circuits of one fragment in one task
    (``units``: the arguments of :func:`_run_fragment` after the template, for
    each), their results as one tuple -- one future for the group, so a
    reconstruction task takes a few group futures instead of one per circuit
    (PyCOMPSs handles each element of a collection separately). The template is
    its own argument so PyCOMPSs ships it once (see :func:`_ship`)."""
    template = _unship(template, values)
    return tuple(_run_fragment_here(template, *unit) for unit in units)


def _ship(template):
    """``(template, values)`` to send to a task: a template with values bound is
    sent as its unbound template -- the same object for every set of values, so
    PyCOMPSs serialises it once -- and the values as a float array (in the
    template's parameter order); otherwise the template as it is and ``None``."""
    unbound = getattr(template, "_unbound", None)
    if unbound is None:
        return template, None
    array = template.__dict__.get("_values_array")
    if array is None:
        array = np.array([template.values[p] for p in unbound.parameters])
        template._values_array = array
    return unbound, array


#: The last few templates :func:`_unship` bound, by the objects they came from
#: (kept alive here, so their ids cannot be reused): without PyCOMPSs every
#: circuit of a set of values gets the same objects, and the values are bound
#: once, not once per circuit. Under PyCOMPSs each task has its own copies.
_UNSHIPPED = {}


def _unship(template, values):
    """The template a task runs: :func:`_ship`'s values bound back."""
    if values is None:
        return template
    key = (id(template), id(values))
    hit = _UNSHIPPED.get(key)
    if hit is not None and hit[0] is template and hit[1] is values:
        return hit[2]
    bound = template.bound(dict(zip(template.parameters, map(float, values))))
    if len(_UNSHIPPED) >= 16:
        _UNSHIPPED.pop(next(iter(_UNSHIPPED)))
    _UNSHIPPED[key] = (template, values, bound)
    return bound


def _bind_circuits(circuits, values):
    """The fragment circuits with ``values`` (``{name: float}``) bound: one bound
    template per fragment, shared by its circuits as the unbound one is."""
    import copy

    templates = {}
    bound = []
    for circuit in circuits:
        key = id(circuit.template)
        if key not in templates:
            templates[key] = circuit.template.bound(values)
        new = copy.copy(circuit)
        new.template = templates[key]
        new._dag = None
        bound.append(new)
    return bound


def _parameter_sets(parameter_values, plans):
    """``parameter_values`` as a list of ``{name: float}``, checked against the
    parameters of the planned circuit."""
    from Qdislib.core.cutting_algorithms.parameters import normalise_values

    parameters = plans[0].parameters if plans else ()
    if hasattr(parameter_values, "items"):
        raise ValueError(
            "parameter_values takes a list of value sets (one result each); for one "
            "set pass [values]"
        )
    return [normalise_values(values, parameters) for values in parameter_values]


class _Batches:
    """What :func:`_run` returns with ``circuits_per_task=``: one future per group of
    circuits (each resolving to a tuple of their results), the fragment of each
    group (a group never spans two fragments) and the positions of its circuits in
    the circuit list (a fragment's groups come in the order of its circuits)."""

    def __init__(self, futures, fragments, positions):
        self.futures = futures
        self.fragments = fragments
        self.positions = positions


def _flatten(groups, positions=None):
    """The per-circuit results of the groups of a :class:`_Batches`: in the
    circuits' order with their ``positions``, else the groups one after another."""
    if positions is None:
        return [result for group in groups for result in group]
    results = [None] * sum(len(p) for p in positions)
    for group, where in zip(groups, positions):
        for result, n in zip(group, where):
            results[n] = result
    return results


def _run_fragment_here(
    template,
    spec,
    runner,
    shots,
    qpu,
    qpu_dict,
    seed,
    mask=None,
    cache=None,
    keep_counts=False,
    merge=None,
):
    """:func:`_run_fragment`'s work, in the calling process."""
    from Qdislib.core.cutting_algorithms.describe import FragmentResult

    corr, original = _fragment_run(
        template, spec, runner, shots, qpu, qpu_dict, seed, mask, cache, keep_counts
    )
    if merge is not None:
        from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

        joint = dict(merge[0].joint)
        for outcome, count in corr.joint.items():
            joint[outcome] = joint.get(outcome, 0) + count
        corr = _JointCounts(joint, corr.prep_sorted, corr.k)
    return FragmentResult(corr, original)


def _fragment_run(
    template,
    spec,
    runner,
    shots,
    qpu,
    qpu_dict,
    seed,
    mask=None,
    cache=None,
    keep_counts=False,
):
    """Build one fragment circuit from its template and run it: ``(corr, original)``.

    The circuit is built here, in the task, so the master only sends the small
    ``spec``; the template is the same object for every circuit of the fragment.
    With ``mask`` (an expectation value) ``corr`` keeps only what is read. With a
    ``cache`` request (:mod:`~Qdislib.core.cutting_algorithms.engine_cache`) the
    whole result is kept -- the counts, or every exact correlator -- and stored.
    """
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _component_original_qubits,
        _prepare_sampled_component,
        _sample_prepared_component,
    )

    dag, _ = template.build(*spec)
    if template.values is None and getattr(template, "parameters", ()):
        from Qdislib.core.cutting_algorithms.parameters import require_bound

        require_bound(dag)
    if cache is not None:
        from Qdislib.core.cutting_algorithms.engine_cache import from_payload

        key = cache.lookup_key(dag, runner)
        if key is not None:
            payload = cache.get(key)
            if payload is not None:  # a semantic hit: the order of its qubits
                corr, _ = from_payload(payload)
                return corr, _component_original_qubits(dag)
        if shots is None:
            result = _exact_fragment(dag, _exact_runner(runner))
        else:
            result = _sample_prepared_component(
                _prepare_sampled_component(dag),
                runner,
                shots,
                qpu,
                qpu_dict,
                seed,
                keep_counts=True,
            )
        if key is not None:
            cache.put(key, result)
        return result
    needed = None
    if keep_counts and shots is not None:  # the counts, for the error bars
        return _sample_prepared_component(
            _prepare_sampled_component(dag),
            runner,
            shots,
            qpu,
            qpu_dict,
            seed,
            keep_counts=True,
        )
    if mask is not None:
        original = _component_original_qubits(dag)
        needed = _needed_masks(template, original, spec[2], mask)
    if shots is None:
        corr, original = _exact_fragment(dag, _exact_runner(runner))
        if needed is None:
            return corr, original
        return {i: float(corr[i]) for i in needed}, original
    return _sample_prepared_component(
        _prepare_sampled_component(dag), runner, shots, qpu, qpu_dict, seed, needed
    )


class _FragmentReads:
    """Fragment ``f``'s circuit results, read for a Pauli on its cut wires.

    The fragment's own output qubits are ``out``; a read returns one value per
    local mask over them (``local_masks``: all ``2**len(out)``, or the one mask of
    an expectation value). A value is the circuit's correlator at the mask whose
    output bits come from the local mask, whose cut-wire bits are the Pauli's
    readout bits, and whose ancilla bits are always set.
    """

    def __init__(self, plan, f, keys, readouts, results, mask=None):
        fragment = plan.fragments[f]
        entries = []
        outputs = set()
        from Qdislib.core.cutting_algorithms.gate_cutting import _JointCounts

        for corr, original in results:
            if isinstance(corr, _JointCounts):  # a cached run: its counts
                corr = corr.correlators()
            if isinstance(corr, dict) and len(corr) == 2 ** len(original):
                corr = np.array([corr[m] for m in range(len(corr))])
            roles = []
            for qubit in original:
                if qubit >= plan.ancilla_base:
                    roles.append(("ancilla", None))
                elif plan.cut_source_index.get(qubit) not in fragment["measures"]:
                    output = plan.prepared_to_orig.get(qubit, qubit)
                    roles.append(("output", output))
                    outputs.add(output)
                else:
                    roles.append(("cut", None))
            positions = [
                original.index(fragment["cut_qubits"][c]) for c in fragment["measures"]
            ]
            if not isinstance(corr, dict):
                corr = np.asarray(corr)
            entries.append((corr, roles, positions))
        self.out = sorted(outputs)
        where = {q: j for j, q in enumerate(self.out)}
        #: True when ``mask`` is the masks of a group of terms: the local masks
        #: are one per term, an index the fragments share (see :func:`_contract`).
        self.shared = isinstance(mask, (tuple, list))
        if mask is None:
            local = np.arange(2 ** len(self.out), dtype=np.int64)
        else:
            local = np.array(
                [
                    sum(((m >> q) & 1) << j for j, q in enumerate(self.out))
                    for m in (mask if self.shared else [mask])
                ],
                dtype=np.int64,
            )
        self.local_masks = local
        self.prepared = []
        for corr, roles, positions in entries:
            index = np.zeros(len(local), dtype=np.int64)
            for i, (role, output) in enumerate(roles):
                if role == "ancilla":
                    index |= 1 << i
                elif role == "output":
                    index |= ((local >> where[output]) & 1) << i
            self.prepared.append((corr, index, positions))
        self.by_key = {}
        for n, key in enumerate(keys):
            self.by_key.setdefault(key, []).append(n)
        self.readouts = readouts
        self.cache = {}

    def read(self, key, letters):
        """The circuits under ``key`` read for the Pauli ``letters`` on the cut
        wires, averaged when several circuits read it."""
        if (key, letters) not in self.cache:
            values = []
            for n in self.by_key.get(key, ()):
                corr, index, positions = self.prepared[n]
                for pauli, sign, bits in self.readouts[n]:
                    if tuple(pauli) == letters:
                        at = index | sum(1 << positions[b] for b in bits)
                        if hasattr(corr, "values_at"):  # a run's counts
                            values.append(sign * corr.values_at(at))
                        elif isinstance(corr, dict):  # only the read correlators
                            values.append(sign * np.array([corr[int(i)] for i in at]))
                        else:
                            values.append(sign * corr[at])
            if not values:
                raise RuntimeError(f"no circuit reads {letters} for {key}")
            self.cache[(key, letters)] = np.mean(values, axis=0)
        return self.cache[(key, letters)]


def _fragment_value(plan, f, reads, digits, config, wire):
    """Fragment ``f``'s value for one configuration (its own part of it): one
    entry per local mask."""
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _SETTING_PAULI,
        _WIRE_PREP,
    )

    fragment = plan.fragments[f]
    leg = plan.gate_legs(f, digits, config)
    letters = tuple(_SETTING_PAULI[wire[i]] for i in fragment["measures"])
    preps = [wire[i] for i in fragment["preps"]]
    if not plan.optimized:
        states = tuple(_WIRE_DIGIT_STATE[d] for d in preps)
        return reads.read((f, leg, states), letters)
    value = 0.0
    for choice in itertools.product(*(_WIRE_PREP[d] for d in preps)):
        coefficient = math.prod(weight for weight, _ in choice)
        states = tuple(state for _, state in choice)
        value = value + coefficient * reads.read((f, leg, states), letters)
    return value


def _fragment_tensor(plan, f, reads):
    """Fragment ``f`` as a tensor: one axis per bond it touches (over the bond's
    terms), then one over its local masks. ``(tensor, labels, out)``."""
    fragment = plan.fragments[f]
    bonds = plan.fragment_bonds(f)
    dims = [len(plan.bond_weights(bond)) for bond in bonds]
    tensor = np.zeros(dims + [len(reads.local_masks)])
    settings = plan._hadamard_settings
    num_mp = len(plan.gate_cuts)
    for index in itertools.product(*(range(d) for d in dims)):
        at = dict(zip(bonds, index))
        digits = [0] * num_mp
        for k, _ in fragment["gate_sides"]:
            digits[k - 1] = at[("gate", k)]
        config = [None] * len(plan.groups)
        for g in fragment["had_groups"]:
            config[g] = settings[g][0][at[("hadamard", g)]]
        wire = {
            i: plan.wire_digits[i - 1][at[("wire", i)]]
            for i in fragment["measures"] + fragment["preps"]
        }
        tensor[index] = _fragment_value(plan, f, reads, digits, config, wire)
    out = ("term",) if reads.shared else ("out", f)
    return tensor, bonds + [out], reads.out


def _contract(plan, parts):
    """Contract the fragment tensors over their bonds, each bond weighted by its
    terms' weights: the tensor over every fragment's local masks, fragment order.

    Pairwise, smallest result first (a chain contracts as a matrix product).
    """
    items, placed = [], set()
    for tensor, labels, _ in parts:
        for axis, label in enumerate(labels[:-1]):
            if label not in placed:  # each bond's weights go on one tensor
                shape = [1] * tensor.ndim
                shape[axis] = -1
                tensor = tensor * plan.bond_weights(label).reshape(shape)
                placed.add(label)
        items.append((tensor, list(labels)))
    # one open index per fragment, or the terms' index they all share
    open_labels = list(dict.fromkeys(labels[-1] for _, labels, _ in parts))

    def einsum(operands, out):
        ids = {}
        args = []
        for tensor, labels in operands:
            args += [tensor, [ids.setdefault(label, len(ids)) for label in labels]]
        return np.einsum(*args, [ids[label] for label in out])

    def kept(labels, others):
        rest = {label for _, other in others for label in other}
        return [l for l in dict.fromkeys(labels) if l in open_labels or l in rest]

    while len(items) > 1:
        best = None
        for i, j in itertools.combinations(range(len(items)), 2):
            (a, la), (b, lb) = items[i], items[j]
            others = [item for n, item in enumerate(items) if n not in (i, j)]
            out = kept(la + lb, others)
            size = dict(zip(la, a.shape)) | dict(zip(lb, b.shape))
            cost = math.prod(size[l] for l in out)
            if best is None or cost < best[0]:
                best = (cost, i, j, out)
        _, i, j, out = best
        merged = (einsum([items[i], items[j]], out), out)
        items = [item for n, item in enumerate(items) if n not in (i, j)] + [merged]
    return einsum(items, open_labels)


def _draws_sum(plan, readers, draws):
    """The sum over drawn configurations (``num_samples``), each weighted
    ``gamma * sign / num_samples`` per draw: the tensor over every fragment's local
    masks, as :func:`_contract` returns."""
    total = 0.0
    for c, (digits, config, wire, weight, count) in enumerate(draws):
        wire = {i + 1: d for i, d in enumerate(wire)}
        contribution = np.array(weight * count)
        for f, reads in enumerate(readers):
            value = _fragment_value(plan, f, reads, digits, config, wire)
            if reads.shared:  # one entry per term, the same in every fragment
                contribution = contribution * value
            else:
                contribution = np.multiply.outer(contribution, value)
        total = total + contribution
    return total


def _split_by_fragment(plan, keys, readouts, results):
    """``keys``, ``readouts`` and ``results`` of each fragment (a key starts with
    its fragment)."""
    parts = [([], [], []) for _ in plan.fragments]
    for key, readout, result in zip(keys, readouts, results):
        part = parts[key[0]]
        part[0].append(key)
        part[1].append(readout)
        part[2].append(result)
    return parts


def _cut_tensor(plan, keys, readouts, results, mask=None, draws=None):
    """The reconstructed Z correlators as a tensor over each fragment's local masks
    (``mask``: only that global mask), and each fragment's output qubits."""
    readers = [
        _FragmentReads(plan, f, k, r, res, mask)
        for f, (k, r, res) in enumerate(
            _split_by_fragment(plan, keys, readouts, results)
        )
    ]
    if draws is not None:
        return _draws_sum(plan, readers, draws), [reads.out for reads in readers]
    parts = [_fragment_tensor(plan, f, reads) for f, reads in enumerate(readers)]
    return _contract(plan, parts), [out for _, _, out in parts]


def _global_correlators(tensor, outs, width):
    """``<Z^A>`` for every global mask ``A`` from the per-fragment tensor (qubits in
    no fragment's output are idle, so their bits do not change the value)."""
    masks = np.arange(2**width, dtype=np.int64)
    index = tuple(
        (
            sum(((masks >> q) & 1) << j for j, q in enumerate(out))
            if out
            else np.zeros_like(masks)
        )
        for out in outs
    )
    return np.asarray(tensor)[index]


def _observable_mask(observable, width):
    """The global Z mask of a Pauli string after its basis change."""
    if observable is None:
        return (1 << width) - 1
    return sum(1 << q for q, pauli in enumerate(observable) if pauli != "I")


def _distribution(z_correlators, width):
    from Qdislib.core.cutting_algorithms.gate_cutting import _fast_walsh_hadamard

    dim = 2**width
    probabilities = _fast_walsh_hadamard(z_correlators) / dim
    keys = ["".join(str((x >> q) & 1) for q in range(width)) for x in range(dim)]
    norm = float(np.abs(probabilities).sum()) or 1.0
    return {key: abs(float(p)) / norm for key, p in zip(keys, probabilities)}


def _expectation(plan, keys, readouts, results, observable, draws):
    """The expectation value of a term (a Pauli string) or of a commuting group
    (``[(coefficient, mask)]``) from the fragment results."""
    if isinstance(observable, list):
        masks = tuple(mask for _, mask in observable)
        tensor, _ = _cut_tensor(plan, keys, readouts, results, masks, draws)
        values = np.asarray(tensor).reshape(-1)
        return float(sum(c * v for (c, _), v in zip(observable, values)))
    mask = _observable_mask(observable, plan.width)
    tensor, _ = _cut_tensor(plan, keys, readouts, results, mask, draws)
    return float(np.asarray(tensor).reshape(-1)[0])


@task(returns=1, results=COLLECTION_IN)
def _reconstruct_cut_task(
    plan,
    keys,
    readouts,
    results,
    mode,
    observable,
    draws=None,
    bootstrap=0,
    seed=None,
    positions=None,
):
    """Expectation (``mode="expectation"``) or distribution from the fragment
    results, as one task fed straight by the fragment futures (or by the group
    futures of a :class:`_Batches`, with their ``positions``). With ``bootstrap``
    resamples, an expectation value is an
    :class:`~Qdislib.core.cutting_algorithms.uncertainty.Estimate` with its
    standard error."""
    if positions is not None:
        results = _flatten(results, positions)
    if mode == "expectation":
        value = _expectation(plan, keys, readouts, results, observable, draws)
        if not bootstrap:
            return value
        from Qdislib.core.cutting_algorithms.uncertainty import Estimate
        from Qdislib.core.cutting_algorithms.uncertainty import bootstrap as resample

        error, shots = resample(
            lambda res, drs: _expectation(plan, keys, readouts, res, observable, drs),
            results,
            draws,
            bootstrap,
            seed,
        )
        return Estimate(value, error, shots)

    def distribution(res, drs):
        tensor, outs = _cut_tensor(plan, keys, readouts, res, None, drs)
        return _distribution(_global_correlators(tensor, outs, plan.width), plan.width)

    value = distribution(results, draws)
    if not bootstrap:
        return value
    from Qdislib.core.cutting_algorithms.uncertainty import Distribution
    from Qdislib.core.cutting_algorithms.uncertainty import bootstrap as resample

    errors, shots = resample(distribution, results, draws, bootstrap, seed)
    return Distribution(value, errors, shots)


@task(returns=1, results=COLLECTION_IN)
def _fragment_tensor_task(plan, f, keys, readouts, results, batched=False):
    """One fragment's tensor over its bonds and all its local masks, reduced where
    its results are (a distribution's fragment results are ``2**width`` each;
    ``batched``: the fragment's group futures of a :class:`_Batches`)."""
    if batched:
        results = _flatten(results)
    return _fragment_tensor(plan, f, _FragmentReads(plan, f, keys, readouts, results))


@task(returns=1, parts=COLLECTION_IN)
def _contract_distribution_task(plan, parts):
    tensor = _contract(plan, parts)
    outs = [out for _, _, out in parts]
    return _distribution(_global_correlators(tensor, outs, plan.width), plan.width)


def _reconstruct(
    plan,
    keys,
    readouts,
    results,
    mode,
    observable,
    draws=None,
    bootstrap=0,
    seed=None,
):
    """The reconstruction as a future. An expectation value is one task (the
    contraction is cheap); a distribution reduces each fragment's results in its own
    task first, then contracts the small tensors."""
    batched = isinstance(results, _Batches)
    if mode == "distribution" and draws is None and not bootstrap:
        if batched:
            split = _split_by_fragment(plan, keys, readouts, [None] * len(keys))
            parts = [
                _fragment_tensor_task(
                    plan,
                    f,
                    k,
                    r,
                    [g for g, gf in zip(results.futures, results.fragments) if gf == f],
                    True,
                )
                for f, (k, r, _) in enumerate(split)
            ]
        else:
            parts = [
                _fragment_tensor_task(plan, f, k, r, res)
                for f, (k, r, res) in enumerate(
                    _split_by_fragment(plan, keys, readouts, results)
                )
            ]
        return _contract_distribution_task(plan, parts)
    if batched:
        return _reconstruct_cut_task(
            plan,
            keys,
            readouts,
            results.futures,
            mode,
            observable,
            draws,
            bootstrap,
            seed,
            results.positions,
        )
    return _reconstruct_cut_task(
        plan, keys, readouts, results, mode, observable, draws, bootstrap, seed
    )


@task(returns=1, values=COLLECTION_IN)
def _weighted_sum(coefficients, values):
    """``sum(c * v)`` of a weighted Pauli sum's terms, as a task (so ``sync=False``
    returns one future); an ``Estimate`` when the values carry errors."""
    from Qdislib.core.cutting_algorithms.uncertainty import combine

    return combine(coefficients, values)


def _prepare_cut(
    circuit,
    cuts,
    observable,
    hadamard,
    optimized,
    joint,
    software,
    qpu,
    qpu_dict,
    num_samples=None,
    seed=None,
    masks=None,
):
    """The plan, the drawn configurations (``None`` when all are run), and one
    :class:`_FragmentCircuit` per fragment circuit to run; ``software`` is what
    runs them (an SDK, registered backend or backend object). No circuit is built
    here: each is built where it runs."""
    plan = _CutPlan(circuit, cuts, observable, hadamard, optimized, joint)
    _warn_cheaper_hadamard(plan)
    draws = None
    if num_samples is not None:
        if int(num_samples) < 1:
            raise ValueError(f"num_samples must be at least 1, got {num_samples}")
        draws = plan.draw(num_samples, next_seed(seed, "draws"))
    templates = [_FragmentTemplate(plan, f) for f in range(len(plan.fragments))]
    # An expectation value reads one mask: its circuits return only that part.
    if masks is not None:  # the terms of a commuting group, read from one run
        mask = tuple(masks)
    else:
        mask = None if observable is None else _observable_mask(observable, plan.width)
    keys, readouts, circuits = [], [], []
    for key, f, gate_config, states, setting, count in plan.units(draws):
        keys.append(key)
        readouts.append(_setting_readout(len(plan.fragments[f]["measures"]), setting))
        circuit_ = _FragmentCircuit(
            templates[f], (gate_config, states, setting), software, qpu, qpu_dict, mask
        )
        circuit_.draws = count
        circuits.append(circuit_)
    return plan, draws, keys, readouts, circuits


def _warn_cheaper_hadamard(plan):
    """Warn when a measure-and-prepare cut of a gate costs clearly more than the
    Hadamard cut would. For a single rotation (``rzz``, ``cp``, ``crz``, ...) they
    are the same; a gate with two or three interactions (``iswap``, ``swap``) is
    cut as a product of rotations, which the Hadamard cut beats."""
    cheaper = [
        (name, mp, gamma)
        for name, (mp, gamma) in plan.general_gate_cuts.items()
        if gamma is not None and gamma < 0.9 * mp
    ]
    if not cheaper:
        return
    import warnings

    detail = "; ".join(
        f"{name}: gamma {mp:g} measure-and-prepare, {had:.3g} Hadamard"
        for name, mp, had in cheaper
    )
    names = [name for name, _, _ in cheaper]
    warnings.warn(
        f"The Hadamard cut is cheaper for {names} ({detail}; the shots grow as "
        f"gamma**2). Cut them with {{'hadamard': {names}}} or hadamard=True.",
        UserWarning,
        stacklevel=_caller_stacklevel(),
    )


def _caller_stacklevel():
    """The ``stacklevel`` that points a warning at the first frame outside Qdislib
    (the user's call), however deep inside the library it is raised."""
    import os
    import sys

    package = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    frame, level = sys._getframe(1), 1
    while frame is not None and frame.f_code.co_filename.startswith(package):
        frame, level = frame.f_back, level + 1
    return level


def _run(circuits, shots, seed, cache=None, keep_counts=False, batch=None):
    """Run every circuit, one task each (``batch``: that many per task); with drawn
    configurations a circuit gets ``shots`` per draw that uses it (the shots follow
    the weights). With a ``cache``
    (:func:`~Qdislib.core.cutting_algorithms.engine_cache.resolve`), a circuit
    already run is read back instead. ``keep_counts`` returns each run's counts
    (for the error bars; a cached run keeps them anyway)."""
    if batch is not None:
        if cache is not None:
            raise ValueError("circuits_per_task= does not combine with cache= yet")
        units = [
            (
                circuit.spec,
                circuit.software,
                None if shots is None else shots * circuit.draws,
                circuit.qpu,
                circuit.qpu_dict,
                next_seed(seed, n),
                None if keep_counts else circuit.mask,
                None,
                keep_counts,
            )
            for n, circuit in enumerate(circuits)
        ]
        # a group stays inside one fragment (the distribution's per-fragment
        # reduction takes only its own fragment's groups)
        by_fragment = {}
        for n, circuit in enumerate(circuits):
            by_fragment.setdefault(id(circuit.template), []).append(n)
        futures, fragments, positions = [], [], []
        for members in by_fragment.values():
            f = circuits[members[0]].template.f
            for start in range(0, len(members), batch):
                where = members[start : start + batch]
                futures.append(
                    _run_fragments(
                        *_ship(circuits[where[0]].template), [units[n] for n in where]
                    )
                )
                fragments.append(f)
                positions.append(where)
        return _Batches(futures, fragments, positions)
    results = []
    for n, circuit in enumerate(circuits):
        unit_shots = None if shots is None else shots * circuit.draws
        if cache is None and keep_counts:
            results.append(
                _run_fragment(
                    _ship(circuit.template)[0],
                    circuit.spec,
                    circuit.software,
                    unit_shots,
                    circuit.qpu,
                    circuit.qpu_dict,
                    next_seed(seed, n),
                    None,
                    None,
                    True,
                    values=_ship(circuit.template)[1],
                )
            )
        elif cache is None:
            results.append(circuit.sample(unit_shots, next_seed(seed, n)))
        else:
            results.append(cache.run(circuit, unit_shots, next_seed(seed, n)))
    return results


def _check_batch(batch):
    if batch is None:
        return
    if isinstance(batch, bool) or not isinstance(batch, int) or batch < 1:
        raise ValueError(f"circuits_per_task= takes a number >= 1, got {batch!r}")


def _check_exact(shots, backend):
    if shots is None and backend is not None and not is_simulator(backend):
        raise ValueError(
            "shots=None (exact) is a statevector simulation: use a simulator backend "
            f"(one of {list(SDKS)}), not {backend!r}."
        )


@standard_arguments
def circuit_cutting(
    circuit,
    cuts,
    observable=None,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    hadamard=False,
    optimized=False,
    joint=True,
    num_samples=None,
    cache=None,
    cache_key=None,
    cache_refresh=False,
    error=False,
    total_shots=None,
    circuits_per_task=None,
    parameter_values=None,
):
    """Reconstruct ``<observable>`` from a circuit cut by any combination of cuts.

    Measure-and-prepare gate cuts, Hadamard gate cuts and wire cuts can be mixed
    freely, in any number; their quasiprobability decompositions tensor. With
    ``optimized=True`` each fragment runs once per distinct setting of the cuts that
    touch it, so the number of circuits grows as a *sum* over fragments, not as a
    product over cuts (see :func:`cutting_overhead`).

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: ``{"gate": [...], "hadamard": [...], "wire": [...]}`` (any subset;
        ``find_cut(..., mixed=True)`` returns one), or a list mixing gate names and
        ``(source, target)`` wire edges.
    :param observable: A Pauli string (``X``/``Y`` are rotated into ``Z``), or a
        weighted Pauli sum ``[(coefficient, string), ...]``, whose qubit-wise
        commuting terms share one run; ``None`` is all ``Z``.
    :param backend: Where the fragments run (see :func:`gate_cutting`).
    :param shots: Shots per fragment circuit (per draw with ``num_samples``);
        ``None`` is exact (statevector).
    :param seed: Seed for the simulators' sampling (each circuit its own stream)
        and for the ``num_samples`` draws.
    :param sync: ``False`` returns the value as a future (PyCOMPSs).
    :param hadamard: Cut the ``"gate"`` cuts (or the gate names of a list) with the
        Hadamard test (:func:`gate_cutting_hadamard`) instead of measure-and-prepare.
    :param optimized: Run each fragment once per distinct setting of its own cuts.
    :param joint: Measure the wires a fragment sends together (with ``optimized``)
        and cut the Hadamard gates crossing one boundary together.
    :param num_samples: Use this many configurations, in proportion to
        ``|weight| / gamma``, instead of running them all (quasiprobability
        Monte-Carlo): only the circuits they use run, each with ``shots`` per draw
        that uses it. The draws are allocated, not drawn at random: every
        configuration gets its share when they all fit in ``num_samples``
        (equivalent to running them all), otherwise the cuts that fit are
        allocated and the rest drawn (see :meth:`_CutPlan.draw`). Unbiased. The
        error has two parts: which configurations are drawn -- at most
        ``gamma / sqrt(num_samples)``, whatever the shots -- and the shot noise. When
        the configurations do not all fit, the first usually dominates and more
        shots per draw do not help: for a given budget use more draws with fewer
        shots each. The way to run cuts with too many configurations to list
        (``6**n`` for gate cuts, ``8**n`` for wire cuts; with ``optimized=True`` the
        circuits per fragment are ``5**n`` and ``3**n + 4**n``), and the shots then
        follow the weights (Hadamard cuts).
    :param cache: Reuse fragment circuits already run: ``True`` (in memory; LMDB
        in the working directory under PyCOMPSs), ``"memory"``, ``"lmdb:<path>"``,
        ``"redis[://host:port/db]"`` or a ``CircuitCache``; ``False`` turns off a
        session cache (:func:`~Qdislib.use_cache`). A hit is the circuit's whole
        result (counts or exact correlators), so it serves any observable.
    :param cache_key: ``"exact"`` (default: the same circuit, backend and shots)
        or ``"semantic"`` (also equivalent circuits, by a ZX hash in the task).
    :param cache_refresh: Run every circuit again and replace what is stored
        (fresh samples, e.g. to repeat a run for error bars).
    :param error: ``True`` (or a number of resamples) returns an
        :class:`~Qdislib.Estimate`: the value -- still a ``float`` -- with
        ``.error``, its standard error by bootstrap over the shots (and the
        ``num_samples`` draws), and ``.shots``. Costs reconstruction time only.
        With ``num_samples`` and more configurations than draws it errs high
        (resampling the draws and their shots separately).
    :param total_shots: One shot budget for the whole call instead of ``shots``
        per circuit: a pilot round (20%) split by each circuit's weight in the
        reconstruction, then the rest by weight x the spread the pilot measured
        (the master waits for the pilot). With ``num_samples``, ``shots`` per draw
        and fragment is ``total_shots // (num_samples * fragments)``, so at most
        the budget is spent.
    :param circuits_per_task: Run this many fragment circuits per task instead of
        one each: fewer, longer tasks when the circuits are small and many (a
        task's overhead outweighs a fast simulation). The results are the same.
    :param parameter_values: For a circuit with free parameters (Qiskit
        ``Parameter`` objects, anywhere but in a cut gate): a list of value sets --
        each ``{Parameter or name: value}`` or a list in ``circuit.parameters``
        order -- with one result for each, in a list. The cut is planned once
        and every set runs its own circuits (each its own seed stream; with
        ``num_samples`` they share the draws).
    :return: The reconstructed expectation value (a list of them with
        ``parameter_values``).
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _observable_terms

    from Qdislib.core.cutting_algorithms.engine_cache import resolve

    from Qdislib.core.cutting_algorithms.uncertainty import resamples

    software = resolve_backend(backend, circuit)
    engine_cache = resolve(cache, cache_key, cache_refresh)
    count = resamples(error)
    budget = None
    if total_shots is not None:
        if int(total_shots) < 1:
            raise ValueError(f"total_shots must be positive, got {total_shots}")
        if num_samples is not None:
            shots = 1  # set per draw below, once the fragments are known
        else:
            if engine_cache is not None:
                raise ValueError(
                    "total_shots= splits the shots between two rounds; it does not "
                    "combine with cache= yet"
                )
            if circuits_per_task is not None:
                raise ValueError(
                    "total_shots= does not combine with circuits_per_task= yet"
                )
            budget, shots = int(total_shots), 1
    _check_exact(shots, backend)
    _check_batch(circuits_per_task)
    terms = _observable_terms(default_observable(observable, circuit))
    # identity terms are constants; the others run in qubit-wise commuting groups
    constant = sum(c for c, term in terms if set(term) <= {"I"})
    groups = _commuting_groups([(c, t) for c, t in terms if not set(t) <= {"I"}])
    prepared = []
    for g, (basis, members) in enumerate(groups):
        group_seed = None if seed is None else next_seed(seed, g)
        masks = [_term_mask(term) for _, term in members]
        plan, draws, keys, readouts, circuits = _prepare_cut(
            circuit,
            cuts,
            basis,
            hadamard,
            optimized,
            joint,
            software,
            False,
            None,
            num_samples,
            group_seed,
            masks=masks,
        )
        prepared.append((plan, draws, keys, readouts, circuits, members, masks))
    if total_shots is not None and num_samples is not None:
        # a draw runs one circuit per fragment, in every group
        per_draw = int(num_samples) * sum(len(p[0].fragments) for p in prepared)
        shots = max(1, int(total_shots) // max(1, per_draw))

    def evaluate(prepared, seed):
        if budget is not None:
            from Qdislib.core.cutting_algorithms.budget import run_with_budget

            all_results = run_with_budget(
                [
                    (plan, circuits, keys, sum(abs(c) for c, _ in members))
                    for plan, _, keys, _, circuits, members, _ in prepared
                ],
                budget,
                seed,
            )
        else:
            all_results = [
                _run(
                    circuits,
                    shots,
                    None if seed is None else next_seed(seed, g),
                    engine_cache,
                    keep_counts=count > 0,
                    batch=circuits_per_task,
                )
                for g, (_, _, _, _, circuits, _, _) in enumerate(prepared)
            ]
        values = []
        for g, ((plan, draws, keys, readouts, _, members, masks), results) in enumerate(
            zip(prepared, all_results)
        ):
            group = [(c, mask) for (c, _), mask in zip(members, masks)]
            values.append(
                _reconstruct(
                    plan,
                    keys,
                    readouts,
                    results,
                    "expectation",
                    group,
                    draws,
                    count,
                    None if seed is None else next_seed(seed, "bootstrap", g),
                )
            )
        if constant:
            values.append(float(constant))
        if not values:
            total = 0.0
        elif len(values) == 1:
            total = values[0]
        else:
            total = _weighted_sum([1.0] * len(values), values)
        if count and not values:
            from Qdislib.core.cutting_algorithms.uncertainty import Estimate

            total = Estimate(total)
        return total

    if parameter_values is None:
        total = evaluate(prepared, seed)
        return total if not sync else compss_wait_on(total)
    totals = [
        evaluate(
            [
                (
                    plan,
                    draws,
                    keys,
                    readouts,
                    _bind_circuits(circuits, values),
                    members,
                    masks,
                )
                for plan, draws, keys, readouts, circuits, members, masks in prepared
            ],
            None if seed is None else next_seed(seed, "parameter_values", i),
        )
        for i, values in enumerate(
            _parameter_sets(parameter_values, [p[0] for p in prepared])
        )
    ]
    return totals if not sync else compss_wait_on(totals)


def _term_mask(term):
    """The Z mask of a Pauli string (after its basis change)."""
    return sum(1 << q for q, pauli in enumerate(term) if pauli != "I")


def _commuting_groups(terms):
    """``[(basis, [(coefficient, term), ...]), ...]``: the terms of a Pauli sum in
    qubit-wise commuting groups (on each qubit they have the same Pauli or ``I``),
    each measured in one basis, so a group runs one set of fragment circuits and
    every term is read from it. Greedy, the terms with most non-identity qubits
    first."""
    order = sorted(range(len(terms)), key=lambda k: -sum(p != "I" for p in terms[k][1]))
    groups = []  # [basis letters, members]
    for k in order:
        coefficient, term = terms[k]
        for basis, members in groups:
            if all(a == "I" or b == "I" or a == b for a, b in zip(basis, term)):
                basis[:] = [b if a == "I" else a for a, b in zip(basis, term)]
                members.append((coefficient, term))
                break
        else:
            groups.append([list(term), [(coefficient, term)]])
    return [("".join(basis), members) for basis, members in groups]


@standard_arguments
def circuit_cutting_sampling(
    circuit,
    cuts,
    *,
    backend=None,
    shots=1024,
    seed=None,
    sync=True,
    hadamard=False,
    optimized=False,
    joint=True,
    num_samples=None,
    cache=None,
    cache_key=None,
    cache_refresh=False,
    error=False,
    circuits_per_task=None,
    parameter_values=None,
):
    """Reconstruct the output **probability distribution** of a circuit cut by any
    combination of cuts, from the same fragment circuits as :func:`circuit_cutting`.
    Arguments as there.

    :param error: ``True`` (or a number of resamples): return a
        :class:`~Qdislib.Distribution` -- still the ``{bitstring: probability}``
        dict -- with ``.errors``, each probability's standard error by bootstrap
        over the shots, and ``.shots``.
    :param circuits_per_task: Fragment circuits per task, as in :func:`circuit_cutting`.
    :param parameter_values: A list of value sets for the circuit's free
        parameters, one distribution each, as in :func:`circuit_cutting`.
    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings), the signed
        estimate abs-normalised like :func:`gate_cutting_sampling` (a list of them
        with ``parameter_values``).
    """
    _check_exact(shots, backend)
    _check_batch(circuits_per_task)
    plan, draws, keys, readouts, circuits = _prepare_cut(
        circuit,
        cuts,
        None,
        hadamard,
        optimized,
        joint,
        resolve_backend(backend, circuit),
        False,
        None,
        num_samples,
        seed,
    )
    from Qdislib.core.cutting_algorithms.engine_cache import resolve

    from Qdislib.core.cutting_algorithms.uncertainty import resamples

    count = resamples(error)
    engine_cache = resolve(cache, cache_key, cache_refresh)

    def evaluate(circuits, seed):
        results = _run(
            circuits,
            shots,
            seed,
            engine_cache,
            keep_counts=count > 0,
            batch=circuits_per_task,
        )
        return _reconstruct(
            plan,
            keys,
            readouts,
            results,
            "distribution",
            None,
            draws,
            count,
            None if seed is None else next_seed(seed, "bootstrap"),
        )

    if parameter_values is None:
        future = evaluate(circuits, seed)
        return compss_wait_on(future) if sync else future
    futures = [
        evaluate(
            _bind_circuits(circuits, values),
            None if seed is None else next_seed(seed, "parameter_values", i),
        )
        for i, values in enumerate(_parameter_sets(parameter_values, [plan]))
    ]
    return compss_wait_on(futures) if sync else futures


_BIT_ORDERS = ("q0_first", "q0_last")


def result_from_counts(subcircuit, counts, bit_order="q0_first", keep_counts=True):
    """What :func:`~Qdislib.expectation_value` would return for ``subcircuit``,
    from measurement counts of its circuit run with your own tools.

    The circuit is the one ``subs.circuits(software, defer_measurements=True)``
    builds, with **every qubit measured at the end** (a gate cut's mid-circuit
    measurement is then an extra qubit, the last ones). ``counts`` maps each
    outcome to its count (or probability): a bitstring with qubit 0 first
    (``bit_order="q0_first"``, Qdislib's order) or last (``"q0_last"``, Qiskit's
    ``get_counts()``; of a key with several registers, the first is read -- the
    one ``measure_all()`` adds). Outcomes may also be tuples of bits.

    :param keep_counts: Keep the counts (default), so ``subs.reconstruct(results,
        observable)`` can read any observable in the list's basis. ``False`` keeps
        only the correlators the list's own observable reads, as Qdislib's own runs
        do: a few numbers instead of the counts, so the reconstruction is quick
        (for a 30-qubit circuit at 400k shots, seconds instead of minutes).
    :return: A :class:`~Qdislib.core.cutting_algorithms.describe.FragmentResult`
        for ``subs.reconstruct``.
    """
    from Qdislib.core.cutting_algorithms.describe import FragmentResult
    from Qdislib.core.cutting_algorithms.gate_cutting import (
        _JointCounts,
        _cut_measurements,
        _fold_cut_bits,
        _prepare_sampled_component,
    )

    if bit_order not in _BIT_ORDERS:
        raise ValueError(f"bit_order must be one of {_BIT_ORDERS}, got {bit_order!r}")
    dag, width, prep_sorted, original = _prepare_sampled_component(subcircuit.dag)
    if dag is None:  # an empty fragment
        return FragmentResult({0: 1.0}, [])
    cut_qubits = _cut_measurements(dag)  # one per cut measurement: the extra qubits
    expected = width + len(cut_qubits)
    mask = getattr(subcircuit, "mask", None)
    if not keep_counts and mask is not None:  # vectorised: no per-outcome Python
        needed = _needed_masks(subcircuit.template, original, subcircuit.spec[2], mask)
        bits = _outcome_bits(counts, bit_order, expected, width, prep_sorted)
        weights = np.fromiter(counts.values(), dtype=float, count=len(counts))
        signed = weights / (weights.sum() or 1.0)
        if (
            bits.shape[1] > width
        ):  # each cut measurement's sign, as the engine folds them
            signed = signed * np.prod(1.0 - 2.0 * (1 - bits[:, width:]), axis=1)
        finals = bits[:, :width].astype(np.int64)
        corr = {}
        for m in needed:
            select = (int(m) >> np.arange(width)) & 1
            parity = (finals @ select) & 1
            corr[int(m)] = float(np.sum(signed * (1 - 2 * parity)))
        return FragmentResult(corr, original)
    joint = {}
    for outcome, count in counts.items():
        if isinstance(outcome, str):
            bits = outcome.split()[0] if " " in outcome.strip() else outcome.strip()
        else:
            bits = "".join(str(int(b)) for b in outcome)
        if bit_order == "q0_last":
            bits = bits[::-1]
        if len(bits) != expected:
            hint = (
                " Its gate cut measures mid-circuit: build it with "
                "defer_measurements=True and measure every qubit at the end."
                if prep_sorted and len(bits) == width
                else ""
            )
            raise ValueError(
                f"an outcome of {len(bits)} bits ({outcome!r}) for a circuit of "
                f"{expected} qubits.{hint}"
            )
        key = (
            tuple(int(b) for b in bits[:width]),
            _fold_cut_bits(bits[width:], cut_qubits, prep_sorted),
        )
        joint[key] = joint.get(key, 0) + count
    return FragmentResult(_JointCounts(joint, prep_sorted, width), original)


def _outcome_bits(counts, bit_order, expected, width, prep_sorted):
    """The outcomes of ``counts`` as a ``(outcomes, bits)`` uint8 matrix, qubit 0
    first (the first register of a Qiskit key, as :func:`result_from_counts`)."""
    strings = []
    for outcome in counts:
        if isinstance(outcome, str):
            text = outcome.strip()
            strings.append(text.split()[0] if " " in text else text)
        else:
            strings.append("".join(str(int(b)) for b in outcome))
    wrong = next((t for t in strings if len(t) != expected), None)
    if wrong is not None:
        hint = (
            " Its gate cut measures mid-circuit: build it with "
            "defer_measurements=True and measure every qubit at the end."
            if prep_sorted and len(wrong) == width
            else ""
        )
        raise ValueError(
            f"an outcome of {len(wrong)} bits ({wrong!r}) for a circuit of "
            f"{expected} qubits.{hint}"
        )
    bits = np.frombuffer("".join(strings).encode("ascii"), dtype=np.uint8)
    bits = (bits - ord("0")).reshape(len(strings), expected)
    return bits[:, ::-1] if bit_order == "q0_last" else bits


@task(returns=1)
def _result_from_counts_task(subcircuit, counts, bit_order):
    """:func:`result_from_counts` as a task (for :meth:`CutCircuits.reconstruct_counts`),
    keeping only what the list's observable reads."""
    return result_from_counts(subcircuit, counts, bit_order, keep_counts=False)


class CutCircuits(list):
    """The fragment circuits of a :func:`circuit_cutting` cut, plus what recombines
    them::

        subs = qd.circuit_cutting_subcircuits(circuit, cuts, "ZZZZ")
        values = [qd.expectation_value(s, shots=4096) for s in subs]
        value = subs.reconstruct(values)

    ``reconstruct`` returns the expectation value, or the output distribution for
    the list :func:`circuit_cutting_sampling_subcircuits` returns. With
    ``num_samples``, each subcircuit's ``draws`` says how many drawn configurations
    use it: give it shots in proportion (``shots=base * s.draws``) to follow the
    weights; any shots keep the estimate unbiased.
    """

    def __init__(
        self,
        circuits,
        *,
        plan,
        keys,
        readouts,
        mode,
        observable=None,
        draws=None,
        basis=None,
    ):
        super().__init__(circuits)
        #: The Pauli string the circuits measure in (a term, or a group's basis).
        self.basis = (
            basis
            if basis is not None
            else (observable if isinstance(observable, str) else None)
        )
        self.plan = plan
        self.keys = keys
        self.readouts = readouts
        self.mode = mode
        self.observable = observable
        self.draws = draws
        #: Cuts applied, of every kind.
        self.number_cuts = (
            len(plan.gate_cuts) + len(plan.had_cuts) + len(plan.wire_cuts)
        )

    def __repr__(self):
        from Qdislib.core.cutting_algorithms.describe import summary

        return summary(self)

    @property
    def parameters(self):
        """The circuit's free parameters (names), in the order :meth:`bind`
        takes a list of values."""
        return self.plan.parameters

    def bind(self, values):
        """The same subcircuits with the circuit's free parameters bound: a new
        list with the same plan, so the cut is planned once for any number of
        values. ``values``: ``{Parameter or name: value}``, or a list in
        :attr:`parameters` order. Run and reconstruct it as usual."""
        from Qdislib.core.cutting_algorithms.parameters import normalise_values

        bound = CutCircuits(
            _bind_circuits(self, normalise_values(values, self.parameters)),
            plan=self.plan,
            keys=self.keys,
            readouts=self.readouts,
            mode=self.mode,
            observable=self.observable,
            draws=self.draws,
            basis=self.basis,
        )
        return bound

    def reconstruct_counts(self, counts, bit_order="q0_first"):
        """:meth:`reconstruct` from the measurement counts of :meth:`circuits`
        (``defer_measurements=True``, every qubit measured at the end) run with
        your own tools, one counts dict per subcircuit in order; see
        :func:`result_from_counts`."""
        counts = list(counts)
        if len(counts) != len(self):
            raise ValueError(
                f"expected {len(self)} counts, one per subcircuit, got {len(counts)}"
            )
        # one task per subcircuit: converting large counts (a 30-qubit circuit at
        # 400k shots) takes seconds each, in parallel under PyCOMPSs
        return self.reconstruct(
            [_result_from_counts_task(s, c, bit_order) for s, c in zip(self, counts)]
        )

    def circuits(self, software=None, *, defer_measurements=False):
        """The native circuits, one per subcircuit, in this list's order: a plain
        list to run, draw or submit yourself. ``software`` builds them in another
        SDK than the list's (``"qiskit"``, ``"qibo"``, ``"cudaq"``, ``"pennylane"``,
        ``"qilisdk"``); ``defer_measurements`` moves a gate cut's mid-circuit
        measurement to an extra qubit read at the end (see
        :meth:`_FragmentCircuit.to_circuit`). Each circuit's qubits are the
        subcircuit's, in order: :meth:`describe` names them."""
        from Qdislib.utils.backends import _first

        return [_first(s.to_circuit(software, defer_measurements)) for s in self]

    def describe(self, values=None):
        """A table of the subcircuits -- fragment, qubits, the cut terms each
        runs -- and, given their results, each one's values named by qubit and
        the reconstructed value. ``.rows`` holds the data."""
        from Qdislib.core.cutting_algorithms.describe import describe

        return describe(self, values)

    def __getattr__(self, name):
        if name in ("number_components", "weights", "config_of"):
            raise AttributeError(
                f"subcircuit lists have no {name!r} since 2.0.0 (they come from "
                "the circuit_cutting engine): combine the results with "
                "subs.reconstruct(results)."
            )
        raise AttributeError(name)

    def reconstruct(self, results, observable=None):
        """Recombine one result per subcircuit, in this order (values or futures).

        :param observable: Another Pauli string to reconstruct from the same
            results, measured in the same basis as this list (a ``Z``/``I`` string
            for a list made with a ``Z``/``I`` observable or for a distribution
            list). The results must hold more than the list's own correlators: the
            counts of :func:`~Qdislib.sample_subcircuit` run with shots, of
            :func:`result_from_counts` or of a cached run. One set of runs then
            serves many observables, e.g. ``<Z>`` of every qubit.
        """
        results = list(results)
        if len(results) != len(self):
            raise ValueError(
                f"expected {len(self)} results, one per subcircuit, got {len(results)}"
            )
        mode, target = self.mode, self.observable
        if observable is not None:
            self._check_same_basis(observable)
            for result in results:
                corr = result[0] if isinstance(result, tuple) else None
                if isinstance(corr, dict) and not hasattr(corr, "values_at"):
                    raise ValueError(
                        "these results hold only the correlators of the list's own "
                        "observable; to read another one, keep the counts: run the "
                        "subcircuits with sample_subcircuit and shots, or bring "
                        "your own counts with result_from_counts (an exact run "
                        "keeps only what the list's observable needs)"
                    )
            mode, target = "expectation", observable
        return compss_wait_on(
            _reconstruct(
                self.plan,
                self.keys,
                self.readouts,
                results,
                mode,
                target,
                self.draws,
            )
        )

    def _check_same_basis(self, observable):
        """Raise unless ``observable`` is measured in this list's basis."""
        width = self.plan.width
        if not isinstance(observable, str) or len(observable) != width:
            raise ValueError(
                f"observable must be a Pauli string of {width} letters, got "
                f"{observable!r}"
            )
        basis = self.basis or "Z" * width
        for q, letter in enumerate(observable.upper()):
            measured = basis[q] if basis[q] in "XY" else "Z"
            if letter not in ("I", measured):
                raise ValueError(
                    f"qubit {q} is measured in {measured} by these subcircuits; "
                    f"{observable!r} asks for {letter} there. Make a list for it."
                )


@standard_arguments
def circuit_cutting_subcircuits(
    circuit,
    cuts,
    observable=None,
    *,
    _software=None,
    hadamard=False,
    optimized=False,
    joint=True,
    num_samples=None,
    seed=None,
):
    """The fragment circuits of :func:`circuit_cutting`, to run yourself.

    Run each with :func:`~Qdislib.expectation_value` and pass the results to
    ``reconstruct``. A weighted Pauli sum gives a
    :class:`~Qdislib.CombinedCutSubcircuits` of one list per qubit-wise commuting
    group of its terms (see :func:`circuit_cutting`).

    :param num_samples: Draw configurations instead of running them all (see
        :func:`circuit_cutting`); each subcircuit's ``draws`` says how many draws
        use it.
    :param seed: Seed for the ``num_samples`` draws.
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _observable_terms
    from Qdislib.utils.backends import CombinedCutSubcircuits

    software = resolve_software(_software, circuit)
    terms = _observable_terms(default_observable(observable, circuit))
    if len(terms) == 1 and terms[0][0] == 1.0:
        groups = [(terms[0][1], None)]  # one term: its own list, as before
    else:  # qubit-wise commuting groups, one list each (see circuit_cutting)
        groups = _commuting_groups(terms)
    parts = []
    for g, (basis, members) in enumerate(groups):
        masks = None if members is None else [_term_mask(t) for _, t in members]
        plan, draws, keys, readouts, circuits = _prepare_cut(
            circuit,
            cuts,
            basis,
            hadamard,
            optimized,
            joint,
            software,
            False,
            None,
            num_samples,
            None if seed is None else next_seed(seed, g),
            masks=masks,
        )
        parts.append(
            CutCircuits(
                circuits,
                plan=plan,
                keys=keys,
                readouts=readouts,
                mode="expectation",
                observable=(
                    basis
                    if members is None
                    else [(c, mask) for (c, _), mask in zip(members, masks)]
                ),
                draws=draws,
                basis=basis,
            )
        )
    if len(parts) == 1:
        return parts[0]
    return CombinedCutSubcircuits(parts, [1.0] * len(parts))


@standard_arguments
def circuit_cutting_sampling_subcircuits(
    circuit,
    cuts,
    *,
    _software=None,
    hadamard=False,
    optimized=False,
    joint=True,
    num_samples=None,
    seed=None,
):
    """The fragment circuits of :func:`circuit_cutting_sampling`, to sample
    yourself with :func:`~Qdislib.sample_subcircuit`; ``reconstruct`` gives the
    output distribution. Arguments as in :func:`circuit_cutting_subcircuits`."""
    plan, draws, keys, readouts, circuits = _prepare_cut(
        circuit,
        cuts,
        None,
        hadamard,
        optimized,
        joint,
        resolve_software(_software, circuit),
        False,
        None,
        num_samples,
        seed,
    )
    return CutCircuits(
        circuits,
        plan=plan,
        keys=keys,
        readouts=readouts,
        mode="distribution",
        draws=draws,
    )


def cutting_overhead(circuit, cuts, *, hadamard=False, optimized=False, joint=True):
    """What a cut costs, for any cut type or mix, without running anything.

    :param cuts: As in :func:`circuit_cutting`: ``{"gate": ..., "hadamard": ...,
        "wire": ...}`` or a list.
    :return: ``{"subcircuits", "configurations", "fragments", "gamma"}``: the
        circuits to run, the terms of the combined decomposition, the independent
        fragments, and the sampling overhead (for the same worst-case precision,
        the shots grow as ``gamma**2``; see :func:`~Qdislib.shots_for`).
    """
    from Qdislib.utils.gates import resolve_cuts

    plan = _CutPlan(
        circuit, resolve_cuts(circuit, cuts), None, hadamard, optimized, joint
    )
    _warn_cheaper_hadamard(plan)
    return {
        "subcircuits": plan.num_units(),
        "configurations": plan.num_configurations,
        "fragments": len(plan.fragments),
        "gamma": round(plan.gamma(), 9),
    }
