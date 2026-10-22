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

import itertools

import numpy as np

from Qdislib.utils.runtime import task, compss_wait_on, COLLECTION_IN

__all__ = [
    "gate_wire_cutting",
    "circuit_cutting",
    "circuit_cutting_subcircuits",
    "circuit_cutting_sampling",
    "circuit_cutting_sampling_subcircuits",
    "circuit_cutting_overhead",
    "CutCircuits",
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
    """Reconstruct ``<observable>`` from a circuit cut by BOTH gate and wire cuts.

    :func:`circuit_cutting` generalises this to any combination of cuts, with the
    output distribution and three-step forms; with ``hadamard``, ``optimized`` or
    ``joint`` set, this call is forwarded to it.

    :param circuit: A qibo / Qiskit / CUDA-Q circuit (or a Qdislib DAG).
    :param gates_cut: Two-qubit gate node names to cut (CZ/CX), as in
        :func:`gate_cutting`.
    :param wire_cut: Wire cuts as gate-pair edges ``[(source, target), ...]``, as in
        :func:`wire_cutting` (each is expanded to one cut per shared qubit).
    :param observable: A Z/I (or general Pauli) string over the circuit's qubits.
    :param software: Fragment-execution backend (``"qibo"``/``"qiskit"``/``"cudaq"``/
        ``"pennylane"``).
    :param shots: Samples per fragment expectation.
    :param method: Backend evaluation method (forwarded to the gate-cut executor).
    :param sync: If ``True`` block and return the value; if ``False`` return the list
        of per-configuration futures (advanced overlap of several cuts).
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :param hadamard: Cut the gates with the Hadamard test (see :func:`circuit_cutting`).
    :param optimized: Run each fragment once per setting of its own cuts.
    :param joint: Joint wire measurement / joint Hadamard cut (see
        :func:`circuit_cutting`).
    :return: The reconstructed expectation value.
    """
    if hadamard or optimized or joint is not None:
        return circuit_cutting(
            circuit,
            {"gate": gates_cut or [], "wire": wire_cut or []},
            observable,
            software=software,
            shots=shots,
            hadamard=hadamard,
            optimized=optimized,
            joint=joint,
            sync=sync,
            qpu=qpu,
            qpu_dict=qpu_dict,
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
    components_before = nx.number_connected_components(dag.to_undirected())

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
    dag = _generate_cut(dag, gates_cut) if gates_cut else _remove_red_edges(dag.copy())

    # Same reasoning as gate_cutting._validate_gate_cuts_partition /
    # wire_cutting._validate_wire_cuts_partition: a cut set that doesn't increase
    # the connected-component count beyond the uncut circuit's baseline pays the
    # combined gate+wire quasiprobability cost for no separation benefit.
    if (gates_cut or wire_cuts) and nx.number_connected_components(
        dag.to_undirected()
    ) <= components_before:
        raise ValueError(
            f"The given gate cut(s) {gates_cut} and wire cut(s) {wire_cuts} do not "
            "partition the circuit into independent fragments: some other, uncut "
            "path still connects the pieces they're meant to separate. Either cut "
            "wires/gates that are the sole connectivity between the pieces you "
            "want to separate, or pass empty cut lists to run the circuit directly."
        )

    n_gate, n_wire = len(gates_cut), len(wire_cuts)
    total_qubits = num_qubits + n_wire  # fresh wire lines are numbered above the rest

    # Enumerate the PRODUCT of the two QPD configuration spaces; each configuration
    # applies both substitutions to a fresh DAG copy, splits into independent fragments,
    # fans them out, and contributes sign * prod(fragments).
    results = []
    for gate_idx in range(6**n_gate):
        for wire_idx in range(8**n_wire):
            config_dag = _remove_red_edges(dag.copy())
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
                for component in nx.connected_components(config_dag.to_undirected())
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
        gate = list(cuts.get("gate") or [])
        had = list(cuts.get("hadamard") or [])
        wire = [tuple(edge) for edge in cuts.get("wire") or []]
    else:
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

    def __init__(self, circuit, cuts, observable, hadamard, optimized, joint):
        import networkx as nx

        from Qdislib.utils.backends import convert_to_dag, change_basis_dag
        from Qdislib.utils.graph_qibo import _remove_red_edges
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _generate_cut,
            _expand_general_gate_cuts,
        )
        from Qdislib.core.cutting_algorithms.hadamard_cut import (
            _generate_cut_hadamard,
            _cut_groups,
            _cut_settings,
            local_unitary_decomposition,
            _num_circuit_qubits,
        )
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _expand_wire_cuts,
            _place_wire_cuts,
            _wire_cut_line_maps,
        )

        gate_cuts, had_cuts, wire_edges = _normalise_cuts(cuts, hadamard)
        overlap = set(gate_cuts) & set(had_cuts)
        if overlap:
            raise ValueError(f"Gate(s) {sorted(overlap)} cut twice.")
        self.optimized = bool(optimized)
        # joint=None: each cut's own default -- joint wire measurement (optimized
        # only), Hadamard gates cut one by one.
        self.joint_wires = True if joint is None else bool(joint)
        joint_gates = False if joint is None else bool(joint)

        dag = convert_to_dag(circuit)
        if observable and set(observable) & {"X", "Y"}:
            dag = change_basis_dag(dag, observable)
        self.num_qubits = 1 + max(
            (q for _, d in dag.nodes(data=True) for q in d.get("qubits") or ()),
            default=0,
        )
        self.width = max(self.num_qubits, _num_circuit_qubits(circuit))
        components_before = nx.number_connected_components(dag.to_undirected())

        # Read the Hadamard decompositions before any surgery touches the gates.
        decomps = [
            local_unitary_decomposition(
                (dag.nodes[name].get("gate") or "").upper(),
                dag.nodes[name].get("parameters") or (),
            )
            for name in had_cuts
        ]
        if gate_cuts:
            dag, gate_cuts = _expand_general_gate_cuts(dag, gate_cuts)
        # Wires first, as gate_wire_cutting does: a wire cut right after a cut gate
        # still finds the gate in place.
        wire_cuts = _expand_wire_cuts(dag, wire_edges) if wire_edges else []
        placed = _place_wire_cuts(dag, wire_cuts, self.num_qubits)
        self.cut_source_index, self.prepared_to_orig = _wire_cut_line_maps(
            placed, self.num_qubits
        )
        if gate_cuts:
            dag = _generate_cut(dag, gate_cuts)
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

        if (gate_cuts or had_cuts or wire_cuts) and nx.number_connected_components(
            graph.to_undirected()
        ) <= components_before:
            raise ValueError(
                f"The cuts {{gate: {gate_cuts}, hadamard: {had_cuts}, wire: "
                f"{wire_cuts}}} do not partition the circuit into independent "
                "fragments: some other, uncut path still connects the pieces they're "
                "meant to separate."
            )
        self.graph = graph
        self.gate_cuts, self.had_cuts, self.wire_cuts = gate_cuts, had_cuts, wire_cuts
        self.ancilla_base = self.num_qubits + len(wire_cuts) + 1

        components = list(nx.connected_components(graph.to_undirected()))
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

        # The gate configurations: measure-and-prepare terms x Hadamard settings.
        mp = [
            (
                digits,
                (-1) ** sum(d in (3, 5) for d in digits) / 2 ** len(gate_cuts),
            )
            for digits in itertools.product(range(6), repeat=len(gate_cuts))
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
        self.gate_configs = [
            (digits, config, w_mp * w_had)
            for digits, w_mp in mp
            for config, w_had in had
        ]
        self.wire_configs = [
            (
                digits,
                (-1) ** sum(d in (3, 5, 7) for d in digits) / 2 ** len(wire_cuts),
            )
            for digits in itertools.product(range(8), repeat=len(wire_cuts))
        ]

    @property
    def num_configurations(self):
        return len(self.gate_configs) * len(self.wire_configs)

    def gamma(self):
        """The sampling overhead: the product of every cut's 1-norm."""
        from Qdislib.core.cutting_algorithms.hadamard_cut import phi_factor

        total = 3.0 ** len(self.gate_cuts) * 4.0 ** len(self.wire_cuts)
        for _, terms in self.groups:
            total *= phi_factor(terms)
        return total

    def gate_legs(self, f, digits, config):
        """Fragment ``f``'s part of a gate configuration."""
        from Qdislib.core.cutting_algorithms.gate_cutting import _GATE_SIDE_OPS

        fragment = self.fragments[f]
        return (
            tuple(_GATE_SIDE_OPS[s][digits[k - 1]] for k, s in fragment["gate_sides"]),
            tuple(config[g] for g in fragment["had_groups"]),
        )

    def units(self):
        """``[(key, f, gate_config, prep_states, setting), ...]``: the fragment
        circuits to run. ``key`` is how the reconstruction finds them."""
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _PREP_STATES,
            _measurement_settings,
        )

        units = []
        if not self.optimized:
            for c, ((digits, config, _), (wire, _)) in enumerate(
                itertools.product(self.gate_configs, self.wire_configs)
            ):
                for f, fragment in enumerate(self.fragments):
                    setting = tuple(
                        _WIRE_DIGIT_BASIS[wire[i - 1]] for i in fragment["measures"]
                    )
                    states = tuple(
                        _WIRE_DIGIT_STATE[wire[i - 1]] for i in fragment["preps"]
                    )
                    units.append(((c, f), f, (digits, config), states, setting))
            return units
        for f, fragment in enumerate(self.fragments):
            legs = {}  # distinct gate legs -> a configuration that has them
            for digits, config, _ in self.gate_configs:
                legs.setdefault(self.gate_legs(f, digits, config), (digits, config))
            for leg, gate_config in legs.items():
                for states in itertools.product(
                    _PREP_STATES, repeat=len(fragment["preps"])
                ):
                    for setting in _measurement_settings(
                        len(fragment["measures"]), self.joint_wires
                    ):
                        units.append(
                            ((f, leg, states), f, gate_config, states, setting)
                        )
        return units

    def build(self, f, gate_config, states, setting):
        """The DAG of fragment ``f`` for one set of legs, and its readout."""
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _GATE_SIDE_OPS,
            _apply_gate_side,
        )
        from Qdislib.core.cutting_algorithms.hadamard_cut import _config_assignments
        from Qdislib.core.cutting_algorithms.wire_cutting import (
            _apply_measurement_setting,
        )

        fragment = self.fragments[f]
        digits, config = gate_config
        sub = self.graph.subgraph(fragment["nodes"]).copy()
        for k, side in fragment["gate_sides"]:
            node = f"GATE{side}_{k}"
            _apply_gate_side(
                sub, node, _GATE_SIDE_OPS[side][digits[k - 1]], f"{node}_rz"
            )
        for cut, state in zip(fragment["preps"], states):
            _prepare_wire_state(sub, cut, state)
        readout = _apply_measurement_setting(
            sub, fragment["measures"], fragment["cut_qubits"], setting
        )
        if fragment["had_groups"]:
            assignments = _config_assignments(self.groups, config)
            _insert_hadamard_gadgets(
                sub,
                {n: a for n, a in assignments.items() if n in sub},
                fragment["ancilla_of"],
            )
        return sub, readout


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
    graph.add_node(flip, gate="X", qubits=graph.nodes[first]["qubits"], parameters=())
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
            append(ancilla, f"{node}_open", "H")
            if g == 1:
                append(ancilla, f"{node}_open_S", "S")
        theta_i, phi_i, lam_i, _ = _u3_angles(u_i)
        theta_m, phi_m, lam_m, gphase = _u3_angles(
            np.asarray(u_j) @ np.asarray(u_i).conj().T
        )
        preds, succs = list(dag.pred[node]), list(dag.succ[node])
        dag.remove_node(node)
        u3, cu3 = f"{node}_U3", f"{node}_CU3"
        dag.add_node(u3, gate="U3", qubits=(qubit,), parameters=(theta_i, phi_i, lam_i))
        dag.add_node(
            cu3, gate="CU3", qubits=(ancilla, qubit), parameters=(theta_m, phi_m, lam_m)
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
            append(ancilla, f"{node}_U1", "U1", (gphase,))
        if last_marker[key] == node:
            append(ancilla, f"{node}_close", "H")


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
            dag.nodes[node].update(gate="CNOT", qubits=(qubit, fresh + deferred))
            deferred += 1
    renumbered, n = _update_qubits(dag)
    if software == "qiskit":
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


class _FragmentCircuit:
    """One fragment circuit of a :func:`circuit_cutting` cut (three-step workflow).

    Run it with :func:`~Qdislib.sample_subcircuit` or
    :func:`~Qdislib.expectation_value` (either works, both return the fragment's
    ``(corr, original)`` Z correlators); ``shots=None`` is exact (statevector, qibo or Qiskit).
    Picklable: a DAG plus metadata.
    """

    _is_optimized_wire_measure = True  # expectation_value calls evaluate()

    def __init__(self, dag, software, qpu=False, qpu_dict=None):
        self.dag = dag
        self.software = software
        self.qpu = qpu
        self.qpu_dict = qpu_dict

    def sample(self, shots=1024, seed=None):
        if shots is None:
            return _exact_fragment(self.dag, self.software)
        from Qdislib.core.cutting_algorithms.gate_cutting import (
            _sample_gate_subcircuit,
        )

        return _sample_gate_subcircuit(
            self.dag, self.software, shots, self.qpu, self.qpu_dict, seed
        )

    def evaluate(self, shots):
        return self.sample(shots)

    def __repr__(self):
        qubits = sorted(
            {q for _, d in self.dag.nodes(data=True) for q in d.get("qubits") or ()}
        )
        return f"_FragmentCircuit(qubits={qubits}, gates={self.dag.number_of_nodes()})"


def _cut_z_correlators(plan, keys, readouts, results, masks):
    """``<Z^A>`` of the uncut circuit for each global mask ``A`` in ``masks``.

    ``results`` holds one ``(corr, original)`` per fragment circuit (``keys`` and
    ``readouts`` parallel to it). Sums ``weight * prod_fragment value`` over every
    configuration; a fragment's value reads the circuit with its legs at the mask
    whose output bits come from ``A``, whose cut-wire bits are the needed Pauli's
    readout bits, and whose ancilla bits are always set.
    """
    from Qdislib.core.cutting_algorithms.wire_cutting import (
        _SETTING_PAULI,
        _WIRE_PREP,
    )

    masks = np.asarray(masks, dtype=np.int64)
    by_key = {}
    for index, key in enumerate(keys):
        by_key.setdefault(key, []).append(index)

    prepared = []  # per circuit: corr array, output index over masks, cut positions
    for index, (corr, original) in enumerate(results):
        f = keys[index][0] if plan.optimized else keys[index][1]
        fragment = plan.fragments[f]
        if isinstance(corr, dict):
            corr = np.array([corr[m] for m in range(len(corr))])
        local = np.zeros(len(masks), dtype=np.int64)
        for i, qubit in enumerate(original):
            if qubit >= plan.ancilla_base:
                local |= 1 << i
            elif plan.cut_source_index.get(qubit) not in fragment["measures"]:
                output = plan.prepared_to_orig.get(qubit, qubit)
                local |= ((masks >> output) & 1) << i
        positions = [
            original.index(fragment["cut_qubits"][c]) for c in fragment["measures"]
        ]
        prepared.append((np.asarray(corr), local, positions))

    cache = {}

    def read(key, letters):
        """The fragment circuits under ``key``, read for the Pauli ``letters`` on
        its cut wires (averaged when several settings read it)."""
        if (key, letters) not in cache:
            values = []
            for index in by_key[key]:
                corr, local, positions = prepared[index]
                for pauli, sign, bits in readouts[index]:
                    if tuple(pauli) == letters:
                        index_bits = sum(1 << positions[b] for b in bits)
                        values.append(sign * corr[local | index_bits])
            if not values:
                raise RuntimeError(f"no circuit reads {letters} for {key}")
            cache[(key, letters)] = np.mean(values, axis=0)
        return cache[(key, letters)]

    total = np.zeros(len(masks))
    for c, ((digits, config, w_gate), (wire, w_wire)) in enumerate(
        itertools.product(plan.gate_configs, plan.wire_configs)
    ):
        contribution = np.full(len(masks), w_gate * w_wire)
        for f, fragment in enumerate(plan.fragments):
            letters = tuple(_SETTING_PAULI[wire[i - 1]] for i in fragment["measures"])
            if not plan.optimized:
                contribution *= read((c, f), letters)
                continue
            leg = plan.gate_legs(f, digits, config)
            value = 0.0
            for choice in itertools.product(
                *(_WIRE_PREP[wire[i - 1]] for i in fragment["preps"])
            ):
                coefficient = float(np.prod([weight for weight, _ in choice]))
                states = tuple(state for _, state in choice)
                value = value + coefficient * read((f, leg, states), letters)
            contribution *= value
        total += contribution
    return total


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


@task(returns=1, results=COLLECTION_IN)
def _reconstruct_cut_task(plan, keys, readouts, results, mode, observable):
    """Expectation (``mode="expectation"``) or distribution from the fragment
    results, as a task fed straight by the fragment futures."""
    if mode == "expectation":
        mask = _observable_mask(observable, plan.width)
        return float(_cut_z_correlators(plan, keys, readouts, results, [mask])[0])
    z = _cut_z_correlators(plan, keys, readouts, results, np.arange(2**plan.width))
    return _distribution(z, plan.width)


def _prepare_cut(
    circuit, cuts, observable, hadamard, optimized, joint, software, qpu, qpu_dict
):
    """The plan plus one :class:`_FragmentCircuit` per fragment circuit to run."""
    if software is None:
        from Qdislib.utils.circuit_cache import _detect_software

        software = _detect_software(circuit) or "qibo"
    plan = _CutPlan(circuit, cuts, observable, hadamard, optimized, joint)
    keys, readouts, circuits = [], [], []
    for key, f, gate_config, states, setting in plan.units():
        dag, readout = plan.build(f, gate_config, states, setting)
        keys.append(key)
        readouts.append(readout)
        circuits.append(_FragmentCircuit(dag, software, qpu, qpu_dict))
    return plan, keys, readouts, circuits


def _run(circuits, shots, seed):
    return [
        circuit.sample(shots, None if seed is None else seed + n)
        for n, circuit in enumerate(circuits)
    ]


def _check_exact(shots, software):
    if shots is None and software not in (None, "qibo", "qiskit"):
        raise ValueError("shots=None (exact) needs software='qibo' or 'qiskit'")


def circuit_cutting(
    circuit,
    cuts,
    observable=None,
    software=None,
    shots=1024,
    hadamard=False,
    optimized=False,
    joint=None,
    sync=True,
    qpu=False,
    qpu_dict=None,
    seed=None,
):
    """Reconstruct ``<observable>`` from a circuit cut by any combination of cuts.

    Measure-and-prepare gate cuts, Hadamard gate cuts and wire cuts can be mixed
    freely, in any number; their quasiprobability decompositions tensor. With
    ``optimized=True`` each fragment runs once per distinct setting of the cuts that
    touch it, so the number of circuits grows as a *sum* over fragments, not as a
    product over cuts (see :func:`circuit_cutting_overhead`).

    :param circuit: A qibo / Qiskit / CUDA-Q / PennyLane circuit (or a Qdislib DAG).
    :param cuts: ``{"gate": [...], "hadamard": [...], "wire": [...]}`` (any subset;
        ``find_cut(..., mixed=True)`` returns one), or a list mixing gate names and
        ``(source, target)`` wire edges.
    :param observable: A Pauli string (``X``/``Y`` are rotated into ``Z``), or a
        weighted Pauli sum ``[(coefficient, string), ...]``; ``None`` is all ``Z``.
    :param software: Fragment backend; ``None`` uses the circuit's own SDK.
    :param shots: Shots per fragment circuit; ``None`` is exact (statevector, on
        qibo or Qiskit).
    :param hadamard: Cut the ``"gate"`` cuts (or the gate names of a list) with the
        Hadamard test (:func:`gate_cutting_hadamard`) instead of measure-and-prepare.
    :param optimized: Run each fragment once per distinct setting of its own cuts.
    :param joint: ``None`` keeps each cut's default (wires measured jointly when
        ``optimized``, Hadamard gates cut one by one); ``True`` / ``False`` turns
        both on / off (:func:`wire_cutting` and :func:`gate_cutting_hadamard`).
    :param sync: ``False`` returns the value as a future.
    :param qpu: Route fragment execution to a QPU.
    :param qpu_dict: QPU connection settings.
    :param seed: Seed for the fragment sampling (``seed + i`` per circuit).
    :return: The reconstructed expectation value.
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _observable_terms

    _check_exact(shots, software)
    terms = _observable_terms(observable) if observable is not None else [(1.0, None)]
    values = []
    for coefficient, term in terms:
        plan, keys, readouts, circuits = _prepare_cut(
            circuit, cuts, term, hadamard, optimized, joint, software, qpu, qpu_dict
        )
        results = _run(circuits, shots, seed)
        values.append(
            (
                coefficient,
                _reconstruct_cut_task(
                    plan, keys, readouts, results, "expectation", term
                ),
            )
        )
    if not sync and len(values) == 1 and values[0][0] == 1.0:
        return values[0][1]
    return sum(c * compss_wait_on(v) for c, v in values)


def circuit_cutting_sampling(
    circuit,
    cuts,
    software=None,
    shots=1024,
    hadamard=False,
    optimized=False,
    joint=None,
    sync=True,
    qpu=False,
    qpu_dict=None,
    seed=None,
):
    """Reconstruct the output **probability distribution** of a circuit cut by any
    combination of cuts -- the distribution counterpart of :func:`circuit_cutting`,
    from the same fragment circuits. Arguments as there.

    :return: ``{bitstring: probability}`` (qubit-0-first bitstrings), the signed
        estimate abs-normalised like :func:`gate_cutting_sampling`.
    """
    _check_exact(shots, software)
    plan, keys, readouts, circuits = _prepare_cut(
        circuit, cuts, None, hadamard, optimized, joint, software, qpu, qpu_dict
    )
    future = _reconstruct_cut_task(
        plan, keys, readouts, _run(circuits, shots, seed), "distribution", None
    )
    return compss_wait_on(future) if sync else future


class CutCircuits(list):
    """The fragment circuits of a :func:`circuit_cutting` cut, plus what recombines
    them::

        subs = qd.circuit_cutting_subcircuits(circuit, cuts, "ZZZZ")
        values = [qd.expectation_value(s, shots=4096) for s in subs]
        value = subs.reconstruct(values)

    ``reconstruct`` returns the expectation value, or the output distribution for
    the list :func:`circuit_cutting_sampling_subcircuits` returns.
    """

    def __init__(self, circuits, *, plan, keys, readouts, mode, observable=None):
        super().__init__(circuits)
        self.plan = plan
        self.keys = keys
        self.readouts = readouts
        self.mode = mode
        self.observable = observable
        #: Cuts applied, of every kind.
        self.number_cuts = (
            len(plan.gate_cuts) + len(plan.had_cuts) + len(plan.wire_cuts)
        )

    def reconstruct(self, results):
        """Recombine one ``(corr, original)`` result per circuit, in this order."""
        results = compss_wait_on(list(results))
        if len(results) != len(self):
            raise ValueError(
                f"expected {len(self)} results, one per subcircuit, got {len(results)}"
            )
        return compss_wait_on(
            _reconstruct_cut_task(
                self.plan, self.keys, self.readouts, results, self.mode, self.observable
            )
        )


def circuit_cutting_subcircuits(
    circuit,
    cuts,
    observable=None,
    software=None,
    hadamard=False,
    optimized=False,
    joint=None,
    qpu=False,
    qpu_dict=None,
):
    """The fragment circuits of :func:`circuit_cutting` (three-step workflow).

    Run each with :func:`~Qdislib.expectation_value` (or
    :func:`~Qdislib.sample_subcircuit`) and pass the results to ``reconstruct``. A
    weighted Pauli sum gives a :class:`~Qdislib.CombinedCutSubcircuits` of one
    :class:`CutCircuits` per term. Arguments as in :func:`circuit_cutting`.
    """
    from Qdislib.core.cutting_algorithms.hadamard_cut import _observable_terms
    from Qdislib.utils.backends import CombinedCutSubcircuits

    terms = _observable_terms(observable) if observable is not None else [(1.0, None)]
    parts = []
    for _, term in terms:
        plan, keys, readouts, circuits = _prepare_cut(
            circuit, cuts, term, hadamard, optimized, joint, software, qpu, qpu_dict
        )
        parts.append(
            CutCircuits(
                circuits,
                plan=plan,
                keys=keys,
                readouts=readouts,
                mode="expectation",
                observable=term,
            )
        )
    if len(parts) == 1 and terms[0][0] == 1.0:
        return parts[0]
    return CombinedCutSubcircuits(parts, [c for c, _ in terms])


def circuit_cutting_sampling_subcircuits(
    circuit,
    cuts,
    software=None,
    hadamard=False,
    optimized=False,
    joint=None,
    qpu=False,
    qpu_dict=None,
):
    """The fragment circuits of :func:`circuit_cutting_sampling` (three-step):
    sample each with :func:`~Qdislib.sample_subcircuit`, then ``reconstruct`` gives
    the output distribution. Arguments as in :func:`circuit_cutting`."""
    plan, keys, readouts, circuits = _prepare_cut(
        circuit, cuts, None, hadamard, optimized, joint, software, qpu, qpu_dict
    )
    return CutCircuits(
        circuits, plan=plan, keys=keys, readouts=readouts, mode="distribution"
    )


def circuit_cutting_overhead(
    circuit, cuts, hadamard=False, optimized=False, joint=None
):
    """What a cut costs, without running anything.

    :return: ``{"subcircuits", "configurations", "fragments", "gamma"}``: the
        fragment circuits to run, the terms of the combined decomposition, the
        independent fragments, and the sampling overhead (shots grow as
        ``gamma**2`` for a given precision).
    """
    plan = _CutPlan(circuit, cuts, None, hadamard, optimized, joint)
    if optimized:
        subcircuits = len(plan.units())
    else:  # every fragment of every configuration
        subcircuits = plan.num_configurations * len(plan.fragments)
    return {
        "subcircuits": subcircuits,
        "configurations": plan.num_configurations,
        "fragments": len(plan.fragments),
        "gamma": round(plan.gamma(), 9),
    }
