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

"""The one place a circuit becomes a Qdislib DAG.

Each SDK's converter only lists its circuit's operations, ``(label, name, qubits,
params)`` in circuit order; :func:`ops_to_dag` does the rest the same way for every
SDK: canonical gate names and parameters (:mod:`Qdislib.utils.gates`), the node
labels, the dependency edges and the ``Observable I`` markers.

Edges, per qubit of a gate:

- **blue**: from the last earlier gate on that qubit (the circuit's wires);
- **red**: from the last earlier gate on that qubit not already linked to this one
  (a second, transitive dependency the cutting code relies on, then drops).

Each qubit's history is kept, so building the DAG is linear in the gates.
"""

import networkx as nx

from Qdislib.utils.gates import canonical

__all__ = ["ops_to_dag"]


def ops_to_dag(ops, obs_I=None, sdk=None):
    """A DAG from an SDK's operations.

    :param ops: ``(label, name, qubits, params)`` per operation, in circuit order;
        ``label`` is the node's name (what a cut names), ``name`` the SDK's gate name.
    :param obs_I: Optional per-qubit observable; an ``"I"`` appends an
        ``Observable I`` marker after the qubit's last gate (internal use).
    :param sdk: The SDK the names come from (see :func:`~Qdislib.utils.gates.canonical`).
    :return: A :class:`networkx.DiGraph`; each node has ``gate`` (canonical name),
        ``qubits`` and ``parameters``.
    """
    dag = nx.DiGraph()
    history = {}
    for label, name, qubits, params in ops:
        qubits = tuple(qubits)
        gate, params = canonical(name, params, sdk)
        dag.add_node(label, gate=gate, qubits=qubits, parameters=params)
        for qubit in qubits:
            earlier = history.get(qubit)
            if earlier:
                dag.add_edge(earlier[-1], label, color="blue")
        for qubit in qubits:
            for earlier in reversed(history.get(qubit, ())):
                if earlier != label and not dag.has_edge(earlier, label):
                    dag.add_edge(earlier, label, color="red")
                    break
        for qubit in qubits:
            history.setdefault(qubit, []).append(label)
    for qubit, pauli in enumerate(obs_I or ()):
        # After the LAST gate on the qubit, whatever its width: requiring a sink
        # acting on (qubit,) alone dropped the marker when a qubit ended on a CZ.
        if pauli == "I" and history.get(qubit):
            marker = f"OBSI_{qubit}"
            dag.add_node(marker, gate="Observable I", qubits=(qubit,), parameters=())
            dag.add_edge(history[qubit][-1], marker, color="blue")
    return dag
