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

"""Regression test: Qiskit multi-register circuits must use global qubit indices."""

from qiskit import QuantumCircuit, QuantumRegister

from tests import BaseTimedTestCase


class MultiRegisterDagTest(BaseTimedTestCase):

    def test_multi_register_global_indices(self):
        from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag

        a = QuantumRegister(2, "a")
        b = QuantumRegister(2, "b")
        qc = QuantumCircuit(a, b)
        qc.cx(a[0], b[0])  # global qubits 0 and 2
        qc.x(b[1])  # global qubit 3

        dag = circuit_qiskit_to_dag(qc)
        qubits_by_node = {n: tuple(d["qubits"]) for n, d in dag.nodes(data=True)}

        cx = next(n for n in qubits_by_node if n.startswith("CX"))
        x = next(n for n in qubits_by_node if n.startswith("X"))

        # Global indices: cx acts on (0, 2); x acts on (3,).
        self.assertEqual(
            qubits_by_node[cx],
            (0, 2),
            f"cx should act on global qubits (0, 2), got {qubits_by_node[cx]}",
        )
        self.assertEqual(
            qubits_by_node[x],
            (3,),
            f"x should act on global qubit (3,), got {qubits_by_node[x]}",
        )
        # And no spurious self-loop introduced by colliding per-register indices.
        self.assertEqual(
            len(set(qubits_by_node[cx])),
            2,
            "cx collapsed to a self-loop (per-register index bug)",
        )
