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

"""The IBM Quantum path (gate_cutting(backend=<IBM backend>)) reconstructs like the simulators.

Runs qiskit-ibm-runtime's SamplerV2 in local testing mode on a noiseless Aer
simulator, so it needs no account or network and checks the logic exactly.
"""

import importlib.util
import unittest

from tests import BaseTimedTestCase

RUNTIME = (
    importlib.util.find_spec("qiskit_ibm_runtime") is not None
    and importlib.util.find_spec("qiskit_aer") is not None
)


def _exact(circuit, pauli):
    from qiskit.quantum_info import SparsePauliOp, Statevector

    return float(
        Statevector(circuit).expectation_value(SparsePauliOp(pauli[::-1])).real
    )


@unittest.skipUnless(RUNTIME, "qiskit-ibm-runtime and qiskit-aer required")
class IBMQPUPathTest(BaseTimedTestCase):
    def _run(self, circuit, cut, observable):
        import Qdislib as qd
        from qiskit_aer import AerSimulator
        from qiskit_ibm_runtime import Batch

        device = AerSimulator()
        with Batch(backend=device) as batch:
            ibm = qd.build_backend(
                "IBM_BACKEND",
                qd.BackendConfig(
                    num_qubits=device.num_qubits,
                    extra={"batch": batch, "backend": device},
                ),
            )
            return qd.gate_cutting(circuit, cut, observable, backend=ibm, shots=8000)

    def test_one_cut_with_identity_on_the_cut_qubit(self):
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.ry(0.7, 0)
        circuit.ry(1.2, 1)
        circuit.cz(0, 1)
        circuit.rx(0.5, 0)
        circuit.ry(0.3, 1)
        value = self._run(circuit, ["CZ_3"], "IZ")
        self.assertLess(abs(value - _exact(circuit, "IZ")), 0.06)

    def test_two_cuts_on_one_pair(self):
        from qiskit import QuantumCircuit

        circuit = QuantumCircuit(2)
        circuit.ry(0.7, 0)
        circuit.rx(1.1, 1)
        circuit.cz(0, 1)
        circuit.ry(0.9, 0)
        circuit.rx(0.6, 1)
        circuit.cz(0, 1)
        circuit.ry(0.4, 0)
        circuit.rx(0.8, 1)
        value = self._run(circuit, ["CZ_3", "CZ_6"], "ZZ")
        self.assertLess(abs(value - _exact(circuit, "ZZ")), 0.08)


if __name__ == "__main__":
    unittest.main()
