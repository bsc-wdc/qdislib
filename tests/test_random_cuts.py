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

"""Random circuits and random cuts against an independent statevector.

Hand-written cases list their cuts in circuit order, on gates that rarely repeat a
qubit pair, and so missed wrong values that only appear otherwise. Here every kind of
cut, in any order, has to give the value Qiskit's statevector gives, or be refused.
"""

import importlib.util
import random
import unittest
import warnings

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

#: Fixed, so a failure can be reproduced. Each circuit runs every kind of cut, twice.
SEED = 5
CIRCUITS = 60
TOLERANCE = 1e-7


def _random_circuit(rng):
    """A few rotations and two-qubit gates on neighbouring lines, so pairs repeat, with
    the circuit's two-qubit gates and the wires between neighbouring gates by name."""
    from qiskit import QuantumCircuit

    width = rng.choice([2, 3, 3, 4])
    circuit = QuantumCircuit(width)
    two_qubit, wires, last = [], [], {}
    for index in range(1, rng.randint(4, 9) + 1):
        if rng.random() < 0.5:
            line = rng.randrange(width)
            kind = rng.choice(["ry", "rx", "h"])
            if kind == "h":
                circuit.h(line)
            else:
                getattr(circuit, kind)(rng.uniform(0.1, 3.0), line)
            name, lines = f"{kind.upper()}_{index}", [line]
        else:
            first = rng.randrange(width - 1)
            kind = rng.choice(["cz", "cz", "cx"])
            getattr(circuit, kind)(first, first + 1)
            name, lines = f"{kind.upper()}_{index}", [first, first + 1]
            two_qubit.append(name)
        for line in lines:
            if line in last:
                wires.append((last[line], name))
            last[line] = name
    return circuit, two_qubit, wires


def _cut_sets(rng, two_qubit, wires):
    """Each kind of cut on the same gates, including out of circuit order."""
    picked = rng.sample(two_qubit, min(len(two_qubit), rng.randint(1, 3)))
    sets = [
        {"gate": picked},
        {"gate": picked[::-1]},
        {"hadamard": picked},
        {"hadamard": picked[::-1]},
        {"gate": picked[:1], "hadamard": picked[1:]},
        {"gate": picked[1:], "hadamard": picked[:1]},
    ]
    if wires:
        sets.append({"gate": picked[:1], "wire": rng.sample(wires, 1)})
    return [{kind: names for kind, names in cuts.items() if names} for cuts in sets]


@unittest.skipUnless(QISKIT, "qiskit required")
class RandomCutsTest(BaseTimedTestCase):
    def test_every_cut_gives_the_circuits_value_or_is_refused(self):
        from qiskit.quantum_info import SparsePauliOp, Statevector

        import Qdislib as qd
        from Qdislib.utils.exceptions import InvalidCutError

        rng = random.Random(SEED)
        compared, wrong = 0, []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for number in range(CIRCUITS):
                circuit, two_qubit, wires = _random_circuit(rng)
                if not two_qubit:
                    continue
                observable = "".join(rng.choice("ZZZI") for _ in circuit.qubits)
                if set(observable) == {"I"}:
                    observable = "Z" * circuit.num_qubits
                exact = float(
                    Statevector(circuit)
                    .expectation_value(SparsePauliOp(observable[::-1]))
                    .real
                )
                for cuts in _cut_sets(rng, two_qubit, wires):
                    for joint in (True, False):
                        try:
                            value = qd.circuit_cutting(
                                circuit, cuts, observable, shots=None, joint=joint
                            )
                        except InvalidCutError:
                            continue  # a cut that does not separate is refused
                        compared += 1
                        if abs(value - exact) > TOLERANCE:
                            wrong.append((number, cuts, joint, value, exact))
        self.assertGreater(compared, 200)  # the test has to keep testing something
        self.assertEqual(wrong, [], f"{len(wrong)} of {compared} gave a wrong value")


if __name__ == "__main__":
    unittest.main()
