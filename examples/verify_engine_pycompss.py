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

"""The cutting engine's options under real PyCOMPSs, against the exact answer.

What the serial tests cannot see: a backend object reaching a task, groups of
circuits in one task, the per-fragment reconstruction tasks, free parameters and
error bars. Prints ``RESULT: PASS`` or ``RESULT: FAIL``.

    runcompss --pythonpath=. examples/verify_engine_pycompss.py
"""

import os
import sys

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp, Statevector

import Qdislib as qd
from Qdislib.utils.runtime import compss_wait_on, pycompss_available

SHOTS = int(os.environ.get("QDISLIB_SHOTS", 20000))
TOL = float(os.environ.get("QDISLIB_TVD_TOL", 0.1))


def circuit(angle=0.7):
    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(angle, 2)
    qc.rx(0.3, 1)
    qc.cx(2, 3)
    return qc


def exact(qc, qubits=None):
    probs = Statevector(qc).probabilities_dict()
    out = {}
    for bits, p in probs.items():
        bits = bits[::-1]
        key = bits if qubits is None else "".join(bits[q] for q in qubits)
        out[key] = out.get(key, 0.0) + p
    return out


def tvd(got, want):
    return 0.5 * sum(
        abs(got.get(k, 0.0) - want.get(k, 0.0)) for k in set(got) | set(want)
    )


def main():
    print("pycompss:", pycompss_available, flush=True)
    qc = circuit()
    aer = qd.build_backend("AER_BACKEND", qd.BackendConfig(num_qubits=4))
    cuts = {
        "gate": {"gate": ["CZ_4"]},
        "hadamard": {"hadamard": ["CZ_4"]},
        "wire": {"wire": [("CZ_4", "RY_5")]},
    }
    futures = {}
    for name, cut in cuts.items():
        run = lambda **k: qd.circuit_cutting_sampling(
            qc, cut, shots=SHOTS, seed=1, sync=False, **k
        )
        futures[f"{name}: backend object"] = (run(backend=aer), exact(qc))
        futures[f"{name}: groups in one task"] = (run(circuits_per_task=4), exact(qc))
        futures[f"{name}: marginal"] = (run(qubits=[3, 1]), exact(qc, [3, 1]))
        futures[f"{name}: error bars"] = (run(error=10), exact(qc))

    failures = []
    for name, (future, want) in futures.items():
        got = compss_wait_on(future)
        distance = tvd(got, want)
        print(f"{name}: TVD {distance:.4f}", flush=True)
        if distance > TOL:
            failures.append(name)

    # free parameters: one cut, two sets of values
    a = Parameter("a")
    template = circuit(a)
    values = [[0.7], [1.9]]
    runs = compss_wait_on(
        qd.circuit_cutting_sampling(
            template,
            {"hadamard": ["CZ_4"]},
            shots=SHOTS,
            seed=2,
            sync=False,
            parameter_values=values,
        )
    )
    for value, got in zip(values, runs):
        distance = tvd(got, exact(circuit(value[0])))
        print(f"parameters {value}: TVD {distance:.4f}", flush=True)
        if distance > TOL:
            failures.append(f"parameters {value}")

    # an expectation value through the backend object
    want = Statevector(qc).expectation_value(SparsePauliOp("ZZZZ")).real
    got = compss_wait_on(
        qd.gate_cutting(qc, ["CZ_4"], "ZZZZ", backend=aer, shots=SHOTS, sync=False)
    )
    print(f"expectation: {float(got):+.4f} (exact {want:+.4f})", flush=True)
    if abs(float(got) - want) > 2 * TOL:
        failures.append("expectation")

    print("RESULT:", "FAIL " + ", ".join(failures) if failures else "PASS", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
