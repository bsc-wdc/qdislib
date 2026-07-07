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

"""Qdislib example application for distributed execution with PyCOMPSs.

Run serially (no PyCOMPSs)::

    python3 examples/compss_app.py

Run distributed with PyCOMPSs::

    runcompss --pythonpath=$(pwd) --python_interpreter=python3 examples/compss_app.py

Under ``runcompss`` the per-subcircuit evaluations execute as COMPSs ``@task``
instances and are synchronised with ``compss_wait_on`` inside Qdislib, so this
app doubles as a smoke test of the parallel execution path (circuit
serialisation across tasks, reconstruction, etc.). It exits non-zero if the
reconstructed values do not match the analytical reference.
"""

import os
import sys

import numpy as np
import qibo
from qibo import models, gates

import Qdislib as qd

TOLERANCE = 0.2

# Backend used to build/execute the subcircuits (override with QDISLIB_SOFTWARE).
SOFTWARE = os.environ.get("QDISLIB_SOFTWARE", "qibo")


def build_circuit():
    """A small connected circuit with several CZ gates (cuttable both ways)."""
    circuit = models.Circuit(5)
    circuit.add(gates.H(0))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.CZ(1, 2))
    circuit.add(gates.RY(2, np.pi / 5))
    circuit.add(gates.CZ(2, 3))
    circuit.add(gates.CZ(3, 4))
    circuit.add(gates.RX(1, np.pi / 4))
    return circuit


def main():
    qibo.set_backend("numpy")
    observable = "Z" * build_circuit().nqubits

    analytical = qd.analytical_solution(build_circuit(), observable)
    print(f"Backend: {SOFTWARE}")
    print(f"Analytical <{observable}> = {analytical:+.4f}")

    gate_cut = qd.find_cut(build_circuit(), gate_cut=True, wire_cut=False)
    gate_value = qd.gate_cutting(build_circuit(), gate_cut, software=SOFTWARE)
    print(
        f"Gate cutting {gate_cut} -> {gate_value:+.4f} "
        f"(|diff| = {abs(gate_value - analytical):.4f})"
    )

    wire_cut = qd.find_cut(build_circuit(), gate_cut=False, wire_cut=True)
    wire_value = qd.wire_cutting(build_circuit(), wire_cut, software=SOFTWARE)
    print(
        f"Wire cutting {wire_cut} -> {wire_value:+.4f} "
        f"(|diff| = {abs(wire_value - analytical):.4f})"
    )

    # Recommended three-step workflow: generate -> evaluate -> reconstruct.
    # Under PyCOMPSs the subcircuits are serializable Subcircuit wrappers, each
    # evaluated as its own task; in serial they are native circuits.
    gate_subs = qd.gate_cutting_subcircuits(
        build_circuit(), gate_cut, software=SOFTWARE
    )
    gate_step = gate_subs.reconstruct([qd.expectation_value(s) for s in gate_subs])
    print(
        f"3-step gate -> {gate_step:+.4f} (|diff| = {abs(gate_step - analytical):.4f})"
    )

    wire_subs = qd.wire_cutting_subcircuits(
        build_circuit(), wire_cut, software=SOFTWARE
    )
    wire_step = wire_subs.reconstruct([qd.expectation_value(s) for s in wire_subs])
    print(
        f"3-step wire -> {wire_step:+.4f} (|diff| = {abs(wire_step - analytical):.4f})"
    )

    # Optimized (fragment-separated, 3**n + 4**n) wire cutting -- the fan-out of
    # fragment @tasks feeds a reconstruction @task, so this exercises that parallel
    # path under runcompss. Falls back to the standard method if not separable.
    wire_opt = qd.wire_cutting(
        build_circuit(), wire_cut, software=SOFTWARE, optimized=True
    )
    print(
        f"Wire cutting (optimized) {wire_cut} -> {wire_opt:+.4f} "
        f"(|diff| = {abs(wire_opt - analytical):.4f})"
    )

    wire_opt_subs = qd.wire_cutting_subcircuits(
        build_circuit(), wire_cut, software=SOFTWARE, optimized=True
    )
    wire_opt_step = wire_opt_subs.reconstruct(
        [qd.expectation_value(s) for s in wire_opt_subs]
    )
    print(
        f"3-step wire (optimized) -> {wire_opt_step:+.4f} "
        f"(|diff| = {abs(wire_opt_step - analytical):.4f})"
    )

    # Hadamard-test gate cut: each (configuration, fragment) is one
    # ``_evaluate_hadamard_fragment`` @task and they fan into a single
    # ``_reduce_hadamard`` @task -- the fan-out/fan-in this app renders under
    # ``runcompss --graph``. Sampled (shot-based) so it exercises the backend tasks.
    had_value = qd.gate_cutting_hadamard(
        build_circuit(), gate_cut, observable, software=SOFTWARE, shots=100000
    )
    print(
        f"Hadamard gate cutting {gate_cut} -> {had_value:+.4f} "
        f"(|diff| = {abs(had_value - analytical):.4f})"
    )

    had_subs = qd.gate_cutting_hadamard_subcircuits(
        build_circuit(), gate_cut, observable, software=SOFTWARE
    )
    had_step = had_subs.reconstruct(
        [qd.expectation_value(s, shots=100000) for s in had_subs]
    )
    print(
        f"3-step Hadamard -> {had_step:+.4f} (|diff| = {abs(had_step - analytical):.4f})"
    )

    ok = all(
        abs(v - analytical) < TOLERANCE
        for v in (
            gate_value,
            wire_value,
            gate_step,
            wire_step,
            wire_opt,
            wire_opt_step,
            had_value,
            had_step,
        )
    )
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
