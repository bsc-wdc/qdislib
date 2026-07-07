#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#
# -*- coding: utf-8 -*-

"""Hadamard-test gate cut under PyCOMPSs -- a focused task-graph showcase.

Runs the ancilla / double-Hadamard-test cut (arXiv:2403.01018) three ways, each of
which fans one ``_evaluate_hadamard_fragment`` @task per (configuration, fragment) into
a single ``_reduce_hadamard`` fan-in @task::

    runcompss -g --pythonpath=$(pwd) --python_interpreter=python3 \\
        examples/hadamard_cut_compss.py

The ``-g`` flag writes the task dependency graph to
``$HOME/.COMPSs/<app>/monitor/complete_graph.dot``. A ``shots`` value (not the exact
statevector path) is used so the fragments really execute as COMPSs tasks. Exits
non-zero unless every variant matches the analytical reference.
"""

import os
import sys

import qibo
from qibo import gates, models

import Qdislib as qd
from Qdislib.utils.runtime import compss_wait_on

TOLERANCE = 0.15
SOFTWARE = os.environ.get("QDISLIB_SOFTWARE", "qibo")
SHOTS = int(os.environ.get("QDISLIB_SHOTS", "8192"))
NUM_SAMPLES = int(os.environ.get("QDISLIB_NUM_SAMPLES", "8"))


def build_circuit():
    """A small circuit with a single cuttable CZ (clean 12-task fan-out for ZZ)."""
    c = models.Circuit(2)
    c.add(gates.H(0))
    c.add(gates.CZ(0, 1))
    c.add(gates.RY(0, 0.8))
    return c


def main():
    qibo.set_backend("numpy")
    circuit, cut, observable = build_circuit(), ["CZ_2"], "ZZ"
    analytical = qd.analytical_solution(build_circuit(), observable)
    print(f"Backend: {SOFTWARE}  shots={SHOTS}  num_samples={NUM_SAMPLES}")
    print(f"Analytical <{observable}> = {analytical:+.4f}")

    # (1) + (2) Submit the enumeration and the Monte-Carlo sampling fan-outs WITHOUT
    # blocking (sync=False), so both clusters' fragment tasks are in flight at once --
    # no barrier between them. The master only synchronises once, at the end.
    enum_future = qd.gate_cutting_hadamard(
        circuit, cut, observable, software=SOFTWARE, shots=SHOTS, sync=False
    )
    samp_future = qd.gate_cutting_hadamard_sampling(
        build_circuit(),
        cut,
        observable,
        software=SOFTWARE,
        shots=SHOTS,
        num_samples=NUM_SAMPLES,
        seed=1,
        sync=False,
    )
    enum = compss_wait_on(enum_future)
    samp = compss_wait_on(samp_future)
    print(
        f"Hadamard (enumerate)  -> {enum:+.4f} (|diff| = {abs(enum - analytical):.4f})"
    )
    print(
        f"Hadamard (sampling)   -> {samp:+.4f} (|diff| = {abs(samp - analytical):.4f})"
    )

    # (3) Three-step: build subcircuits, evaluate each as a task, reconstruct.
    subs = qd.gate_cutting_hadamard_subcircuits(
        build_circuit(), cut, observable, software=SOFTWARE
    )
    step = subs.reconstruct([qd.expectation_value(s, shots=SHOTS) for s in subs])
    print(
        f"Hadamard (3-step)     -> {step:+.4f} (|diff| = {abs(step - analytical):.4f})"
    )

    # Enumeration and the three-step evaluate every configuration, so they must be
    # accurate. The Monte-Carlo estimator is deliberately run with few samples here so
    # the task graph stays readable, which leaves high variance (~phi/sqrt(N)); the
    # unit tests validate its convergence, so here we only guard against a gross error.
    ok = all(abs(v - analytical) < TOLERANCE for v in (enum, step))
    ok = ok and abs(samp - analytical) < 1.0  # sign / order-of-magnitude sanity only
    print("RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
