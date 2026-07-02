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

"""Verify ``gate_cutting_sampling`` under distributed (PyCOMPSs) execution.

The sampling reconstruction fans a task out per subcircuit
(``_sample_gate_subcircuit``) and gathers with ``compss_wait_on``. This app
runs it and checks the reconstructed distribution against the *exact*
statevector distribution, so the same script doubles as a correctness smoke
test of the parallel path: task-boundary serialisation of the component DAGs,
order-independent (positional) reconstruction, and -- crucially -- that no
subcircuit helper is itself a ``@task`` (PyCOMPSs has no nested tasks, so a
``@task`` calling a ``@task`` would return an unresolved future on the worker).

Run serially (no PyCOMPSs; the runtime shim runs every task inline)::

    python3 examples/verify_sampling_pycompss.py

Run distributed with PyCOMPSs, generating the task-dependency graph (``-g``)::

    runcompss -g --pythonpath=$(pwd) --python_interpreter=python3 \
        examples/verify_sampling_pycompss.py

``-g`` writes ``complete_graph.dot`` into ``$HOME/.COMPSs/<app>/monitor/``.
Turn it into a PDF with::

    compss_gengraph $HOME/.COMPSs/verify_sampling_pycompss_*/monitor/complete_graph.dot

A correct run shows the sampling tasks fanning out in parallel (one box per
subcircuit per quasiprobability configuration) with **no** task nested inside
another. The app exits non-zero if any backend's reconstruction is off.

Environment knobs:
  QDISLIB_SOFTWARE   comma-separated backends to check (default "qiskit";
                     e.g. "qiskit,qibo,cudaq" -- skipped if not importable).
  QDISLIB_SHOTS      shots per subcircuit (default 8192; sampling has
                     quasiprobability overhead so it wants more than expectation).
  QDISLIB_TVD_TOL    max total-variation distance to still PASS (default 0.06).
"""

import os
import sys

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

import Qdislib as qd

SOFTWARES = [s.strip() for s in os.environ.get("QDISLIB_SOFTWARE", "qiskit").split(",")]
SHOTS = int(os.environ.get("QDISLIB_SHOTS", "8192"))
TVD_TOL = float(os.environ.get("QDISLIB_TVD_TOL", "0.06"))


def build_circuit():
    """A 4-qubit circuit whose min-cost gate cut is *non-terminal*.

    The ``ry`` gates after the ``cz`` chain keep the cut qubits evolving past
    the cut point, which is exactly the case that used to reconstruct wrong --
    so it exercises the collapse-and-continue reconstruction, not just a
    trivial terminal cut.
    """
    qc = QuantumCircuit(4)
    for q in range(4):
        qc.ry(0.6 + 0.1 * q, q)
    qc.cz(0, 1)
    qc.cz(1, 2)
    qc.cz(2, 3)
    qc.ry(0.4, 1)  # post-cut gate -> non-terminal cut on the middle qubits
    qc.ry(0.4, 2)
    return qc


def exact_distribution(qc):
    """Exact output distribution as ``{qubit-0-first bitstring: probability}``."""
    probs = Statevector(qc).probabilities()
    n = qc.num_qubits
    # Statevector index i is little-endian (bit q of i == qubit q). Emit the
    # bitstring qubit-0-first to match gate_cutting_sampling's convention.
    dist = {}
    for i, p in enumerate(probs):
        if p < 1e-12:
            continue
        bitstring = "".join(str((i >> q) & 1) for q in range(n))
        dist[bitstring] = float(p)
    return dist


def total_variation(exact, recon):
    """TVD = 1/2 * sum_x |exact(x) - recon(x)| over the union of supports."""
    keys = set(exact) | set(recon)
    return 0.5 * sum(abs(exact.get(k, 0.0) - recon.get(k, 0.0)) for k in keys)


def main():
    qc = build_circuit()
    gates_cut = qd.find_cut(qc, gate_cut=True, wire_cut=False)
    exact = exact_distribution(qc)

    print(f"Circuit: {qc.num_qubits} qubits, gate cut = {gates_cut}")
    print(f"Backends: {SOFTWARES} | shots/subcircuit: {SHOTS} | TVD tol: {TVD_TOL}")
    print(f"Exact distribution has {len(exact)} non-zero bitstrings.\n")

    all_ok = True
    for software in SOFTWARES:
        # (a) All-in-one gate_cutting_sampling.
        try:
            recon = qd.gate_cutting_sampling(
                qc, gates_cut, shots=SHOTS, software=software
            )
        except ImportError as exc:
            print(f"  {software:7s}  SKIP (backend not importable: {exc})")
            continue
        except Exception as exc:  # noqa: BLE001 -- surface any parallel-path failure
            print(f"  {software:7s}  all-in-one FAIL (raised: {exc})")
            all_ok = False
            continue
        tvd_all = total_variation(exact, recon)

        # (b) Three-step workflow: generate -> sample each -> reconstruct.
        try:
            subs = qd.gate_cutting_sampling_subcircuits(
                qc, gates_cut, software=software
            )
            sampled = [qd.sample_subcircuit(s, shots=SHOTS) for s in subs]
            recon3 = subs.reconstruct(sampled)
            tvd_3step = total_variation(exact, recon3)
        except Exception as exc:  # noqa: BLE001
            print(f"  {software:7s}  3-step FAIL (raised: {exc})")
            all_ok = False
            continue

        ok = tvd_all < TVD_TOL and tvd_3step < TVD_TOL
        all_ok = all_ok and ok
        print(
            f"  {software:7s}  all-in-one TVD {tvd_all:.4f} | "
            f"3-step TVD {tvd_3step:.4f}  ->  {'PASS' if ok else 'FAIL'}"
        )

    print("\nRESULT:", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
