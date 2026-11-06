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

"""``sync=False`` returns a future everywhere, and nothing on the master opens one.

Runs in a subprocess under the strict fake PyCOMPSs of :mod:`tests._fake_compss`,
where a task result is an opaque future: a reconstruction done on the master
(``sum(results)``), a function that waits although ``sync=False``, or a task-built
converter called outside a task all fail here as they would under PyCOMPSs.
"""

import importlib.util
import os
import subprocess
import sys
import textwrap
import unittest

from tests import BaseTimedTestCase

QISKIT = importlib.util.find_spec("qiskit") is not None

_SCRIPT = textwrap.dedent("""
    import sys, warnings
    warnings.simplefilter("ignore")
    from tests import _fake_compss as fc
    fc.install()
    import numpy as np
    from qiskit import QuantumCircuit
    import Qdislib as qd

    def the_store():
        return store

    def chain(blocks, kinds):
        # blocks of two qubits in a line; kinds[j] joins block j and j + 1
        rng = np.random.default_rng(1)
        n = 2 * blocks
        qc = QuantumCircuit(n)
        names = lambda: f"{qc.data[-1].operation.name.upper()}_{len(qc.data)}"
        last = {}
        def u(*qs):
            for q in qs:
                qc.ry(rng.uniform(0, 3), q)
                last[q] = names()
        u(*range(n))
        for j in range(blocks):
            qc.cz(2 * j, 2 * j + 1)
        u(*range(n))
        cuts = {"gate": [], "hadamard": [], "wire": []}
        for j, kind in enumerate(kinds):
            a, b = 2 * j + 1, 2 * j + 2
            if kind == "wire":
                u(a)
                source = last[a]
                qc.cz(a, b)
                cuts["wire"].append((source, names()))
            else:
                qc.cz(a, b)
                cuts[kind].append(names())
        u(*range(n))
        return qc, cuts

    gate, gc = chain(3, ["gate", "gate"])
    had, hc = chain(3, ["hadamard", "hadamard"])
    wire, wc = chain(3, ["wire", "wire"])
    mixed, mc = chain(4, ["gate", "hadamard", "wire"])
    ham = lambda n: [(0.5, "Z" * n), (-0.3, "X" + "Z" * (n - 1))]
    exact = [
        (qd.gate_cutting, gate, gc["gate"], {}),
        (qd.gate_cutting, gate, gc["gate"], {"observable": ham(6)}),
        (qd.gate_cutting_sampling, gate, gc["gate"], {}),
        (qd.gate_cutting_hadamard, had, hc["hadamard"], {}),
        (qd.gate_cutting_hadamard_sampling, had, hc["hadamard"], {}),
        (qd.wire_cutting, wire, wc["wire"], {}),
        (qd.wire_cutting, wire, wc["wire"], {"optimized": True}),
        (qd.wire_cutting_sampling, wire, wc["wire"], {}),
        (qd.circuit_cutting, mixed, mc, {"observable": ham(8)}),
        (qd.circuit_cutting, mixed, mc, {"optimized": True}),
        (qd.circuit_cutting_sampling, mixed, mc, {}),
        (qd.circuit_cutting, mixed, mc, {"num_samples": 40}),
    ]
    # with shots, the paths that do not go through circuit_cutting
    sampled = [
        (qd.gate_cutting, gate, gc["gate"], {}),
        (qd.gate_cutting, gate, gc["gate"], {"observable": ham(6)}),
        (qd.gate_cutting_hadamard, had, hc["hadamard"], {}),
        (qd.gate_cutting_hadamard, had, hc["hadamard"], {"observable": ham(6)}),
        (qd.gate_cutting_hadamard_sampling, had, hc["hadamard"], {}),
        (qd.wire_cutting, wire, wc["wire"], {}),
        (qd.wire_cutting, wire, wc["wire"], {"observable": ham(6)}),
        (qd.wire_cutting, wire, wc["wire"], {"optimized": True}),
        (qd.wire_cutting_sampling, wire, wc["wire"], {}),
    ]
    cases = [(f, c, k, dict(kw, shots=None, seed=3)) for f, c, k, kw in exact]
    cases += [(f, c, k, dict(kw, shots=100, seed=3)) for f, c, k, kw in sampled]
    failures = []
    for fn, qc, cuts, kw in cases:
        name = f"{fn.__name__}({sorted(kw.items())})"
        try:
            reference = fn(qc, cuts, sync=True, **kw)
            out = fn(qc, cuts, sync=False, **kw)
            if not isinstance(out, fc.Future):
                failures.append(f"{name}: returned {type(out).__name__}, not a future")
                continue
            value = fc.compss_wait_on(out)
            if isinstance(value, dict):
                diff = max(abs(value.get(b, 0) - reference.get(b, 0))
                           for b in set(value) | set(reference))
            else:
                diff = abs(value - reference)
            if diff > 1e-9:
                failures.append(f"{name}: differs from sync=True by {diff}")
        except Exception as e:
            failures.append(f"{name}: {type(e).__name__}: {e}")
    # the engine's cache: a miss stores from inside the task, a hit is read on the
    # master and mixed with the futures of the circuits still to run
    # a store every worker shares, as LMDB or Redis are on a cluster (an in-memory one
    # reaches a task as a copy)
    class SharedStore(qd.InMemoryCircuitCache):
        def __reduce__(self):
            return (the_store, ())
    store = SharedStore()
    for fn, qc, cuts, kw in (
        (qd.circuit_cutting, mixed, mc, {"observable": ham(8)}),
        (qd.gate_cutting, gate, gc["gate"], {"shots": 100, "seed": 3}),
        (qd.wire_cutting_sampling, wire, wc["wire"], {"shots": None}),
    ):
        name = f"cached {fn.__name__}"
        try:
            kw = dict({"shots": None}, **kw)
            first = fc.compss_wait_on(fn(qc, cuts, sync=False, cache=store, **kw))
            second = fn(qc, cuts, sync=False, cache=store, **kw)
            if not isinstance(second, fc.Future):
                failures.append(f"{name}: returned {type(second).__name__}")
                continue
            second = fc.compss_wait_on(second)
            if isinstance(first, dict):
                diff = max(abs(first.get(b, 0) - second.get(b, 0))
                           for b in set(first) | set(second))
            else:
                diff = abs(first - second)
            if diff > 1e-9:
                failures.append(f"{name}: a hit differs by {diff}")
        except Exception as e:
            failures.append(f"{name}: {type(e).__name__}: {e}")
    if not store.stats()["hits"]:
        failures.append("cached runs never hit")
    # error bars and a shot budget through futures
    try:
        value = fc.compss_wait_on(qd.circuit_cutting(
            mixed, mc, ham(8), shots=100, seed=3, sync=False, error=20))
        if not value.error > 0:
            failures.append("error=: no error bar")
        value = fc.compss_wait_on(qd.circuit_cutting(
            gate, gc, "Z" * 6, total_shots=5000, seed=3, sync=False, error=20))
        if value.shots != 5000:
            failures.append(f"total_shots=: spent {value.shots}")
        dist = fc.compss_wait_on(qd.circuit_cutting_sampling(
            gate, gc, shots=100, seed=3, sync=False, error=20))
        if not dist.errors or not max(dist.errors.values()) > 0:
            failures.append("error= on a distribution: no error bars")
    except Exception as e:
        failures.append(f"error/total_shots: {type(e).__name__}: {e}")
    # parameter_values=: one future per set of values, the cut planned once
    try:
        from qiskit.circuit import Parameter

        t = Parameter("t")
        pq = QuantumCircuit(4)
        pq.ry(t, 0); pq.cx(0, 1); pq.cz(1, 2); pq.ry(0.4, 2); pq.cx(2, 3)
        sets = [[0.3], [1.2], [2.5]]
        for per_task in (None, 3):
            futures = qd.circuit_cutting(
                pq, {"gate": ["CZ_3"]}, "ZIII", shots=None, sync=False,
                parameter_values=sets, circuits_per_task=per_task)
            if not isinstance(futures, list) or len(futures) != len(sets):
                failures.append("parameter_values: not one future per set")
            for future, v in zip(futures, sets):
                exact = qd.analytical_solution(pq.assign_parameters(v), "ZIII")
                if abs(fc.compss_wait_on(future) - exact) > 1e-9:
                    failures.append(f"parameter_values({per_task}): wrong value")
    except Exception as e:
        failures.append(f"parameter_values: {type(e).__name__}: {e}")
    # circuits_per_task=: the same values from fewer tasks, every cut type, and a
    # last batch of a single circuit (PyCOMPSs returns it as one future, not a list)
    for fn, qc, cuts, kw in (
        (qd.circuit_cutting, mixed, mc, {"observable": ham(8), "shots": None}),
        (qd.circuit_cutting_sampling, mixed, mc, {"shots": 100, "seed": 3}),
        (qd.gate_cutting, gate, gc["gate"], {"shots": 100, "seed": 3}),
        (qd.gate_cutting_hadamard_sampling, had, hc["hadamard"], {"shots": None}),
        (qd.wire_cutting, wire, wc["wire"], {"shots": 100, "seed": 3, "error": 20}),
    ):
        name = f"circuits_per_task {fn.__name__}"
        try:
            before = fc.calls["tasks"]
            reference = fc.compss_wait_on(fn(qc, cuts, sync=False, **kw))
            plain = fc.calls["tasks"] - before
            for per_task in (1, 3, 1000):
                before = fc.calls["tasks"]
                value = fc.compss_wait_on(
                    fn(qc, cuts, sync=False, circuits_per_task=per_task, **kw))
                used = fc.calls["tasks"] - before
                if isinstance(value, dict):
                    diff = max(abs(value.get(b, 0) - reference.get(b, 0))
                               for b in set(value) | set(reference))
                else:
                    diff = abs(value - reference)
                    if getattr(reference, "error", 0) and abs(value.error - reference.error) > 1e-12:
                        failures.append(f"{name}({per_task}): other error bar")
                if diff > 1e-9:
                    failures.append(f"{name}({per_task}): differs by {diff}")
                if per_task > 1 and not used < plain:
                    failures.append(f"{name}({per_task}): {used} tasks, unbatched {plain}")
        except Exception as e:
            failures.append(f"{name}: {type(e).__name__}: {e}")
    for bad in (0, True, 2.5):
        try:
            qd.circuit_cutting(mixed, mc, shots=None, circuits_per_task=bad)
            failures.append(f"circuits_per_task={bad!r} accepted")
        except ValueError:
            pass
    for make, qc, cuts in (
        (qd.gate_cutting_subcircuits, gate, gc["gate"]),
        (qd.gate_cutting_hadamard_subcircuits, had, hc["hadamard"]),
        (qd.wire_cutting_subcircuits, wire, wc["wire"]),
        (qd.circuit_cutting_subcircuits, mixed, mc),
    ):
        try:  # the run step's futures straight into reconstruct
            subs = make(qc, cuts)
            expected = qd.analytical_solution(qc, "Z" * qc.num_qubits)
            value = subs.reconstruct([qd.expectation_value(s, shots=None) for s in subs])
            if abs(fc.compss_wait_on(value) - expected) > 1e-9:
                failures.append(f"{make.__name__}: wrong value")
        except Exception as e:
            failures.append(f"{make.__name__}: {type(e).__name__}: {e}")
    try:  # counts from your own runs: one conversion task per subcircuit
        from qiskit_aer import AerSimulator
        subs = qd.circuit_cutting_subcircuits(gate, gc)
        circuits = []
        for c in subs.circuits("qiskit", defer_measurements=True):
            c = c.copy(); c.measure_all(); circuits.append(c)
        result = AerSimulator().run(circuits, shots=20000, seed_simulator=1).result()
        before = fc.calls["tasks"]
        value = subs.reconstruct_counts(
            [result.get_counts(i) for i in range(len(circuits))], bit_order="q0_last")
        if fc.calls["tasks"] - before < len(subs):
            failures.append("reconstruct_counts: conversions not run as tasks")
        value = fc.compss_wait_on(value)
        if abs(value - qd.analytical_solution(gate, "Z" * gate.num_qubits)) > 0.05:
            failures.append(f"reconstruct_counts: wrong value {value}")
    except Exception as e:
        failures.append(f"reconstruct_counts: {type(e).__name__}: {e}")
    print("\\n".join(failures) or "ALL OK")
    """)


_WAITS_SCRIPT = textwrap.dedent("""
    import warnings
    warnings.simplefilter("ignore")
    from tests import _fake_compss as fc
    fc.install()
    from qiskit import QuantumCircuit
    import Qdislib as qd

    qc = QuantumCircuit(2)
    qc.x(0)
    # no cut, on the 1.x path (a deprecated option): a value, not a future
    value = qd.wire_cutting(qc, [], "ZI", gpu=True, shots=10, seed=1)
    future = qd.wire_cutting(qc, [], "ZI", gpu=True, shots=10, seed=1, sync=False)
    print(isinstance(value, float), isinstance(future, fc.Future), value)
    """)


@unittest.skipUnless(QISKIT, "qiskit not installed")
class SyncFalseTest(BaseTimedTestCase):
    def test_sync_true_waits_without_cuts_on_the_deprecated_path(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(
            os.environ, PYTHONPATH=root + os.pathsep + os.environ.get("PYTHONPATH", "")
        )
        run = subprocess.run(
            [sys.executable, "-c", _WAITS_SCRIPT],
            cwd=root,
            env=env,
            capture_output=True,
            text=True,
            timeout=300,
        )
        self.assertEqual(run.returncode, 0, run.stderr[-3000:])
        self.assertEqual(run.stdout.strip().splitlines()[-1], "True True -1.0")

    def test_every_function_returns_a_future(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = dict(
            os.environ, PYTHONPATH=root + os.pathsep + os.environ.get("PYTHONPATH", "")
        )
        run = subprocess.run(
            [sys.executable, "-c", _SCRIPT],
            cwd=root,
            env=env,
            capture_output=True,
            text=True,
            timeout=900,
        )
        self.assertEqual(run.returncode, 0, run.stderr[-3000:])
        self.assertEqual(run.stdout.strip().splitlines()[-1], "ALL OK", run.stdout)


if __name__ == "__main__":
    unittest.main()
