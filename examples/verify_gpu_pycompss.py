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

"""GPU tasks under real PyCOMPSs, on a GPU node: a backend on a GPU gets one per task,
a CPU backend none, and CUDA-Q on the GPU gives the exact distribution. Needs
``QDISLIB_PROBE_FILE`` (a path the workers can write) and ``examples/`` on the
PyCOMPSs path. Prints ``RESULT: PASS`` or ``RESULT: FAIL``."""

import collections
import json
import os
import sys

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

import Qdislib as qd
from Qdislib.utils.runtime import compss_wait_on

from _gpu_probe import ProbeCpu, ProbeCudaq, ProbeGpu


def main():
    qc = QuantumCircuit(4)
    qc.ry(0.4, 0)
    qc.rx(1.1, 2)
    qc.cx(0, 1)
    qc.cz(1, 2)
    qc.ry(0.7, 2)
    qc.cx(2, 3)
    want = {k[::-1]: v for k, v in Statevector(qc).probabilities_dict().items()}
    cut = {"hadamard": ["CZ_4"]}
    config = qd.BackendConfig(num_qubits=4)
    runs = [
        (ProbeGpu(config), {}),
        (ProbeGpu(config), {"circuits_per_task": 4}),
        (ProbeCpu(config), {}),
        (ProbeCudaq(qd.BackendConfig(num_qubits=4, extra={"target": "nvidia"})), {}),
    ]
    failures = []
    for backend, options in runs:
        got = compss_wait_on(
            qd.circuit_cutting_sampling(
                qc, cut, backend=backend, shots=20000, seed=1, **options
            )
        )
        distance = 0.5 * sum(
            abs(got.get(k, 0.0) - want.get(k, 0.0)) for k in set(got) | set(want)
        )
        print(type(backend).__name__, options, f"TVD {distance:.4f}", flush=True)
        if distance > 0.05:
            failures.append(f"{type(backend).__name__} TVD {distance:.3f}")
    with open(os.environ["QDISLIB_PROBE_FILE"]) as f:
        rows = [json.loads(line) for line in f if line.strip()]
    by_kind = collections.defaultdict(list)
    for row in rows:
        by_kind[row["kind"]].append(row)
    for kind, group in by_kind.items():
        print(
            kind,
            len(group),
            "tasks; GPUs:",
            dict(collections.Counter(str(r["gpus"]) for r in group)),
        )
    for kind in ("gpu", "cudaq"):
        bad = [
            r["gpus"] for r in by_kind[kind] if not (r["gpus"] and "," not in r["gpus"])
        ]
        if bad or not by_kind[kind]:
            failures.append(f"{kind}: not exactly one GPU per task")
    if any(r["gpus"] for r in by_kind["cpu"]):
        failures.append("cpu: COMPSs gave a CPU task a GPU")
    if any(r["target"] != "nvidia" for r in by_kind["cudaq"]):
        failures.append("cudaq: not on the nvidia target")
    print("RESULT:", "FAIL " + "; ".join(failures) if failures else "PASS", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
