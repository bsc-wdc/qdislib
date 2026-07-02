#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#

# -*- coding: utf-8 -*-

"""Generate the Qdislib learning-path notebooks.

The notebooks under ``notebooks/`` are generated from this single script so the
whole learning path stays consistent, valid and easy to regenerate:

    python scripts/build_learning_path.py

Each notebook is standalone (its own imports + example circuit) but they are
numbered to be read in order. Cells that need a runtime we can't assume in the
docs build (PyCOMPSs, CUDA-Q, IBM Quantum, the BSC QPU bridge) are tagged
``illustrative`` -- ``nbsphinx_execute = "never"`` means none are executed at
build time anyway.
"""

import os

import nbformat as nbf
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

NB_DIR = os.path.join(os.path.dirname(__file__), "..", "notebooks")

# A compact example circuit reused (defined inline) across the notebooks: the CZ
# chain makes it cuttable, and the RY rotations keep <Z...Z> clearly non-zero (an
# H on qubit 0 would drive <Z> to 0 and make the printed values vacuous).
EXAMPLE = '''\
from qibo import models, gates


def example_circuit():
    """A small 5-qubit circuit: cuttable (CZ chain), non-trivial <ZZZZZ> ~= 0.44."""
    circuit = models.Circuit(5)
    circuit.add(gates.RY(0, 0.8))
    circuit.add(gates.CZ(0, 1))
    circuit.add(gates.RY(1, 0.5))
    circuit.add(gates.CZ(1, 2))
    circuit.add(gates.RY(2, 0.6))
    circuit.add(gates.CZ(2, 3))
    circuit.add(gates.RY(3, 0.4))
    circuit.add(gates.CZ(3, 4))
    circuit.add(gates.RX(4, 0.3))
    return circuit


circuit = example_circuit()
print(circuit.draw())'''


def md(text):
    return new_markdown_cell(text)


def code(text):
    return new_code_cell(text)


def illus(text):
    """A code cell that documents a runtime we don't execute here (tagged)."""
    cell = new_code_cell(text)
    cell.metadata["tags"] = ["illustrative"]
    return cell


def img_md(rel_path, alt):
    """Markdown cell embedding a PNG (base64) so the notebook is self-contained.

    Lets a notebook show a rendered result (e.g. the expected reconstructed image)
    without executing anything or shipping a separate file.
    """
    import base64

    full = os.path.join(NB_DIR, rel_path)
    with open(full, "rb") as handle:
        data = base64.b64encode(handle.read()).decode()
    return new_markdown_cell(f"![{alt}](data:image/png;base64,{data})")


NOTEBOOKS = {}

# --------------------------------------------------------------------------- #
# 01 — Quickstart
# --------------------------------------------------------------------------- #
NOTEBOOKS["01_quickstart"] = [
    md(
        "# Qdislib — 01 · Quickstart\n"
        "\n"
        "**Qdislib** runs quantum circuits that are too large for one device by\n"
        "**cutting** them into smaller subcircuits, executing those independently\n"
        "(optionally distributed with PyCOMPSs, or on real QPUs), and\n"
        "**reconstructing** the original expectation value from the pieces.\n"
        "\n"
        "This is a learning path — read the notebooks in order:\n"
        "\n"
        "1. **Quickstart** (this notebook) — your first cut.\n"
        "2. **The cutting workflow** — generate → evaluate → reconstruct.\n"
        "3. **Backends** — Qibo / Qiskit / CUDA-Q, and QASM interop.\n"
        "4. **Observables** — Pauli strings and weighted Hamiltonians.\n"
        "5. **Caching & performance** — the semantic cache and logging.\n"
        "6. **Distributed & QPU** — PyCOMPSs and real quantum hardware.\n"
        "\n"
        "Install: `pip install Qdislib[qibo]` (a Qibo-only install is enough to\n"
        "follow along; add `qiskit`, `cudaq`, `cache` as needed)."
    ),
    md(
        "## 1. Build a circuit\n"
        "\n"
        "We use Qibo here, but Qdislib is backend-agnostic (notebook 3)."
    ),
    code(EXAMPLE),
    code("import Qdislib as qd"),
    md(
        "## 2. The one-liner: `estimate`\n"
        "\n"
        "`qd.estimate` does the whole pipeline — find a cut, run the subcircuits,\n"
        "reconstruct — in a single call. Here we estimate `<ZZZZZ>`, asking for\n"
        "subcircuits of at most 3 qubits."
    ),
    code(
        'value = qd.estimate(example_circuit(), "ZZZZZ", max_qubits=3, shots=8192)\n'
        'print("estimate  :", value)'
    ),
    md(
        "For a small circuit we can check against the exact (statevector) value —\n"
        "the two agree up to shot noise."
    ),
    code(
        'exact = qd.analytical_solution(example_circuit(), "ZZZZZ")\n'
        'print("analytical:", exact)'
    ),
    md(
        "## 3. Under the hood: find a cut, then cut\n"
        "\n"
        "`estimate` is a convenience wrapper. The two explicit steps are\n"
        "`find_cut` (pick where to cut) and `gate_cutting` / `wire_cutting` (cut,\n"
        "run and reconstruct). **Gate cutting** splits a two-qubit gate:"
    ),
    code(
        "circuit = example_circuit()\n"
        "gate_cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        'print("gate cut:", gate_cut)\n'
        "\n"
        "value = qd.gate_cutting(example_circuit(), gate_cut, shots=8192)\n"
        'print("reconstructed <ZZZZZ>:", value)'
    ),
    md(
        "## 4. Wire cutting\n"
        "\n"
        "**Wire cutting** severs a qubit wire between two gates instead of a gate.\n"
        "Same call shape:"
    ),
    code(
        "circuit = example_circuit()\n"
        "wire_cut = qd.find_cut(circuit, gate_cut=False, wire_cut=True)\n"
        'print("wire cut:", wire_cut)\n'
        "\n"
        "value = qd.wire_cutting(example_circuit(), wire_cut, shots=8192)\n"
        'print("reconstructed <ZZZZZ>:", value)'
    ),
    md(
        "---\n"
        "Next: **02 · The cutting workflow** — the recommended three-step API that\n"
        "lets you inspect and distribute the subcircuits."
    ),
]

# --------------------------------------------------------------------------- #
# 02 — The cutting workflow
# --------------------------------------------------------------------------- #
NOTEBOOKS["02_cutting_workflow"] = [
    md(
        "# Qdislib — 02 · The cutting workflow\n"
        "\n"
        "`gate_cutting` / `wire_cutting` do everything in one call. The\n"
        "**recommended workflow** breaks that into three steps you control:\n"
        "\n"
        "1. **generate** the subcircuits,\n"
        "2. **evaluate** each one (locally, distributed, or on a QPU),\n"
        "3. **reconstruct** the original expectation value.\n"
        "\n"
        "This is what you want for distributed runs and for inspecting the pieces."
    ),
    code(EXAMPLE),
    code("import Qdislib as qd"),
    md("## 1. Generate the subcircuits\n\n`find_cut` picks the cut; `gate_cutting_subcircuits` produces the pieces."),
    code(
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        'print("cut:", cut)\n'
        "\n"
        "subcircuits = qd.gate_cutting_subcircuits(example_circuit(), cut, software=\"qibo\")\n"
        'print("number of subcircuits:", len(subcircuits))'
    ),
    md(
        "## 2. Evaluate each subcircuit\n"
        "\n"
        "`expectation_value` runs one subcircuit and returns its expectation. This\n"
        "is the loop you'd replace to dispatch to a cluster or a QPU."
    ),
    code("values = [qd.expectation_value(s, shots=8192) for s in subcircuits]"),
    md(
        "## 3. Reconstruct\n"
        "\n"
        "The object returned by `*_subcircuits` remembers the cut bookkeeping, so\n"
        "`reconstruct` needs nothing else."
    ),
    code(
        "result = subcircuits.reconstruct(values)\n"
        'print("reconstructed:", result)\n'
        'print("analytical   :", qd.analytical_solution(example_circuit(), "ZZZZZ"))'
    ),
    md("## Visualise the cut\n\n`draw_cut` shows the circuit DAG with the cut highlighted."),
    illus("qd.draw_cut(example_circuit(), cut)  # needs matplotlib; opens a figure"),
    md(
        "## Wire cutting — same three steps\n"
        "\n"
        "Only the generator function changes."
    ),
    code(
        "circuit = example_circuit()\n"
        "wire_cut = qd.find_cut(circuit, gate_cut=False, wire_cut=True)\n"
        "subcircuits = qd.wire_cutting_subcircuits(example_circuit(), wire_cut, software=\"qibo\")\n"
        "values = [qd.expectation_value(s, shots=8192) for s in subcircuits]\n"
        'print("wire reconstruction:", subcircuits.reconstruct(values))'
    ),
    md(
        "**Gate vs wire.** Gate cutting splits a two-qubit gate (cost grows fast\n"
        "with the number of cuts); wire cutting severs a qubit wire. `find_cut`\n"
        "can weigh both — see notebook 3 onwards and the docs for guidance.\n"
        "\n"
        "---\n"
        "Next: **03 · Backends** — the same workflow on Qibo, Qiskit and CUDA-Q."
    ),
]

# --------------------------------------------------------------------------- #
# 03 — Backends
# --------------------------------------------------------------------------- #
NOTEBOOKS["03_backends"] = [
    md(
        "# Qdislib — 03 · Backends and interop\n"
        "\n"
        "Qdislib is **backend-agnostic**: the same workflow drives Qibo, Qiskit and\n"
        "CUDA-Q — only the `software=` argument changes. Any framework that exports\n"
        "OpenQASM 2.0 can also feed in via `from_qasm`."
    ),
    code(EXAMPLE),
    code(
        "import Qdislib as qd\n"
        "\n"
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        'print("cut:", cut)'
    ),
    md("## Qibo"),
    code('print("qibo  :", qd.gate_cutting(example_circuit(), cut, software="qibo", shots=8192))'),
    md("## Qiskit\n\nSame circuit and cut, executed on Qiskit Aer (`pip install Qdislib[qiskit]`)."),
    code('print("qiskit:", qd.gate_cutting(example_circuit(), cut, software="qiskit", shots=8192))'),
    md(
        "## CUDA-Q\n"
        "\n"
        "With `pip install Qdislib[cudaq]` you can build the circuit as a CUDA-Q\n"
        "kernel and run with `software=\"cudaq\"` (GPU when available, CPU otherwise).\n"
        "`find_cut` works directly on a kernel too."
    ),
    illus(
        "import cudaq\n"
        "import Qdislib as qd\n"
        "\n"
        "@cudaq.kernel\n"
        "def kernel():\n"
        "    q = cudaq.qvector(5)\n"
        "    ry(0.8, q[0])\n"
        "    z.ctrl(q[0], q[1])\n"
        "    ry(0.5, q[1])\n"
        "    z.ctrl(q[1], q[2])\n"
        "    ry(0.6, q[2])\n"
        "    z.ctrl(q[2], q[3])\n"
        "    ry(0.4, q[3])\n"
        "    z.ctrl(q[3], q[4])\n"
        "    rx(0.3, q[4])\n"
        "\n"
        "cut = qd.find_cut(kernel, gate_cut=True, wire_cut=False)\n"
        'print(qd.gate_cutting(kernel, cut, software="cudaq", shots=8192))'
    ),
    md(
        "## Interop: `from_qasm`\n"
        "\n"
        "Bring circuits from Cirq, PennyLane, Braket, ... through their OpenQASM 2.0\n"
        "export. Here we round-trip a Qibo circuit through QASM into a Qiskit one\n"
        "and cut that."
    ),
    code(
        "qasm = example_circuit().to_qasm()\n"
        'qc = qd.from_qasm(qasm, software="qiskit")\n'
        "cut = qd.find_cut(qc, gate_cut=True, wire_cut=False)\n"
        'print("cut on the imported circuit:", cut)\n'
        'print("value:", qd.gate_cutting(qd.from_qasm(qasm, software="qiskit"), cut, software="qiskit", shots=8192))'
    ),
    md(
        "---\n"
        "Next: **04 · Observables** — beyond `<ZZZZZ>`."
    ),
]

# --------------------------------------------------------------------------- #
# 04 — Observables
# --------------------------------------------------------------------------- #
NOTEBOOKS["04_observables"] = [
    md(
        "# Qdislib — 04 · Observables\n"
        "\n"
        "So far we measured `<ZZZZZ>`. Qdislib accepts any per-qubit Pauli string\n"
        "over `I, X, Y, Z`, and **weighted sums** of them (a Hamiltonian)."
    ),
    code(EXAMPLE),
    code(
        "import Qdislib as qd\n"
        "\n"
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)"
    ),
    md(
        "## Pauli strings\n"
        "\n"
        "`Z` measures in the computational basis; `I` ignores a qubit (including one\n"
        "on a cut — the reconstruction handles it); `X`/`Y` apply the basis rotation\n"
        "for you. Compare each to the exact value."
    ),
    code(
        'for obs in ["ZZZZZ", "ZIZIZ", "XZIZI", "XZZZY"]:\n'
        "    value = qd.gate_cutting(example_circuit(), cut, observables=obs, shots=20000)\n"
        "    exact = qd.analytical_solution(example_circuit(), obs)\n"
        '    print(f"{obs}:  cut={value:+.3f}   exact={exact:+.3f}")'
    ),
    md(
        "> **Sampling variance.** Each cut multiplies the estimator's variance (the\n"
        "> quasiprobability *sampling overhead*), and `X`/`Y` observables add more\n"
        "> on top, so measurements involving them need **more shots** to converge —\n"
        "> the estimate stays unbiased, it just gets noisier. Bump `shots` (or make\n"
        "> fewer cuts) if a value looks off."
    ),
    md(
        "## Weighted sums (Hamiltonians)\n"
        "\n"
        "Pass a list of `(coefficient, pauli)` pairs or a `{pauli: coefficient}`\n"
        "mapping; Qdislib evaluates `<Σ cᵢ Pᵢ> = Σ cᵢ <Pᵢ>` term by term."
    ),
    code(
        "hamiltonian = [(0.5, \"ZZIII\"), (-1.0, \"IIZZI\"), (0.25, \"ZIZII\")]\n"
        "value = qd.gate_cutting(example_circuit(), cut, observables=hamiltonian, shots=20000)\n"
        "exact = qd.analytical_solution(example_circuit(), hamiltonian)\n"
        'print(f"<H>:  cut={value:+.3f}   exact={exact:+.3f}")'
    ),
    md(
        "The same `observables=` argument works for `wire_cutting`, `estimate` and\n"
        "the three-step `*_subcircuits` workflow.\n"
        "\n"
        "---\n"
        "Next: **05 · Caching & performance**."
    ),
]

# --------------------------------------------------------------------------- #
# 05 — Caching & performance
# --------------------------------------------------------------------------- #
NOTEBOOKS["05_caching"] = [
    md(
        "# Qdislib — 05 · Caching & performance\n"
        "\n"
        "Cutting produces many *functionally equivalent* subcircuits (repeated\n"
        "fragments, or the same ansatz across iterations). The **semantic cache**\n"
        "keys each subcircuit by a ZX-reduced Weisfeiler–Lehman hash and simulates\n"
        "an equivalent fragment only once. Install with `pip install Qdislib[cache]`."
    ),
    code(EXAMPLE),
    code("import Qdislib as qd"),
    md(
        "## Enable a cache\n"
        "\n"
        "You choose the store: `\"memory\"` (here), `\"lmdb\"` (persistent on disk),\n"
        "or `\"redis\"` (which becomes the BSC COMPSs/PSCO store under PyCOMPSs).\n"
        "Pass it to any cutting call."
    ),
    code(
        'cache = qd.default_circuit_cache("memory")\n'
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        "\n"
        "value = qd.gate_cutting(example_circuit(), cut, cache=cache, shots=8192)\n"
        'print("value:", value)\n'
        'print("cache:", cache.stats())'
    ),
    md(
        "`stats()` reports hits, misses, hit-rate and size. Caching never changes\n"
        "the result — a circuit that can't be keyed just falls back to a normal\n"
        "simulation."
    ),
    md(
        "## Cache granularity\n"
        "\n"
        "`cache_mode=\"expectation\"` (default) caches the per-subcircuit value keyed\n"
        "by hash + observable. `cache_mode=\"counts\"` caches raw measurement counts\n"
        "keyed by the hash alone, so one simulation serves several observables."
    ),
    code(
        "cache = qd.default_circuit_cache(\"memory\")\n"
        "for obs in [\"ZZZZZ\", \"ZIZIZ\"]:\n"
        "    qd.gate_cutting(example_circuit(), cut, observables=obs,\n"
        "                    cache=cache, cache_mode=\"counts\", shots=8192)\n"
        'print("counts-mode cache:", cache.stats())'
    ),
    md(
        "## Diagnostics\n"
        "\n"
        "There are no `verbose=` flags; turn on Qdislib's diagnostic logging with\n"
        "one switch (or the `QDISLIB_LOG_LEVEL` env var)."
    ),
    code(
        'qd.set_log_level("INFO")   # "DEBUG" for candidate cut scores, backend, ...\n'
        "qd.find_cut(example_circuit(), gate_cut=True, wire_cut=False)\n"
        'qd.set_log_level("WARNING")  # back to quiet'
    ),
    md(
        "Under PyCOMPSs the cache is shared across tasks (Redis/PSCO store) and a\n"
        "placeholder/refresh scheme keeps hit rates up under high concurrency.\n"
        "\n"
        "---\n"
        "Next: **06 · Distributed & QPU**."
    ),
]

# --------------------------------------------------------------------------- #
# 06 — Distributed & QPU (illustrative)
# --------------------------------------------------------------------------- #
NOTEBOOKS["06_distributed_and_qpu"] = [
    md(
        "# Qdislib — 06 · Distributed execution & real QPUs\n"
        "\n"
        "The whole point of cutting is to run the pieces **in parallel** — across a\n"
        "cluster with PyCOMPSs, and/or on real quantum hardware. The cells below are\n"
        "**illustrative**: they need a PyCOMPSs runtime, IBM credentials, or the BSC\n"
        "QPU environment, so they aren't executed in the docs build."
    ),
    md(
        "## 1. Distributed with PyCOMPSs\n"
        "\n"
        "The API is identical — under PyCOMPSs each `expectation_value` runs as a\n"
        "COMPSs task and `*_subcircuits` returns serializable `Subcircuit` wrappers\n"
        "automatically. In a notebook you start the runtime interactively:"
    ),
    illus(
        "import pycompss.interactive as ipycompss\n"
        "ipycompss.start(graph=True, monitor=1000)\n"
        "\n"
        "import Qdislib as qd\n"
        "from qibo import models, gates\n"
        "\n"
        "circuit = models.Circuit(5)\n"
        "circuit.add(gates.RY(0, 0.8)); circuit.add(gates.CZ(0, 1)); circuit.add(gates.RY(1, 0.5))\n"
        "circuit.add(gates.CZ(1, 2)); circuit.add(gates.RY(2, 0.6)); circuit.add(gates.CZ(2, 3))\n"
        "circuit.add(gates.RY(3, 0.4)); circuit.add(gates.CZ(3, 4)); circuit.add(gates.RX(4, 0.3))\n"
        "\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        "subcircuits = qd.gate_cutting_subcircuits(circuit, cut, software=\"qibo\")\n"
        "values = [qd.expectation_value(s) for s in subcircuits]   # scheduled as tasks\n"
        "print(subcircuits.reconstruct(values))\n"
        "\n"
        "ipycompss.stop(sync=True)"
    ),
    md(
        "For batch jobs, run a plain script under `runcompss` instead — see\n"
        "[`examples/compss_app.py`](https://gitlab.bsc.es/wdc/quantum/qdislib) and\n"
        "`examples/README.md`. Nothing in your Qdislib code changes."
    ),
    md(
        "## 2. IBM Quantum\n"
        "\n"
        "With `Qdislib[qiskit]` and an IBM Quantum account (get a token at\n"
        "[quantum.ibm.com](https://quantum.ibm.com/)), set your credentials as\n"
        "environment variables and run the cut subcircuits on real IBM hardware\n"
        "with `qpu=True`."
    ),
    illus(
        "import os\n"
        "\n"
        "# Required:\n"
        'os.environ["IBM_QUANTUM_TOKEN"] = "<your_token>"\n'
        "# Optional:\n"
        'os.environ["IBM_QUANTUM_INSTANCE"] = "<your_instance>"\n'
        'os.environ["IBM_QUANTUM_CHANNEL"] = "ibm_cloud"        # or "ibm_quantum"\n'
        'os.environ["IBM_QUANTUM_MAX_TIME"] = "60"              # batch session, seconds\n'
        'os.environ["IBM_QUANTUM_QPU_NAME"] = "ibm_marrakesh"   # omit -> least-busy backend'
    ),
    md(
        "If `IBM_QUANTUM_QPU_NAME` is unset, Qdislib runs on the least-busy real\n"
        "backend (`service.least_busy(operational=True, simulator=False)`)."
    ),
    illus(
        "import Qdislib as qd\n"
        "from qiskit import QuantumCircuit\n"
        "\n"
        "circuit = QuantumCircuit(10)\n"
        "for i in range(9):\n"
        "    circuit.cz(i, i + 1)\n"
        "\n"
        "# find_cut's IBM path uses Qiskit's hardware-aware cut finder:\n"
        'cut = qd.find_cut(circuit, max_qubits=5, implementation="ibm-ckt")\n'
        "\n"
        "# qpu_dict maps a QPU to its qubit count (the chip's maximum).\n"
        "value = qd.gate_cutting(circuit, cut, software=\"qiskit\",\n"
        "                        qpu=True, qpu_dict={\"IBM_Quantum\": 127})\n"
        "print(value)"
    ),
    md(
        "## 3. BSC QPUs (cross-cluster bridge)\n"
        "\n"
        "At BSC the classical (PyCOMPSs) job and the quantum-control job run on\n"
        "**separate clusters**. They communicate through a shared filesystem\n"
        "directory: a **listener** on the quantum cluster executes submitted\n"
        "subcircuits on a QPU (Blue or Red), while the cut job submits them from the\n"
        "classical side. Both agree on `$QDISLIB_SUBCIRCUITS_DIR`.\n"
        "\n"
        "On the **quantum cluster** (see `examples/qpu_listener.py`):"
    ),
    illus(
        "import Qdislib as qd\n"
        "\n"
        "# $QDISLIB_SUBCIRCUITS_DIR and $RUNCARD_RED must be set.\n"
        'qd.run_qpu_listener(qpu="red", nshots=1000)   # runs until a STOP file / Ctrl-C'
    ),
    md("On the **classical cluster**, point a cut at the QPU (same shared directory):"),
    illus(
        "import Qdislib as qd\n"
        "from qibo import models, gates\n"
        "\n"
        "circuit = models.Circuit(5)\n"
        "circuit.add(gates.RY(0, 0.8)); circuit.add(gates.CZ(0, 1)); circuit.add(gates.RY(1, 0.5))\n"
        "circuit.add(gates.CZ(1, 2)); circuit.add(gates.RY(2, 0.6)); circuit.add(gates.CZ(2, 3))\n"
        "circuit.add(gates.RY(3, 0.4)); circuit.add(gates.CZ(3, 4)); circuit.add(gates.RX(4, 0.3))\n"
        "\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        "value = qd.gate_cutting(circuit, cut, software=\"qibo\",\n"
        "                        qpu=True, qpu_dict={\"MN_Ona\": 5})\n"
        "print(value)"
    ),
    md(
        "## You've finished the learning path 🎉\n"
        "\n"
        "You can now cut a circuit, run it (serially, distributed, or on a QPU),\n"
        "choose a backend, measure arbitrary observables, and cache results. See\n"
        "the API reference for the full surface."
    ),
]


# --------------------------------------------------------------------------- #
# 07 — Hands-on: reconstruct a quantum image by cutting it
# --------------------------------------------------------------------------- #
NOTEBOOKS["07_image_reconstruction"] = [
    md(
        "# Qdislib — 07 · Cutting a quantum image (hands-on)\n"
        "\n"
        "**The key idea: a handful of qubits holds a whole image.** With QPIE\n"
        "(Quantum Probability Image Encoding), an `n`-qubit state is a probability\n"
        "distribution over `2**n` bitstrings — read each bitstring as a **pixel\n"
        "address** and its probability as the pixel **brightness**. So **4 qubits\n"
        "index all 16 rows** (`2**4 = 16`), 4 more index all 16 columns, and just\n"
        "**8 qubits encode the entire 16×16 = 256-pixel image** — an *exponential*\n"
        "packing (8 numbers of hardware for 256 pixels).\n"
        "\n"
        "That packing is exactly *why you'd cut*: an image that needs 8 qubits can be\n"
        "run as **two independent 4-qubit halves** — the rows and the columns — and\n"
        "then reassembled. In this notebook you'll encode an image, apply a quantum\n"
        "transform, **cut the 8-qubit circuit into two 4-qubit pieces**, and rebuild\n"
        "the full image with `gate_cutting_sampling` (which reconstructs the whole\n"
        "output *distribution*, not just an expectation value).\n"
        "\n"
        "Here's the whole story you'll build — encode → transform → cut & reconstruct:"
    ),
    img_md("images/07_image_cutting.png", "encode, transform, cut and reconstruct"),
    md("Look for the **🔧 Your turn** cells: change a line, re-run, watch it change."),
    code(
        "import numpy as np\n"
        "import matplotlib.pyplot as plt\n"
        "from qiskit import QuantumCircuit, transpile\n"
        "from qiskit.circuit.library import StatePreparation\n"
        "import qiskit.quantum_info as qi\n"
        "import Qdislib as qd\n"
        "\n"
        "N_ROW = N_COL = 4          # 4 row qubits + 4 column qubits\n"
        "N = N_ROW + N_COL          # -> 16 x 16 = 256 pixels\n"
        "\n"
        "\n"
        "def brightness_profile(n_bits, center, width):\n"
        "    \"\"\"Amplitudes whose |.|**2 form a Gaussian bump over 2**n_bits values.\"\"\"\n"
        "    x = np.arange(2 ** n_bits)\n"
        "    p = np.exp(-0.5 * ((x - center) / width) ** 2)\n"
        "    return np.sqrt(p / p.sum())\n"
        "\n"
        "\n"
        "def to_image(distribution):\n"
        "    \"\"\"{bitstring: probability}  ->  16x16 array (row bits | column bits).\"\"\"\n"
        "    image = np.zeros((2 ** N_ROW, 2 ** N_COL))\n"
        "    for bitstring, prob in distribution.items():\n"
        "        row = int(bitstring[:N_ROW][::-1], 2)\n"
        "        col = int(bitstring[N_ROW:][::-1], 2)\n"
        "        image[row, col] += prob\n"
        "    return image\n"
        "\n"
        "\n"
        "def show(images, titles):\n"
        "    fig, axes = plt.subplots(1, len(images), figsize=(4 * len(images), 4))\n"
        "    for ax, img, title in zip(np.atleast_1d(axes), images, titles):\n"
        "        ax.imshow(img, cmap='inferno'); ax.set_title(title); ax.axis('off')\n"
        "    plt.show()"
    ),
    md(
        "## 1. Encode the image  ·  *4 qubits → 16 positions*\n"
        "\n"
        "We load a **row** brightness profile onto the 4 row qubits and a **column**\n"
        "profile onto the 4 column qubits. Each 4-qubit register spans all\n"
        "`2**4 = 16` positions along its axis. Preparing the two registers\n"
        "*independently* makes the image **separable** — a pure outer product\n"
        "`row_profile ⊗ column_profile`, i.e. a smooth blob."
    ),
    code(
        "def encode_image(transform):\n"
        "    qc = QuantumCircuit(N)\n"
        "    # QPIE: load the row and column brightness profiles (a smooth Gaussian bump).\n"
        "    qc.append(StatePreparation(brightness_profile(N_ROW, center=7.5, width=2.2)),\n"
        "              range(N_ROW))\n"
        "    qc.append(StatePreparation(brightness_profile(N_COL, center=7.5, width=2.2)),\n"
        "              range(N_ROW, N))\n"
        "    qc = transpile(qc, basis_gates=['ry', 'rz', 'cx', 'cz', 'h', 'x'],\n"
        "                   optimization_level=1)\n"
        "    if transform:\n"
        "        # TRANSFORM (step 2): fold rows into columns with two CX gates --\n"
        "        # each ties a column bit to a row bit  (col_bit ^= row_bit).\n"
        "        qc.cx(0, N_ROW)\n"
        "        qc.cx(N_ROW - 1, N - 1)\n"
        "    return qc\n"
        "\n"
        "\n"
        "def exact_image(qc):\n"
        "    \"\"\"The true image, from the exact statevector (our ground truth).\"\"\"\n"
        "    probs = np.abs(qi.Statevector(qc).data) ** 2\n"
        "    dist = {''.join(str((i >> q) & 1) for q in range(N)): probs[i]\n"
        "            for i in range(2 ** N)}\n"
        "    return to_image(dist)\n"
        "\n"
        "\n"
        "show([exact_image(encode_image(transform=False))], ['separable image (a blob)'])"
    ),
    md(
        "## 2. Transform: fold rows into columns\n"
        "\n"
        "Now we couple the two registers with two `CX` (controlled-NOT) gates, each\n"
        "from a **row qubit to a column qubit**. A `CX` sets `col_bit ^= row_bit`, so\n"
        "the **column address is XOR-shifted by the row address** — each row moves\n"
        "horizontally by an amount that depends on the row, which visibly **folds**\n"
        "the image (the blob wraps into an arc).\n"
        "\n"
        "Crucially, this coupling **entangles** the row and column registers: the\n"
        "image is no longer a separable outer product, so it genuinely needs all 8\n"
        "qubits *together* — which is what makes cutting it non-trivial. (These two\n"
        "`CX` gates double as our cut points in step 3.)"
    ),
    code(
        "show([exact_image(encode_image(transform=False)),\n"
        "      exact_image(encode_image(transform=True))],\n"
        "     ['before (separable blob)', 'after (CX fold: rows -> columns)'])"
    ),
    md(
        "### 🔧 Your turn (1)\n"
        "\n"
        "Change the **blob**: move its center or make it wider/narrower, and re-run\n"
        "the cell below. Try `center=4.0` (top-left) or `width=1.0` (sharper)."
    ),
    code(
        "def brightness_profile(n_bits, center, width):        # <-- edit center / width\n"
        "    x = np.arange(2 ** n_bits)\n"
        "    p = np.exp(-0.5 * ((x - center) / width) ** 2)\n"
        "    return np.sqrt(p / p.sum())\n"
        "\n"
        "show([exact_image(encode_image(transform=True))], ['your transformed image'])"
    ),
    md(
        "## 3. The image is an 8-qubit circuit — cut it in two\n"
        "\n"
        "Suppose our device only has ~4 qubits. `find_cut` looks for a good place to\n"
        "split the circuit; here it finds the two row↔column `CZ` gates, which\n"
        "separate the **row half** from the **column half**."
    ),
    code(
        "transformed = encode_image(transform=True)\n"
        "cut = qd.find_cut(transformed, gate_cut=True, wire_cut=False, max_qubits=N_ROW)\n"
        "print('cut gates:', cut)   # the two row<->column entanglers"
    ),
    md(
        "## 4. Reconstruct the image by sampling the pieces\n"
        "\n"
        "`gate_cutting_sampling` samples each subcircuit independently and recombines\n"
        "them (via the gate-cut quasiprobability decomposition) into the **full\n"
        "output distribution** — the reconstructed image. Two cuts means `6**2 = 36`\n"
        "subcircuit configurations, each a small 4-qubit circuit."
    ),
    code(
        "distribution = qd.gate_cutting_sampling(transformed, cut, shots=200_000,\n"
        "                                        software='qiskit')\n"
        "reconstructed = to_image(distribution)\n"
        "\n"
        "show([exact_image(transformed), reconstructed],\n"
        "     ['transformed (exact)', 'reconstructed from the 2 halves'])\n"
        "\n"
        "tvd = 0.5 * np.abs(exact_image(transformed) - reconstructed).sum()\n"
        "print(f'total-variation distance vs exact: {tvd:.4f}   (~shot noise)')"
    ),
    md(
        "### 🔧 Your turn (2)\n"
        "\n"
        "Fewer shots ⇒ a noisier reconstruction. Re-run the reconstruction with\n"
        "`shots=20_000` and compare the image and the TVD. (Sampling reconstruction\n"
        "carries a quasiprobability overhead, so it needs more shots than a plain\n"
        "expectation value.)"
    ),
    code(
        "distribution = qd.gate_cutting_sampling(transformed, cut, shots=20_000,   # <-- change\n"
        "                                        software='qiskit')\n"
        "noisy = to_image(distribution)\n"
        "show([reconstructed, noisy], ['200k shots', '20k shots'])\n"
        "print('TVD:', round(0.5 * np.abs(exact_image(transformed) - noisy).sum(), 4))"
    ),
    md(
        "## 5. Run it distributed with PyCOMPSs\n"
        "\n"
        "Each of the 36 subcircuit configurations is an **independent task**, so on a\n"
        "cluster they sample **in parallel** — the call is exactly the same. Put the\n"
        "code in a script and launch it with `runcompss`; Qdislib submits each\n"
        "subcircuit as a PyCOMPSs `@task` and waits on the results for you."
    ),
    illus(
        "# my_image_app.py  --  run with:  runcompss --python my_image_app.py\n"
        "import Qdislib as qd\n"
        "\n"
        "transformed = encode_image(transform=True)\n"
        "cut = qd.find_cut(transformed, gate_cut=True, wire_cut=False, max_qubits=4)\n"
        "# The 36 subcircuits sample in parallel across the available workers:\n"
        "distribution = qd.gate_cutting_sampling(transformed, cut, shots=200_000,\n"
        "                                        software='qiskit')"
    ),
    md(
        "## Recap & bonus\n"
        "\n"
        "You encoded an image as a quantum distribution, transformed it by\n"
        "entangling rows and columns, **cut the 8-qubit circuit into two 4-qubit\n"
        "halves**, and rebuilt the image by sampling — reconstructing the full\n"
        "distribution, not just one number.\n"
        "\n"
        "**🔧 Bonus:** bump the resolution to `32×32` by setting `N_ROW = N_COL = 5`\n"
        "(re-run from the setup cell). The circuit grows to 10 qubits and each half\n"
        "to 5 — but it's still just **two cuts** (the cost of cutting scales with the\n"
        "number of cuts, `6**n_cuts`, not the circuit size)."
    ),
]


def main():
    metadata = {
        "kernelspec": {
            "display_name": "Python 3",
            "language": "python",
            "name": "python3",
        },
        "language_info": {"name": "python"},
    }
    for name, cells in NOTEBOOKS.items():
        notebook = new_notebook(cells=cells, metadata=metadata)
        path = os.path.join(NB_DIR, f"{name}.ipynb")
        with open(path, "w") as handle:
            nbf.write(notebook, handle)
        print(f"wrote {path}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
