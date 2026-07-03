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


NOTEBOOKS = {}  # the numbered learning path -> notebooks/
APPLICATIONS = {}  # standalone application notebooks -> notebooks/applications/

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
    md(
        "## 1. Generate the subcircuits\n\n`find_cut` picks the cut; `gate_cutting_subcircuits` produces the pieces."
    ),
    code(
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        'print("cut:", cut)\n'
        "\n"
        'subcircuits = qd.gate_cutting_subcircuits(example_circuit(), cut, software="qibo")\n'
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
    md(
        "## Visualise the cut\n\n`draw_cut` shows the circuit DAG with the cut highlighted."
    ),
    illus("qd.draw_cut(example_circuit(), cut)  # needs matplotlib; opens a figure"),
    md(
        "## Wire cutting — same three steps\n"
        "\n"
        "Only the generator function changes."
    ),
    code(
        "circuit = example_circuit()\n"
        "wire_cut = qd.find_cut(circuit, gate_cut=False, wire_cut=True)\n"
        'subcircuits = qd.wire_cutting_subcircuits(example_circuit(), wire_cut, software="qibo")\n'
        "values = [qd.expectation_value(s, shots=8192) for s in subcircuits]\n"
        'print("wire reconstruction:", subcircuits.reconstruct(values))'
    ),
    md(
        "**Gate vs wire.** Gate cutting splits a two-qubit gate (cost grows fast\n"
        "with the number of cuts); wire cutting severs a qubit wire. `find_cut`\n"
        "can weigh both — see notebook 3 onwards and the docs for guidance.\n"
        "\n"
        "---\n"
        "Next: **03 · Finding cuts** — let Qdislib choose where to cut."
    ),
]

# --------------------------------------------------------------------------- #
# 03 — Finding cuts
# --------------------------------------------------------------------------- #
NOTEBOOKS["03_finding_cuts"] = [
    md(
        "# Qdislib — 03 · Finding cuts\n"
        "\n"
        "`find_cut` decides **where** to cut. You steer it with a few knobs — the device\n"
        "size (`max_qubits`), whether to look for **gate** or **wire** cuts, a cap on the\n"
        "number of cuts (`max_cuts`), and the partitioning `implementation` — and it\n"
        "returns the gates to cut, e.g. `['CZ_3']`.\n"
        "\n"
        "```python\n"
        "find_cut(circuit, max_qubits=None, max_cuts=None,\n"
        "         wire_cut=True, gate_cut=True, implementation='qdislib', weights=None)\n"
        "```"
    ),
    code(EXAMPLE),
    code("import Qdislib as qd"),
    md(
        "## Gate cuts vs wire cuts\n"
        "\n"
        "A **gate cut** splits a two-qubit gate; a **wire cut** severs a qubit wire\n"
        "between gates. `find_cut` searches both by default — restrict it with\n"
        "`gate_cut=` / `wire_cut=` (at least one must be enabled, or it raises)."
    ),
    code(
        "circuit = example_circuit()\n"
        'print("gate only:", qd.find_cut(circuit, gate_cut=True, wire_cut=False))\n'
        'print("wire only:", qd.find_cut(circuit, gate_cut=False, wire_cut=True))\n'
        'print("both     :", qd.find_cut(circuit, gate_cut=True, wire_cut=True))'
    ),
    md(
        "## Fit your device: `max_qubits`\n"
        "\n"
        "`max_qubits` caps the qubits per subcircuit — the key hardware constraint. A\n"
        "smaller device forces the circuit into more (smaller) pieces: `find_cut` aims\n"
        "for about `ceil(n_qubits / max_qubits)` components."
    ),
    code(
        "for mq in (4, 3, 2):\n"
        "    cut = qd.find_cut(example_circuit(), gate_cut=True, wire_cut=False, max_qubits=mq)\n"
        '    print(f"max_qubits={mq}: {cut}")'
    ),
    md(
        "## Cap the overhead: `max_cuts`\n"
        "\n"
        "Each gate cut multiplies the reconstruction cost by 6 (`6**n_cuts`), so\n"
        "`max_cuts` bounds how many cuts `find_cut` is allowed to use."
    ),
    code(
        "gc = lambda **kw: qd.find_cut(example_circuit(), gate_cut=True, wire_cut=False, **kw)\n"
        'print("max_cuts=1:      ", gc(max_cuts=1))\n'
        'print("max_qubits=2:    ", gc(max_qubits=2))    # wants 2 cuts (two 2-3 qubit pieces)\n'
        'print("both (too tight):", gc(max_qubits=2, max_cuts=1))  # no cut satisfies both -> []'
    ),
    md(
        "## Partitioner: `implementation`\n"
        "\n"
        "`implementation='qdislib'` (default) scores several graph-partitioning\n"
        "strategies — Kernighan–Lin, Girvan–Newman, spectral, METIS — with a loss\n"
        "function and keeps the best. `implementation='ibm-ckt'` instead uses Qiskit's\n"
        "hardware-aware cutter (gate cuts only; needs the `qiskit` extra and a Qiskit\n"
        "circuit):"
    ),
    illus(
        "from qiskit import QuantumCircuit\n"
        "qc = QuantumCircuit(10)\n"
        "for i in range(9):\n"
        "    qc.cz(i, i + 1)\n"
        "cut = qd.find_cut(qc, max_qubits=5, gate_cut=True, wire_cut=False,\n"
        "                  implementation='ibm-ckt')\n"
        "print(cut)"
    ),
    md(
        "### 🔧 Your turn\n"
        "\n"
        "Turn on `qd.set_log_level('DEBUG')` and re-run `find_cut` — it logs the\n"
        "candidate partitions, their loss scores, and which strategy won. (Advanced:\n"
        "pass `weights={...}` to retune that loss, which trades off the number of cuts,\n"
        "the number of components, and the largest component's size.)"
    ),
    code(
        'qd.set_log_level("DEBUG")\n'
        "cut = qd.find_cut(example_circuit(), gate_cut=True, wire_cut=False, max_qubits=3)\n"
        'qd.set_log_level("WARNING")   # quiet again\n'
        'print("chosen cut:", cut)'
    ),
    md(
        "## Visualise the cut: `draw_cut`\n"
        "\n"
        "`draw_cut(circuit, cut, filename)` renders the circuit's DAG with the cut gates\n"
        "highlighted, so you can see exactly where the pieces split."
    ),
    code(
        "cut = qd.find_cut(example_circuit(), gate_cut=True, wire_cut=False)\n"
        'qd.draw_cut(example_circuit(), cut, "cut.png")   # writes cut.png\n'
        'print("wrote cut.png with the cut highlighted")'
    ),
    md(
        "---\n"
        "Next: **04 · Backends & QASM interop** — run the pieces on Qibo, Qiskit or\n"
        "CUDA-Q, and bring circuits in from other frameworks."
    ),
]

# --------------------------------------------------------------------------- #
# 04 — Backends
# --------------------------------------------------------------------------- #
NOTEBOOKS["04_backends"] = [
    md(
        "# Qdislib — 04 · Backends and interop\n"
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
    code(
        'print("qibo  :", qd.gate_cutting(example_circuit(), cut, software="qibo", shots=8192))'
    ),
    md(
        "## Qiskit\n\nSame circuit and cut, executed on Qiskit Aer (`pip install Qdislib[qiskit]`)."
    ),
    code(
        'print("qiskit:", qd.gate_cutting(example_circuit(), cut, software="qiskit", shots=8192))'
    ),
    md(
        "## CUDA-Q\n"
        "\n"
        "With `pip install Qdislib[cudaq]` you can build the circuit as a CUDA-Q\n"
        'kernel and run with `software="cudaq"` (GPU when available, CPU otherwise).\n'
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
    md("---\n" "Next: **05 · Observables** — beyond `<ZZZZZ>`."),
]

# --------------------------------------------------------------------------- #
# 05 — Observables
# --------------------------------------------------------------------------- #
NOTEBOOKS["05_observables"] = [
    md(
        "# Qdislib — 05 · Observables\n"
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
        "> **Identity on a *cut* qubit.** An `I` on a qubit that a gate cut passes\n"
        "> through is the subtle case: there is no basis measurement on that qubit to\n"
        "> fix the sign of each reconstruction term, so Qdislib recovers the sign from\n"
        "> the cut's own mid-circuit measurement. You don't do anything special — it\n"
        "> still returns the exact value — but it's why measuring identities across a\n"
        "> cut works at all."
    ),
    md(
        "## Weighted sums (Hamiltonians)\n"
        "\n"
        "Pass a list of `(coefficient, pauli)` pairs or a `{pauli: coefficient}`\n"
        "mapping; Qdislib evaluates `<Σ cᵢ Pᵢ> = Σ cᵢ <Pᵢ>` term by term."
    ),
    code(
        'hamiltonian = [(0.5, "ZZIII"), (-1.0, "IIZZI"), (0.25, "ZIZII")]\n'
        "value = qd.gate_cutting(example_circuit(), cut, observables=hamiltonian, shots=20000)\n"
        "exact = qd.analytical_solution(example_circuit(), hamiltonian)\n"
        'print(f"<H>:  cut={value:+.3f}   exact={exact:+.3f}")'
    ),
    md(
        "The same `observables=` argument works for `wire_cutting`, `estimate` and\n"
        "the three-step `*_subcircuits` workflow.\n"
        "\n"
        "---\n"
        "Next: **06 · Sampling the output distribution**."
    ),
]

# --------------------------------------------------------------------------- #
# 06 — Sampling the output distribution
# --------------------------------------------------------------------------- #
NOTEBOOKS["06_sampling"] = [
    md(
        "# Qdislib — 06 · Sampling the output distribution\n"
        "\n"
        "So far a cut gave back one **expectation value** `<O>`. Sometimes you want the\n"
        "whole **output distribution** `{bitstring: probability}` — to read out an image\n"
        "(application 1) or inspect the state. `gate_cutting_sampling` and\n"
        "`wire_cutting_sampling` reconstruct that full distribution from the pieces.\n"
        "\n"
        "The trade-off: rebuilding a distribution carries a **quasiprobability overhead**,\n"
        "so it needs more shots than a plain expectation value for the same accuracy."
    ),
    code(EXAMPLE),
    code(
        "import numpy as np\n"
        "import qibo\n"
        "import Qdislib as qd\n"
        "\n"
        "qibo.set_backend('numpy')\n"
        "\n"
        "\n"
        "def exact_distribution(circuit):\n"
        '    """The true output distribution |amplitude|**2 (qubit-0-first bitstrings)."""\n'
        "    probs = np.abs(circuit().state()) ** 2\n"
        "    n = circuit.nqubits\n"
        "    return {format(i, f'0{n}b'): float(probs[i]) for i in range(2 ** n)}\n"
        "\n"
        "\n"
        "def tvd(a, b):\n"
        '    """Total-variation distance between two {bitstring: prob} distributions."""\n'
        "    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))"
    ),
    md(
        "## Expectation value vs full distribution\n"
        "\n"
        "`gate_cutting` returns a single number; `gate_cutting_sampling` returns the\n"
        "whole distribution over output bitstrings."
    ),
    code(
        "circuit = example_circuit()\n"
        "cut = qd.find_cut(circuit, gate_cut=True, wire_cut=False)\n"
        "\n"
        "value = qd.gate_cutting(example_circuit(), cut, observables='ZZZZZ', shots=20000)\n"
        "print('expectation <ZZZZZ>:', round(value, 3))\n"
        "\n"
        "dist = qd.gate_cutting_sampling(example_circuit(), cut, shots=200_000, software='qibo')\n"
        "top = sorted(dist.items(), key=lambda kv: -kv[1])[:4]\n"
        "print('top bitstrings:   ', [(b, round(p, 3)) for b, p in top])"
    ),
    md(
        "## It matches the exact distribution\n"
        "\n"
        "Reconstructing from the pieces reproduces the full statevector distribution, up\n"
        "to shot noise — the total-variation distance is ~0."
    ),
    code(
        "exact = exact_distribution(example_circuit())\n"
        "print('gate-cut sampling TVD vs exact:', round(tvd(exact, dist), 3))"
    ),
    md(
        "## The three-step sampling workflow\n"
        "\n"
        "Sampling has the same three-step form — generate the sampling subcircuits,\n"
        "`sample_subcircuit` each (a parallel COMPSs task under PyCOMPSs), then\n"
        "`reconstruct` the distribution:"
    ),
    code(
        "subs = qd.gate_cutting_sampling_subcircuits(example_circuit(), cut, software='qibo')\n"
        "sampled = [qd.sample_subcircuit(s, shots=200_000) for s in subs]\n"
        "print('3-step TVD vs exact:', round(tvd(exact, subs.reconstruct(sampled)), 3))"
    ),
    md(
        "## Wire-cut sampling\n"
        "\n"
        "Wire cuts reconstruct the distribution too — `wire_cutting_sampling` (and its\n"
        "`wire_cutting_sampling_subcircuits` three-step form):"
    ),
    code(
        "wcut = qd.find_cut(example_circuit(), gate_cut=False, wire_cut=True)\n"
        "wdist = qd.wire_cutting_sampling(example_circuit(), wcut, shots=200_000, software='qibo')\n"
        "print('wire-cut sampling TVD vs exact:', round(tvd(exact, wdist), 3))"
    ),
    md(
        "### 🔧 Your turn\n"
        "\n"
        "Fewer shots ⇒ a noisier distribution. Re-run with `shots=20_000` and watch the\n"
        "TVD grow — the quasiprobability overhead is why sampling wants more shots than\n"
        "an expectation value."
    ),
    code(
        "noisy = qd.gate_cutting_sampling(example_circuit(), cut, shots=20_000, software='qibo')\n"
        "print('TVD at 20k shots:', round(tvd(exact, noisy), 3))"
    ),
    md(
        "---\n"
        "Next: **07 · Caching & performance** — reuse repeated subcircuits with the\n"
        "semantic cache."
    ),
]

# --------------------------------------------------------------------------- #
# 07 — Caching & performance
# --------------------------------------------------------------------------- #
NOTEBOOKS["07_caching"] = [
    md(
        "# Qdislib — 07 · Caching & performance\n"
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
        'You choose the store: `"memory"` (here), `"lmdb"` (persistent on disk),\n'
        'or `"redis"` (which becomes the BSC COMPSs/PSCO store under PyCOMPSs).\n'
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
        "## Why it's *semantic*\n"
        "\n"
        "The key isn't a hash of the gate list — `semantic_key` hashes the circuit's\n"
        "**ZX-reduced** form, so two circuits that realise the *same unitary* collide\n"
        "even when written differently. `H·H·H` reduces to a single `H`, so the two\n"
        "share a key and the second is served from cache."
    ),
    code(
        "import qibo\n"
        "from qibo import Circuit, gates\n"
        "qibo.set_backend('numpy')\n"
        "\n"
        "one = Circuit(2)\n"
        "one.add(gates.H(0)); one.add(gates.CZ(0, 1)); one.add(gates.M(0, 1))\n"
        "\n"
        "three = Circuit(2)\n"
        "for _ in range(3):\n"
        "    three.add(gates.H(0))          # H·H·H == H\n"
        "three.add(gates.CZ(0, 1)); three.add(gates.M(0, 1))\n"
        "\n"
        "print('same semantic key:', qd.semantic_key(one) == qd.semantic_key(three))"
    ),
    md(
        "## Cache granularity\n"
        "\n"
        '`cache_mode="expectation"` (default) caches the per-subcircuit value keyed\n'
        'by hash + observable. `cache_mode="counts"` caches raw measurement counts\n'
        "keyed by the hash alone, so one simulation serves several observables."
    ),
    code(
        'cache = qd.default_circuit_cache("memory")\n'
        'for obs in ["ZZZZZ", "ZIZIZ"]:\n'
        "    qd.gate_cutting(example_circuit(), cut, observables=obs,\n"
        '                    cache=cache, cache_mode="counts", shots=8192)\n'
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
        "Next: **08 · Backends under the hood** — how a backend is chosen, and how to\n"
        "add your own."
    ),
]

# --------------------------------------------------------------------------- #
# 08 — Backends under the hood
# --------------------------------------------------------------------------- #
NOTEBOOKS["08_backend_internals"] = [
    md(
        "# Qdislib — 08 · Backends under the hood\n"
        "\n"
        "Every execution target — a CPU simulator (Qibo / Qiskit / CUDA-Q), a GPU, or a\n"
        "real QPU — is a `BaseBackend`. The cutting algorithms never branch on the\n"
        "backend: they hand each subcircuit to a backend object that\n"
        "`build_execution_backend` chose from your `software=` / `gpu=` / `qpu=`\n"
        "arguments. This notebook opens that box — and shows how to plug in **your own**\n"
        "backend."
    ),
    code("import Qdislib as qd"),
    md(
        "## The registry\n"
        "\n"
        "Backends are registered by name and imported only on first use, so a missing\n"
        "SDK never breaks `import Qdislib`."
    ),
    code("print(qd.available_backends())"),
    md(
        "## How `software=` picks a backend\n"
        "\n"
        "`build_execution_backend` turns `(software, gpu, qpu, qpu_dict, num_qubits)` into\n"
        "one concrete backend — the single place that decides. `software=` maps to a CPU\n"
        "backend; `qpu=True` + a big-enough QPU in `qpu_dict` routes to that QPU; a\n"
        "large-enough `gpu=True` routes to GPU Aer."
    ),
    code(
        "for sw in ('qiskit', 'qibo', 'cudaq'):\n"
        "    backend = qd.build_execution_backend(sw, num_qubits=4)\n"
        "    print(f'{sw:7s} -> {backend.name}')\n"
        "\n"
        "# qpu=True + a big-enough QPU in qpu_dict routes to the QPU backend instead:\n"
        "qpu = qd.build_execution_backend('qibo', num_qubits=4, qpu=True, qpu_dict={'MN_Ona': 8})\n"
        "print('qpu     ->', qpu.name)"
    ),
    md(
        "## The `BaseBackend` interface\n"
        "\n"
        "A backend implements a few hooks: `transpile`, `_execute_samples`,\n"
        "`_execute_estimator` (the generic sample/estimator path), plus the **cutting\n"
        "seams** the algorithms call — `evaluate_expectation` (⟨O⟩ of a subcircuit),\n"
        "`evaluate_samples` (counts), and `evaluate_joint_counts` (joint final +\n"
        "cut-point counts, for sampling reconstruction). Any backend that implements\n"
        "these can run cut subcircuits."
    ),
    md(
        "## Plug in your own backend\n"
        "\n"
        "Subclass `BaseBackend`, implement the hooks (wrap an in-house simulator or a\n"
        "cloud QPU's API), register it by name, then select it with `software='...'` —\n"
        "with no change to the cutting code:"
    ),
    illus(
        "from Qdislib import BaseBackend, register_backend\n"
        "\n"
        "class MyBackend(BaseBackend):\n"
        "    def __init__(self, config):\n"
        "        super().__init__('MY_BACKEND', config)\n"
        "    def transpile(self, circuit):\n"
        "        ...\n"
        "    def _execute_samples(self, circuit, transpile, shots, job):\n"
        "        ...   # -> SampleJobResult(counts=...)\n"
        "    def _execute_estimator(self, circuit, transpile, observable, shots, job):\n"
        "        ...   # -> ExpectationJobResult(expectation_value=...)\n"
        "    def evaluate_expectation(self, dag, num_qubits, prep_qubits=(), shots=1024,\n"
        "                             method='automatic'):\n"
        "        ...   # build the native circuit from the DAG, run it, return <O>\n"
        "\n"
        "# Registered by import path + class name (imported lazily on first use):\n"
        "register_backend('MY_BACKEND', 'my_pkg.my_backend', 'MyBackend')\n"
        "value = qd.gate_cutting(circuit, cut, software='MY_BACKEND', shots=8192)"
    ),
    md(
        "## The unified front-end: `QdBackEnd`\n"
        "\n"
        "`QdBackEnd` wraps a pool of backends and load-balances by shortest queue and\n"
        "qubit fit — handy when you have several devices. It's optional; the cutting API\n"
        "calls `build_execution_backend` directly."
    ),
    illus(
        "backend = qd.QdBackEnd(qpu_dict={'IBM_Quantum': 127, 'MN_Ona': 5})\n"
        "counts = backend.sample_circuit(circuit, shots=1024)          # routed to the best-fit device\n"
        "value = backend.get_expectation_value(circuit, 'ZZ', shots=1024)"
    ),
    md(
        "---\n"
        "Next: **09 · Distributed execution & real QPUs** — run the pieces in parallel\n"
        "with PyCOMPSs, and on IBM / BSC hardware."
    ),
]

# --------------------------------------------------------------------------- #
# 09 — Distributed & QPU (illustrative)
# --------------------------------------------------------------------------- #
NOTEBOOKS["09_distributed_and_qpu"] = [
    md(
        "# Qdislib — 09 · Distributed execution & real QPUs\n"
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
        'subcircuits = qd.gate_cutting_subcircuits(circuit, cut, software="qibo")\n'
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
        'value = qd.gate_cutting(circuit, cut, software="qiskit",\n'
        '                        qpu=True, qpu_dict={"IBM_Quantum": 127})\n'
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
        'value = qd.gate_cutting(circuit, cut, software="qibo",\n'
        '                        qpu=True, qpu_dict={"MN_Ona": 5})\n'
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
APPLICATIONS["01_image_reconstruction"] = [
    md(
        "# Qdislib — Application 1 · Cutting a quantum image, distributed with PyCOMPSs\n"
        "\n"
        "**A handful of qubits holds a whole image.** With QPIE (Quantum Probability\n"
        "Image Encoding) an `n`-qubit state is a probability distribution over `2**n`\n"
        "bitstrings — read each bitstring as a **pixel address** and its probability as\n"
        "the pixel **brightness**. So **4 qubits index all 16 rows** (`2**4 = 16`), 4\n"
        "more index all 16 columns, and just **8 qubits encode the whole 16×16 = 256\n"
        "pixel image** — an *exponential* packing (8 qubits for 256 pixels).\n"
        "\n"
        "That packing is exactly *why you'd cut*: an image that needs 8 qubits can be\n"
        "run as **two independent 4-qubit halves** — the rows and the columns — and\n"
        "reassembled. In this notebook you encode an image, **entangle** the two\n"
        "registers, **cut the 8-qubit circuit into two 4-qubit pieces**, and rebuild the\n"
        "full image with `gate_cutting_sampling` (which reconstructs the whole output\n"
        "*distribution*, not just one number). **The whole notebook runs under\n"
        "PyCOMPSs**, so the subcircuits sample in parallel.\n"
        "\n"
        "Here's the story you'll build — encode → entangle → cut & reconstruct:"
    ),
    img_md(
        "applications/images/image_cutting.png",
        "encode, entangle, cut and reconstruct an image",
    ),
    md("Look for the **🔧 Your turn** cells: change a line, re-run, watch it change."),
    md(
        "## Start the PyCOMPSs runtime\n"
        "\n"
        "This starts the COMPSs runtime **for the whole notebook**: from here on every\n"
        "`gate_cutting_sampling` call submits its subcircuits as parallel COMPSs\n"
        "`@task`s and gathers the results for you — the code is identical to running it\n"
        "serially. `graph=True` records the task graph under `$HOME/.COMPSs/…/monitor/`\n"
        "(render it with `compss_gengraph`). Run this notebook top to bottom; the very\n"
        "last cell stops the runtime."
    ),
    code(
        "import pycompss.interactive as ipycompss\n"
        "ipycompss.start(graph=True, monitor=1000)"
    ),
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
        '    """Amplitudes whose |.|**2 form a Gaussian bump over 2**n_bits values."""\n'
        "    x = np.arange(2 ** n_bits)\n"
        "    p = np.exp(-0.5 * ((x - center) / width) ** 2)\n"
        "    return np.sqrt(p / p.sum())\n"
        "\n"
        "\n"
        "def to_image(distribution):\n"
        '    """{bitstring: probability}  ->  16x16 array (row bits | column bits)."""\n'
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
        "profile onto the 4 column qubits (`StatePreparation`). Each 4-qubit register\n"
        "spans all `2**4 = 16` positions along its axis — so **8 qubits hold the whole\n"
        "16×16 image**. Preparing the two registers independently makes a smooth,\n"
        "separable blob (`row_profile ⊗ column_profile`)."
    ),
    code(
        "def encode_image(entangle, transform=None, center=7.5, profile='blob'):\n"
        "    qc = QuantumCircuit(N)\n"
        "    # QPIE: load a row profile and a column profile. 'blob' is the smooth\n"
        "    # Gaussian bump; 'ramp' is a directional gradient -- asymmetric, so a\n"
        "    # reflection visibly reverses it instead of merely relocating a symmetric blob.\n"
        "    def amps(n_bits):\n"
        "        if profile == 'ramp':\n"
        "            p = np.arange(2 ** n_bits) + 1.0        # dark -> bright gradient\n"
        "            return np.sqrt(p / p.sum())\n"
        "        return brightness_profile(n_bits, center=center, width=2.2)\n"
        "    qc.append(StatePreparation(amps(N_ROW)), range(N_ROW))\n"
        "    qc.append(StatePreparation(amps(N_COL)), range(N_ROW, N))\n"
        "    qc = transpile(qc, basis_gates=['ry', 'rz', 'cx', 'cz', 'h', 'x'],\n"
        "                   optimization_level=1)\n"
        "    # GEOMETRIC TRANSFORM (step 2b): an X on an address qubit flips that\n"
        "    # coordinate bit; flipping every bit maps index k -> (2**n - 1) - k, a clean\n"
        "    # mirror. X acts inside one register, so it never crosses the row<->column\n"
        "    # cut -- the transform reconstructs faithfully from the two halves.\n"
        "    if transform in ('flip_rows', 'rotate_180'):\n"
        "        for q in range(N_ROW):        # mirror rows: top <-> bottom\n"
        "            qc.x(q)\n"
        "    if transform in ('flip_cols', 'rotate_180'):\n"
        "        for q in range(N_ROW, N):     # mirror columns: left <-> right\n"
        "            qc.x(q)\n"
        "    if entangle:\n"
        "        # ENTANGLE (step 2): couple the two registers with two CZ gates, each\n"
        "        # tying a row qubit to a column qubit. These are our cut points.\n"
        "        qc.cz(0, N_ROW)\n"
        "        qc.cz(N_ROW - 1, N - 1)\n"
        "    return qc\n"
        "\n"
        "\n"
        "def exact_image(qc):\n"
        '    """The true image, from the exact statevector (our ground truth)."""\n'
        "    probs = np.abs(qi.Statevector(qc).data) ** 2\n"
        "    dist = {''.join(str((i >> q) & 1) for q in range(N)): probs[i]\n"
        "            for i in range(2 ** N)}\n"
        "    return to_image(dist)\n"
        "\n"
        "\n"
        "show([exact_image(encode_image(entangle=False))], ['encoded image (a blob)'])"
    ),
    md(
        "## 2. Entangle the rows and columns\n"
        "\n"
        "Now we couple the two registers with two `CZ` (controlled-Z) gates, each from\n"
        "a **row qubit to a column qubit**. This genuinely **entangles** the 8-qubit\n"
        "state — it is no longer a product of an independent row state and column\n"
        "state, so you can't just run the two 4-qubit halves separately and multiply.\n"
        "\n"
        "A `CZ` only adds a **phase**, so it leaves the measurement probabilities — the\n"
        "**image** — exactly unchanged. That's the point of this demo: the picture\n"
        "stays put while the circuit becomes entangled, so when we cut and rebuild it,\n"
        "a correct reconstruction must return the **identical** blob. (These two `CZ`\n"
        "gates are the cut points in step 3.)"
    ),
    code(
        "show([exact_image(encode_image(entangle=False)),\n"
        "      exact_image(encode_image(entangle=True))],\n"
        "     ['before (separable)', 'after CZ (entangled -- image unchanged)'])"
    ),
    md(
        "### 🔧 Your turn (1)\n"
        "\n"
        "Change the **blob**: move its center or make it wider/narrower, then re-run.\n"
        "Try `center=4.0` (top-left) or `width=1.0` (sharper)."
    ),
    code(
        "def brightness_profile(n_bits, center, width):        # <-- edit center / width\n"
        "    x = np.arange(2 ** n_bits)\n"
        "    p = np.exp(-0.5 * ((x - center) / width) ** 2)\n"
        "    return np.sqrt(p / p.sum())\n"
        "\n"
        "show([exact_image(encode_image(entangle=True))], ['your image'])"
    ),
    md(
        "## 2b. A geometric transform that survives the cut\n"
        "\n"
        "The `CZ` left the image untouched — a pure fidelity check. But we can also\n"
        "**transform** the image with gates and watch the transform come back through\n"
        "the cut. The clean transforms in QPIE are the **geometric** ones: an `X` on an\n"
        "address qubit flips that coordinate bit, and flipping *every* bit of a\n"
        "coordinate maps index `k → (2**n − 1) − k` — a **mirror reflection**. Each `X`\n"
        "acts inside a single register, so it never crosses the row↔column cut: the\n"
        "reflected image is reconstructed from the two halves just as faithfully. (We\n"
        "switch to a directional **gradient** here — a mirror of the *symmetric* blob\n"
        "would only relocate it, so you'd never actually see the flip.)\n"
        "\n"
        "*(Why reflect and not shrink/zoom? Scaling a QPIE image is a **resampling** —\n"
        "shrinking is 2-to-1 (non-unitary) and zooming needs an extra qubit — so it is\n"
        "not a clean gate. The easy way to resize the blob is at encode time, via\n"
        "`width`. Reflections, translations and 90° rotations are the clean ones.)*"
    ),
    code(
        "# A directional gradient makes the reflection unmistakable -- the bright side\n"
        "# reverses. (A symmetric blob would just relocate, so the flip would be\n"
        "# invisible.) transform: 'flip_rows', 'flip_cols' or 'rotate_180' (both axes).\n"
        "ramp = dict(entangle=False, profile='ramp')\n"
        "show([exact_image(encode_image(**ramp)),\n"
        "      exact_image(encode_image(**ramp, transform='flip_rows')),\n"
        "      exact_image(encode_image(**ramp, transform='flip_cols')),\n"
        "      exact_image(encode_image(**ramp, transform='rotate_180'))],\n"
        "     ['original (gradient)', 'flip rows', 'flip cols', 'rotate 180°'])"
    ),
    md(
        "### 🔧 Your turn (1b)\n"
        "\n"
        "**Cut and reconstruct a transformed image.** Pick a `TRANSFORM` and re-run: the\n"
        "reconstruction from the two 4-qubit halves matches the *transformed* image, not\n"
        "the original — the transform rode through the cut."
    ),
    code(
        "TRANSFORM = 'flip_cols'      # <-- 'flip_rows' / 'flip_cols' / 'rotate_180' / None\n"
        "\n"
        "moved = encode_image(entangle=True, profile='ramp', transform=TRANSFORM)\n"
        "cut_moved = qd.find_cut(moved, gate_cut=True, wire_cut=False, max_qubits=N_ROW)\n"
        "dist_moved = qd.gate_cutting_sampling(moved, cut_moved, shots=200_000,\n"
        "                                      software='qiskit')\n"
        "show([exact_image(moved), to_image(dist_moved)],\n"
        "     ['transformed image (exact)', 'reconstructed from the 2 halves'])"
    ),
    md(
        "## 3. The image is an 8-qubit circuit — cut it in two\n"
        "\n"
        "Suppose our device only has ~4 qubits. `find_cut` looks for a good place to\n"
        "split the circuit; here it finds the two row↔column `CZ` gates, which separate\n"
        "the **row half** from the **column half**."
    ),
    code(
        "entangled = encode_image(entangle=True)\n"
        "cut = qd.find_cut(entangled, gate_cut=True, wire_cut=False, max_qubits=N_ROW)\n"
        "print('cut gates:', cut)   # the two row<->column entanglers"
    ),
    md(
        "### The circuit we're cutting\n"
        "\n"
        "Here is the actual **transpiled** 8-qubit circuit `find_cut` worked on. Each\n"
        "`StatePreparation` has compiled down to `ry` / `rz` / `cx` **inside its\n"
        "register**, and the only gates joining the row register (`q0–q3`) to the\n"
        "column register (`q4–q7`) are the two `CZ` — the cut points. It is a concrete\n"
        "gate list, not an opaque block, which is exactly what lets it be cut."
    ),
    code(
        "print('gate counts:', dict(entangled.count_ops()), '| depth:', entangled.depth())\n"
        "# Text drawing (no extra deps). The two CZ (q0-q4 and q3-q7) are the only gates\n"
        "# spanning the row and column registers -- the cut points.\n"
        "print(entangled.draw(fold=74))"
    ),
    md(
        "## 4. Reconstruct the image by sampling the pieces — in parallel\n"
        "\n"
        "We use the **three-step workflow**, which exposes each piece (so you can\n"
        "inspect, batch, cache or place them individually):\n"
        "\n"
        "1. `gate_cutting_sampling_subcircuits` splits the cut into the sampling\n"
        "   subcircuits — two cuts ⇒ `6**2 = 36` configurations, each a small 4-qubit\n"
        "   circuit;\n"
        "2. `sample_subcircuit` samples one subcircuit — **under PyCOMPSs each call is a\n"
        "   COMPSs task**, so the pieces sample in parallel (and you decide where each\n"
        "   one runs — e.g. a QPU, step 5);\n"
        "3. `subs.reconstruct(...)` recombines them (via the gate-cut quasiprobability\n"
        "   decomposition) into the **full output distribution** — the reconstructed\n"
        "   image.\n"
        "\n"
        "The one-liner `gate_cutting_sampling(entangled, cut, ...)` just wraps these\n"
        "three steps."
    ),
    code(
        "subs = qd.gate_cutting_sampling_subcircuits(entangled, cut, software='qiskit')\n"
        "print(len(subs), 'sampling subcircuits (6**2 configs x components)')\n"
        "\n"
        "# Each sample_subcircuit call is a COMPSs task -> the pieces sample in parallel.\n"
        "sampled = [qd.sample_subcircuit(s, shots=200_000) for s in subs]\n"
        "distribution = subs.reconstruct(sampled)\n"
        "reconstructed = to_image(distribution)\n"
        "\n"
        "show([exact_image(entangled), reconstructed],\n"
        "     ['exact image', 'reconstructed from the 2 halves'])\n"
        "\n"
        "tvd = 0.5 * np.abs(exact_image(entangled) - reconstructed).sum()\n"
        "print(f'total-variation distance vs exact: {tvd:.4f}   (~shot noise)')"
    ),
    md(
        "### 🔧 Your turn (2)\n"
        "\n"
        "Fewer shots ⇒ a noisier reconstruction. Re-run with `shots=20_000` and compare\n"
        "the image and the TVD. (Sampling reconstruction carries a quasiprobability\n"
        "overhead, so it needs more shots than a plain expectation value.)"
    ),
    code(
        "distribution = qd.gate_cutting_sampling(entangled, cut, shots=20_000,   # <-- change\n"
        "                                        software='qiskit')\n"
        "noisy = to_image(distribution)\n"
        "show([reconstructed, noisy], ['200k shots', '20k shots'])\n"
        "print('TVD:', round(0.5 * np.abs(exact_image(entangled) - noisy).sum(), 4))"
    ),
    md(
        "## 5. Send the pieces to a real BSC quantum computer\n"
        "\n"
        "The subcircuits are only **4 qubits** — small enough for today's hardware. To\n"
        "sample them on a **real BSC QPU** instead of the local simulator, pass\n"
        "`qpu=True` and name the machine with `qpu_dict` (QPU → its qubit count). The\n"
        "call is otherwise identical: Qdislib ships each subcircuit to the QPU, gathers\n"
        "the counts, and rebuilds the image for you. (This needs access to the BSC\n"
        "quantum environment, so it isn't run in this docs build.)"
    ),
    illus(
        "distribution = qd.gate_cutting_sampling(\n"
        "    entangled, cut, shots=4096,\n"
        "    software='qibo', qpu=True, qpu_dict={'MN_Ona': N_ROW},\n"
        ")\n"
        "show([exact_image(entangled), to_image(distribution)],\n"
        "     ['exact image', 'reconstructed on a BSC QPU'])"
    ),
    md(
        "### What a real (noisy) run looks like\n"
        "\n"
        "On hardware the reconstruction still recovers the blob, but gate and readout\n"
        "errors blur it and add a background haze — and cutting's quasiprobability\n"
        "overhead amplifies shot noise, so hardware runs want **more shots** and **error\n"
        "mitigation**. Here is the same image reconstructed on an ideal simulator versus\n"
        "a noisy QPU run:"
    ),
    img_md(
        "applications/images/noisy_execution.png",
        "the reconstructed image on an ideal simulator versus a noisy QPU",
    ),
    md(
        "## Recap & bonus\n"
        "\n"
        "You encoded a 256-pixel image in **8 qubits**, entangled the row and column\n"
        "registers with `CZ`, **cut the 8-qubit circuit into two 4-qubit halves**, and\n"
        "rebuilt the *full image distribution* by sampling — all running distributed\n"
        "under PyCOMPSs — recovering the identical blob.\n"
        "\n"
        "**🔧 Bonus:** bump the resolution to `32×32` with `N_ROW = N_COL = 5` (re-run\n"
        "from the setup cell). The circuit grows to 10 qubits and each half to 5 — but\n"
        "it's still just **two cuts** (the cost scales with the number of cuts,\n"
        "`6**n_cuts`, not the circuit size)."
    ),
    md(
        "## Stop the PyCOMPSs runtime\n"
        "\n"
        "Run this **last**, once you're done re-running the cells above — it shuts the\n"
        "COMPSs runtime down cleanly, syncing any outstanding results."
    ),
    code("ipycompss.stop(sync=True)"),
]


# --------------------------------------------------------------------------- #
# 08 — Application: VQE ground-state energy of H2 by cutting the ansatz
# --------------------------------------------------------------------------- #
APPLICATIONS["02_vqe_h2"] = [
    md(
        "# Qdislib — Application 2 · VQE for a molecule (H₂), distributed with PyCOMPSs\n"
        "\n"
        "**Find a molecule's ground-state energy on a device too small for the "
        "circuit.** The variational quantum eigensolver (VQE) prepares a trial state "
        "`|ψ(θ)⟩` with a parameterized **ansatz** and minimizes the energy "
        "`E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩` toward the ground state. For H₂ (STO-3G basis, "
        "Jordan–Wigner mapping) the qubit Hamiltonian `H` is a weighted sum of "
        "**15 Pauli strings on 4 qubits**.\n"
        "\n"
        "The idea: **cut the 4-qubit ansatz into two 2-qubit halves**, so it runs on "
        "~2-qubit hardware. The energy is `⟨H⟩ = Σ_i c_i ⟨P_i⟩`, and each `⟨P_i⟩` is "
        "reconstructed from small subcircuits. That's **15 Hamiltonian terms × 36 cut "
        "configurations ≈ 500 independent subcircuits** — and **under PyCOMPSs they "
        "run in parallel**. Circuits are built and executed with **CUDA-Q**.\n"
        "\n"
        "Look for the **🔧 Your turn** cells."
    ),
    md(
        "## Start the PyCOMPSs runtime\n"
        "\n"
        "This starts the COMPSs runtime **for the whole notebook**: every "
        "`expectation_value` below runs as a parallel COMPSs `@task`. Run the notebook "
        "top to bottom; the last cell stops the runtime."
    ),
    code(
        "import pycompss.interactive as ipycompss\n"
        "ipycompss.start(graph=True, monitor=1000)"
    ),
    code(
        "import numpy as np\n"
        "import matplotlib.pyplot as plt\n"
        "import cudaq\n"
        "from cudaq import spin\n"
        "import Qdislib as qd\n"
        "\n"
        "# H2 / STO-3G, Jordan-Wigner, bond length 0.7414 A. The qubit Hamiltonian is\n"
        "# a weighted sum of Pauli strings (coefficients from a standard\n"
        "# quantum-chemistry calculation). Qdislib takes it directly as\n"
        "# observables=[(coeff, pauli), ...] (q0 = first character).\n"
        "HAMILTONIAN = [\n"
        '    (-0.04207897, "IIII"),\n'
        '    ( 0.17771287, "ZIII"), ( 0.17771287, "IZII"),\n'
        '    (-0.24274281, "IIZI"), (-0.24274281, "IIIZ"),\n'
        '    ( 0.17059738, "ZZII"),\n'
        '    ( 0.12293305, "ZIZI"), ( 0.16768319, "ZIIZ"),\n'
        '    ( 0.16768319, "IZZI"), ( 0.12293305, "IZIZ"),\n'
        '    ( 0.17627641, "IIZZ"),\n'
        '    ( 0.04475014, "YXXY"), (-0.04475014, "YYXX"),\n'
        '    (-0.04475014, "XXYY"), ( 0.04475014, "XYYX"),\n'
        "]\n"
        "\n"
        "\n"
        "def spin_operator(hamiltonian):\n"
        '    """The Hamiltonian as a CUDA-Q SpinOperator (for the exact energy)."""\n'
        "    paulis = {'X': spin.x, 'Y': spin.y, 'Z': spin.z}\n"
        "    total = 0\n"
        "    for coeff, pauli in hamiltonian:\n"
        "        term = None\n"
        "        for qubit, op in enumerate(pauli):\n"
        "            if op == 'I':\n"
        "                continue\n"
        "            factor = paulis[op](qubit)\n"
        "            term = factor if term is None else term * factor\n"
        "        total += coeff * (term if term is not None else spin.i(0))\n"
        "    return total\n"
        "\n"
        "\n"
        "H_OP = spin_operator(HAMILTONIAN)"
    ),
    md(
        "## 1. The ansatz — a 4-qubit CUDA-Q kernel\n"
        "\n"
        "H₂'s ground state lives in the two-configuration space spanned by the "
        "Hartree–Fock state `|1100⟩` (both electrons in the low orbital) and the "
        "doubly-excited `|0011⟩`. A **single parameter** `θ` interpolates between them "
        "— `cos(θ/2)|0011⟩ + sin(θ/2)|1100⟩` — which is enough to reach the exact "
        "ground state. We build it as a CUDA-Q kernel: an `RY(θ)` sets the amplitude, "
        "and the two `CX` gates from qubit 0 into the `{2,3}` register are the only "
        "gates that **cross the two halves** — the cut points."
    ),
    code(
        "def ansatz(theta):\n"
        "    @cudaq.kernel\n"
        "    def kernel():\n"
        "        q = cudaq.qvector(4)\n"
        "        ry(theta, q[0])       # amplitude of the excitation\n"
        "        x.ctrl(q[0], q[1])    # correlate the {0,1} half\n"
        "        x(q[2]); x(q[3])      # start the {2,3} half at |11>\n"
        "        x.ctrl(q[0], q[2])    # <-- crosses the cut\n"
        "        x.ctrl(q[0], q[3])    # <-- crosses the cut\n"
        "    return kernel\n"
        "\n"
        "# The ground-state manifold is exactly the two occupation patterns |0011>/|1100>:\n"
        "print(cudaq.sample(ansatz(1.0), shots_count=2000))"
    ),
    md(
        "## 2. The exact energy landscape\n"
        "\n"
        "Sweep `θ` and evaluate `E(θ) = ⟨H⟩` exactly with `cudaq.observe` (no cutting "
        "yet). The minimum of the curve is the **ground-state energy** — for H₂ it "
        "should land near the textbook **−1.137 Ha**."
    ),
    code(
        "thetas = np.linspace(0, 2 * np.pi, 200)\n"
        "energies = [cudaq.observe(ansatz(t), H_OP).expectation() for t in thetas]\n"
        "\n"
        "theta_opt = float(thetas[int(np.argmin(energies))])\n"
        "E_ground = float(min(energies))\n"
        "\n"
        "plt.figure(figsize=(6, 4))\n"
        "plt.plot(thetas, energies)\n"
        "plt.axvline(theta_opt, ls='--', c='grey')\n"
        "plt.xlabel('theta'); plt.ylabel('energy (Ha)')\n"
        "plt.title('H2 VQE energy landscape'); plt.show()\n"
        "\n"
        "print(f'ground-state energy: {E_ground:.4f} Ha  at theta = {theta_opt:.3f}"
        "   (H2 ~ -1.137 Ha)')"
    ),
    img_md(
        "applications/images/vqe_landscape.png",
        "H2 VQE energy landscape: E(theta) dips to the exact ground energy",
    ),
    md(
        "### 🔧 Your turn (1)\n"
        "\n"
        "The ansatz has **one** parameter, so the landscape is a single curve. Increase "
        "the scan resolution (`200` → `400`) or zoom into `theta ∈ [3.0, 3.7]` and read "
        "off the minimum more precisely."
    ),
    md(
        "## 3. Cut the ansatz into two 2-qubit halves\n"
        "\n"
        "`find_cut` looks at the ansatz and finds the two `CX` gates that couple the "
        "`{0,1}` and `{2,3}` registers — cutting them splits the 4-qubit circuit into "
        "two independent 2-qubit pieces."
    ),
    code(
        "optimal = ansatz(theta_opt)\n"
        "cut = qd.find_cut(optimal, gate_cut=True, wire_cut=False, max_qubits=2)\n"
        "print('cut gates:', cut)   # the two cross-register CX"
    ),
    md(
        "## 4. Reconstruct the energy by cutting — ~500 subcircuits in parallel\n"
        "\n"
        "We use the **three-step workflow** for each Hamiltonian term `P_i`:\n"
        "\n"
        "1. `gate_cutting_subcircuits` → the term's cut subcircuits (2 cuts ⇒ "
        "`6**2 = 36` configurations);\n"
        "2. `expectation_value(s)` → evaluate each one — **a COMPSs task** on CUDA-Q;\n"
        "3. `subs.reconstruct(values)` → the term's expectation `⟨P_i⟩`.\n"
        "\n"
        "We **submit every term's evaluations first** (Phase 1) so all "
        "`14 × 36 = 504` subcircuits (the identity term is free) schedule together and "
        "PyCOMPSs runs them in parallel; then we **reconstruct and take the weighted "
        "sum** `⟨H⟩ = Σ_i c_i ⟨P_i⟩` (Phase 2)."
    ),
    code(
        "# Phase 1 -- generate every term's subcircuits and SUBMIT all evaluations.\n"
        "#   Each expectation_value(...) is a COMPSs task; building the list schedules\n"
        "#   all subcircuits at once, so they run in parallel across workers.\n"
        "plans = []\n"
        "for coeff, pauli in HAMILTONIAN:\n"
        "    if set(pauli) == {'I'}:            # identity term: <I> = 1, no circuit\n"
        "        plans.append((coeff, None, None))\n"
        "        continue\n"
        "    subs = qd.gate_cutting_subcircuits(optimal, cut, observables=pauli,\n"
        "                                       software='cudaq')\n"
        "    values = [qd.expectation_value(s, shots=100_000) for s in subs]\n"
        "    plans.append((coeff, subs, values))\n"
        "\n"
        "# Phase 2 -- reconstruct each term and take the weighted sum <H> = sum_i c_i <P_i>.\n"
        "energy = sum(c if subs is None else c * subs.reconstruct(v)\n"
        "             for c, subs, v in plans)\n"
        "\n"
        "print(f'cut-reconstructed energy: {energy:+.4f} Ha')\n"
        "print(f'exact (statevector):      {E_ground:+.4f} Ha')\n"
        "print(f'error: {abs(energy - E_ground):.4f} Ha   (~shot noise x quasiprobability overhead)')"
    ),
    md(
        "## Recap & bonus\n"
        "\n"
        "You estimated H₂'s ground-state energy by **cutting a 4-qubit VQE ansatz into "
        "two 2-qubit halves** and reconstructing `⟨H⟩` from ~500 small subcircuits — "
        "the 15 Hamiltonian terms × 36 cut configurations — **sampled in parallel "
        "under PyCOMPSs, on CUDA-Q**, matching the exact energy to shot noise.\n"
        "\n"
        "**🔧 Your turn (2):** a full VQE just wraps this energy in a classical "
        "minimizer. Sharpen the result — raise `shots` (the error is shot-noise × the "
        "quasiprobability overhead) or refine `theta_opt` with a finer scan near the "
        "minimum — and watch the cut energy approach the textbook **−1.137 Ha**.\n"
        "\n"
        "**On real hardware:** the 2-qubit subcircuits are QPU-sized, but here the cut "
        "qubits keep evolving after the cut, so the terms that observe them with "
        "identity need the cut-point **sign decoupling** the CPU / CUDA-Q evaluators "
        "provide — aggregated QPU counts can't (yet). See **notebook 06** for the QPU "
        "API on circuits whose cut is terminal."
    ),
    md(
        "## Stop the PyCOMPSs runtime\n"
        "\n"
        "Run this **last** — it shuts the COMPSs runtime down cleanly, syncing any "
        "outstanding results."
    ),
    code("ipycompss.stop(sync=True)"),
]


# --------------------------------------------------------------------------- #
# 09 — Scaling up: a bigger molecule (N2) and a 25-qubit lattice, by cutting
# --------------------------------------------------------------------------- #
APPLICATIONS["03_scaling"] = [
    md(
        "# Qdislib — Application 3 · Scaling up: a bigger molecule and a 25-qubit lattice\n"
        "\n"
        "Application 2 cut **H₂ (4 qubits)**. The real payoff of cutting is running "
        "circuits **far too big for your device** — so here we scale up on two fronts:\n"
        "\n"
        "1. **N₂**, a real molecule — 20 qubits at full size, which we bring to **12 "
        "qubits** with an *active space*, and cut into two 6-qubit halves;\n"
        "2. a **25-qubit spin lattice**, cut into **four blocks with 3 cuts**.\n"
        "\n"
        "In both cases the energy `⟨H⟩ = Σ_i c_i ⟨P_i⟩` is reconstructed from small "
        "subcircuits that **run in parallel under PyCOMPSs**. (We use `software='qiskit'` "
        "here for speed at scale; `'cudaq'` / `'qibo'` are one-argument swaps — cutting "
        "is backend-agnostic.)\n"
        "\n"
        "*Cost intuition:* the work is `#Hamiltonian-terms × 6**(#cuts)` subcircuits — "
        "it grows with the **number of cuts and terms**, not the circuit's qubit count. "
        "That's why a 25-qubit *lattice* (few terms) is cheaper than a 24-qubit "
        "*molecule* (~15,000 terms)."
    ),
    img_md(
        "applications/images/scaling.png",
        "cut-reconstructed vs exact energy for H2 (4 qubits) and N2 (12 qubits)",
    ),
    md("## Start the PyCOMPSs runtime"),
    code(
        "import pycompss.interactive as ipycompss\n"
        "ipycompss.start(graph=True, monitor=1000)"
    ),
    code(
        "import os, json\n"
        "import numpy as np\n"
        "from qiskit import QuantumCircuit\n"
        "from qiskit.quantum_info import Statevector, SparsePauliOp\n"
        "import Qdislib as qd\n"
        "\n"
        "\n"
        "def hardware_efficient(n, seed=1):\n"
        '    """A generic n-qubit ansatz: RY layer, linear CZ entanglers, RY layer.\n'
        "    The linear entangler means a single bond crosses each block boundary, so\n"
        '    find_cut splits it with the fewest possible cuts."""\n'
        "    rng = np.random.default_rng(seed)\n"
        "    qc = QuantumCircuit(n)\n"
        "    for i in range(n):\n"
        "        qc.ry(float(rng.uniform(0, np.pi)), i)\n"
        "    for i in range(n - 1):\n"
        "        qc.cz(i, i + 1)\n"
        "    for i in range(n):\n"
        "        qc.ry(0.3, i)\n"
        "    return qc\n"
        "\n"
        "\n"
        "def exact_expectation(qc, hamiltonian):\n"
        '    """Exact <H> of a circuit via statevector (the reference to match)."""\n'
        "    op = SparsePauliOp.from_list([(p[::-1], c) for c, p in hamiltonian])\n"
        "    return float(Statevector(qc).expectation_value(op).real)"
    ),
    md(
        "## 1. A real molecule: N₂ at 12 qubits (active space)\n"
        "\n"
        "N₂ in the STO-3G basis is **20 qubits** with a Hamiltonian of *thousands* of "
        "Pauli terms — too big to reconstruct term-by-term. Chemists tame this with an "
        "**active space**: freeze the core orbitals and keep only the chemically active "
        "ones. **CAS(6,6)** (6 electrons in 6 orbitals) gives a **12-qubit**, 383-term "
        "Hamiltonian — a genuine N₂ calculation, three times H₂'s size. We load a "
        "precomputed copy (generation cell below)."
    ),
    code(
        "n2 = json.load(open(os.path.join('data', 'n2_cas66_hamiltonian.json')))\n"
        "HAM_N2 = [(coeff, pauli) for coeff, pauli in n2['terms']]\n"
        "N_N2 = n2['n_qubits']\n"
        'print(f"N2 CAS(6,6): {N_N2} qubits, {len(HAM_N2)} Pauli terms | "\n'
        "      f\"FCI reference {n2['fci_energy']:.4f} Ha\")"
    ),
    illus(
        "# How data/n2_cas66_hamiltonian.json was generated (one-off):\n"
        "#   pip install pyscf openfermion openfermionpyscf\n"
        "from openfermion import MolecularData, jordan_wigner\n"
        "from openfermionpyscf import run_pyscf\n"
        "\n"
        "geometry = [('N', (0, 0, 0)), ('N', (0, 0, 1.10))]\n"
        "molecule = run_pyscf(MolecularData(geometry, 'sto-3g', 1, 0), run_scf=1, run_fci=1)\n"
        "# CAS(6,6): freeze the 4 lowest spatial orbitals, keep 6 active -> 12 qubits.\n"
        "hamiltonian = molecule.get_molecular_hamiltonian(occupied_indices=range(4),\n"
        "                                                 active_indices=range(4, 10))\n"
        "qubit_op = jordan_wigner(hamiltonian)   # -> 383 Pauli terms"
    ),
    md(
        "Build a 12-qubit ansatz and cut it. With a linear entangler only the middle "
        "bond crosses, so `find_cut` splits it into **two 6-qubit halves with one cut** "
        "(`6**1 = 6` configurations). One energy is then `383 terms × 6 = ~2300` small "
        "subcircuits — all independent, so PyCOMPSs runs them in parallel."
    ),
    code(
        "ansatz_n2 = hardware_efficient(N_N2)\n"
        "cut_n2 = qd.find_cut(ansatz_n2, gate_cut=True, wire_cut=False, max_qubits=N_N2 // 2)\n"
        "print('N2 cut:', cut_n2, '-> two 6-qubit halves')\n"
        "\n"
        "energy_n2 = qd.gate_cutting(ansatz_n2, cut_n2, observables=HAM_N2,\n"
        "                            shots=20_000, software='qiskit')\n"
        "print(f'N2 <H>  cut-reconstructed: {energy_n2:+.3f} Ha')\n"
        "print(f'N2 <H>  exact (this ansatz): {exact_expectation(ansatz_n2, HAM_N2):+.3f} Ha')\n"
        "# The cut value tracks the exact <H> up to ~0.2 Ha: that residual is shot noise\n"
        "# summed over 383 terms x the quasiprobability overhead -- VQE's measurement cost,\n"
        "# which more shots, Pauli-term grouping and error mitigation shrink. (The ansatz\n"
        "# is not optimised; a VQE loop would minimise <H> toward the FCI energy above.)"
    ),
    md(
        "## 2. A 25-qubit lattice: cut into four blocks with 3 cuts\n"
        "\n"
        "A **transverse-field Ising** chain `H = -Σ Z_iZ_{i+1} - h Σ X_i` on **25 "
        "qubits** — way past any small device. Because it's a lattice, it has only "
        "`O(N)` terms (~49), so unlike a molecule we can push the qubit count hard. We "
        "cut the 25-qubit ansatz into **four ~6-qubit blocks with 3 cuts** "
        "(`6**3 = 216` configurations)."
    ),
    code(
        "def ising_hamiltonian(n, h=1.0):\n"
        "    H = []\n"
        "    for i in range(n - 1):                      # -Z_i Z_{i+1}\n"
        "        s = ['I'] * n; s[i] = s[i + 1] = 'Z'; H.append((-1.0, ''.join(s)))\n"
        "    for i in range(n):                          # -h X_i\n"
        "        s = ['I'] * n; s[i] = 'X'; H.append((-h, ''.join(s)))\n"
        "    return H\n"
        "\n"
        "N_SPIN = 25\n"
        "HAM_SPIN = ising_hamiltonian(N_SPIN)\n"
        "ansatz_spin = hardware_efficient(N_SPIN, seed=2)\n"
        "cut_spin = qd.find_cut(ansatz_spin, gate_cut=True, wire_cut=False, max_qubits=7)\n"
        "print(f'25-qubit Ising: {len(HAM_SPIN)} terms | {len(cut_spin)} cuts '\n"
        "      f'-> {len(cut_spin) + 1} blocks: {cut_spin}')"
    ),
    md(
        "Reconstruct the energy. This is `49 terms × 216 configs ≈ 10⁴` subcircuits (each "
        "≤ 7 qubits) — heavy serially, but embarrassingly parallel under PyCOMPSs. The "
        "exact reference is a 25-qubit statevector (~512 MB of RAM; lower `N_SPIN` if "
        "that's too much)."
    ),
    code(
        "energy_spin = qd.gate_cutting(ansatz_spin, cut_spin, observables=HAM_SPIN,\n"
        "                              shots=20_000, software='qiskit')\n"
        "print(f'25-qubit lattice <H>  cut-reconstructed: {energy_spin:+.3f}')\n"
        "print(f'25-qubit lattice <H>  exact:             {exact_expectation(ansatz_spin, HAM_SPIN):+.3f}')"
    ),
    md(
        "## Recap\n"
        "\n"
        "You reconstructed the energy of a **12-qubit N₂** molecule (cut into 2×6) and a "
        "**25-qubit lattice** (cut into 4×~6 with 3 cuts) from small subcircuits, in "
        "parallel under PyCOMPSs — circuits well beyond a single small device.\n"
        "\n"
        "The lesson from the two examples: cutting's cost is `#terms × 6**(#cuts)`. A "
        "lattice keeps `#terms` small so you can chase **qubit count** (25+); a molecule's "
        "`O(N⁴)` terms mean you chase **accuracy** at a modest size (active spaces). The "
        "frontier — full N₂ at 20 qubits, or larger molecules — is exactly where cutting, "
        "Pauli-term grouping, and distributed execution have to work together."
    ),
    md("## Stop the PyCOMPSs runtime"),
    code("ipycompss.stop(sync=True)"),
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

    def write_all(collection, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        for name, cells in collection.items():
            notebook = new_notebook(cells=cells, metadata=metadata)
            path = os.path.join(out_dir, f"{name}.ipynb")
            with open(path, "w") as handle:
                nbf.write(notebook, handle)
            print(f"wrote {path}  ({len(cells)} cells)")

    write_all(NOTEBOOKS, NB_DIR)  # the numbered learning path
    write_all(APPLICATIONS, os.path.join(NB_DIR, "applications"))  # standalone apps


if __name__ == "__main__":
    main()
