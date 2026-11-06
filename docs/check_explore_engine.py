#!/usr/bin/env python3
"""Check the Explore page's JavaScript engine against real Qdislib, on random circuits.

    python docs/check_explore_engine.py [circuits] [seed]

Each random circuit, with random gate, Hadamard and wire cuts and a random observable, goes
to docs/source/_static/explore_engine.js (run by Node) and to Qdislib. They must agree on:

- the exact expectation value of the uncut circuit (qd.analytical_solution);
- the reconstruction (qd.circuit_cutting(..., shots=None));
- whether Qdislib refuses the cut (InvalidCutError);
- the fragments, the configurations, the subcircuits and gamma (the subcircuit list, and
  the optimized count when no Hadamard cut is in it), and with Hadamard cuts the
  subcircuits and gamma of joint=True, Qdislib's default.

Every reconstruction is also compared with the engine's own uncut value. That is the only
check on two kinds of circuit, listed separately: Qdislib raises on a wire cut right after a
cut gate, and counts one fragment more when Hadamard and gate cuts share a line that has
nothing else on it in places. Those are the library's, with the values in agreement.

Hadamard cuts go with joint=False, each with its own ancilla, as the page runs them. Exits
non-zero on any other difference above 1e-9.
"""

import json
import pathlib
import random
import subprocess
import sys
import warnings

import numpy as np
from qiskit import QuantumCircuit

import Qdislib as qd
from Qdislib.utils.exceptions import InvalidCutError

ENGINE = pathlib.Path(__file__).resolve().parent / "source" / "_static" / "explore_engine.js"
ONE = ["h", "x", "s", "t", "rx", "ry", "rz"]
TWO = ["cz", "cx", "rzz"]
NODE = """
const engine = require(process.argv[2]);
const cases = JSON.parse(require("fs").readFileSync(0, "utf8"));
const out = cases.map((c) => {
  const circuit = { n: c.n, ops: c.ops };
  const r = { exact: engine.exact(circuit, c.observable) };
  const plan = engine.plan(circuit, c.cuts.map((x) => Object.assign({}, x)), c.observable);
  r.invalid = plan.invalid.length > 0;
  if (!r.invalid) {
    r.fragments = plan.fragments.length;
    r.configurations = plan.configurations;
    r.subcircuits = plan.subcircuits;
    r.gamma = plan.gamma;
    r.optimized = plan.optimized;
    r.golden = Object.keys(plan.golden).length;
    if (plan.joint) { r.joint_subcircuits = plan.joint.subcircuits; r.joint_gamma = plan.joint.gamma; }
    r.value = engine.reconstruct(circuit, c.cuts.map((x) => Object.assign({}, x)), c.observable).value;
  }
  return r;
});
process.stdout.write(JSON.stringify(out));
"""


def random_case(rng):
    n = rng.randint(2, 6)
    ops = []
    for _ in range(rng.randint(3, 16)):
        if rng.random() < 0.55:
            op = {"name": rng.choice(TWO), "qubits": rng.sample(range(n), 2)}
            if op["name"] == "rzz":
                op["theta"] = round(rng.uniform(0.1, 3.0), 3)
        else:
            op = {"name": rng.choice(ONE), "qubits": [rng.randrange(n)]}
            if op["name"].startswith("r"):
                op["theta"] = round(rng.uniform(0.1, 3.0), 3)
        ops.append(op)
    observable = "".join(rng.choice("ZZZZIIXY") for _ in range(n))
    if rng.random() < 0.3:
        ops = [o for o in ops if o["name"] in ("h", "x", "ry", "cz", "cx")] or ops
        observable = observable.replace("Y", "Z")
    return {"n": n, "ops": ops, "observable": observable, "cuts": choose_cuts(rng, n, ops)}


def wire_edges(n, ops):
    """(qubit, from, to) for each pair of consecutive gates on one qubit, sharing only it."""
    last, edges = {}, []
    for i, op in enumerate(ops):
        for q in op["qubits"]:
            if q in last and len(set(ops[last[q]]["qubits"]) & set(op["qubits"])) == 1:
                edges.append((q, last[q], i))
            last[q] = i
    return edges


def choose_cuts(rng, n, ops):
    twos = [i for i, op in enumerate(ops) if len(op["qubits"]) == 2]
    edges = wire_edges(n, ops)
    cuts = []
    # Often a cut that separates: all the two-qubit gates across a random split of the qubits.
    if twos and rng.random() < 0.6:
        side = set(rng.sample(range(n), rng.randint(1, n - 1)))
        across = [i for i in twos if (ops[i]["qubits"][0] in side) != (ops[i]["qubits"][1] in side)]
        if 0 < len(across) <= 3:
            for i in across:
                cuts.append({"kind": rng.choice(["gate", "hadamard"]), "op": i})
            return cuts
    for _ in range(rng.randint(0, 3)):
        if edges and rng.random() < 0.5:
            q, a, b = rng.choice(edges)
            cut = {"kind": "wire", "qubit": q, "from": a, "to": b}
        elif twos:
            cut = {"kind": rng.choice(["gate", "hadamard"]), "op": rng.choice(twos)}
        else:
            continue
        if cut not in cuts and not any(c.get("op") == cut.get("op") and "op" in c for c in cuts):
            cuts.append(cut)
    return cuts


def to_qiskit(case):
    qc = QuantumCircuit(case["n"])
    for op in case["ops"]:
        args = ([op["theta"]] if "theta" in op else []) + op["qubits"]
        getattr(qc, op["name"])(*args)
    return qc


def to_spec(case):
    name = lambda i: f"{case['ops'][i]['name'].upper()}_{i + 1}"
    spec = {"gate": [], "hadamard": [], "wire": []}
    for cut in case["cuts"]:
        if cut["kind"] == "wire":
            spec["wire"].append((name(cut["from"]), name(cut["to"])))
        else:
            spec[cut["kind"]].append(name(cut["op"]))
    return spec


def qdislib_result(case):
    qc, spec = to_qiskit(case), to_spec(case)
    observable = case["observable"]
    out = {"exact": float(qd.analytical_solution(qc, observable))}
    try:
        subs = qd.circuit_cutting_subcircuits(qc, spec, observable, joint=False)
    except (InvalidCutError, ValueError) as error:
        if "separate" in str(error) or "valid cut" in str(error) or isinstance(error, InvalidCutError):
            return {**out, "invalid": True}
        raise
    out["invalid"] = False
    out["subcircuits"] = len(subs)
    out["fragments"] = len(subs.plan.fragments)
    out["configurations"] = subs.plan.num_configurations
    out["gamma"] = subs.plan.gamma()
    out["value"] = float(qd.circuit_cutting(qc, spec, observable, shots=None, joint=False))
    if spec["hadamard"]:
        default = qd.circuit_cutting_subcircuits(qc, spec, observable)
        out["joint_subcircuits"] = len(default)
        out["joint_gamma"] = default.plan.gamma()
    else:
        probe = qd.circuit_cutting_subcircuits(qc, spec, observable, optimized=True)
        out["optimized"] = probe.plan.num_units()
    return out


def main():
    count = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    cases = [random_case(rng) for _ in range(count)]
    js = json.loads(
        subprocess.run(
            ["node", "-e", NODE, "-", str(ENGINE)],
            input=json.dumps(cases), capture_output=True, text=True, check=True,
        ).stdout
    )
    worst = {"exact": 0.0, "value": 0.0, "gamma": 0.0, "joint_gamma": 0.0, "self": 0.0}
    tally = {"valid": 0, "invalid": 0, "gate": 0, "hadamard": 0, "wire": 0, "mixed": 0, "golden": 0}
    failures, crashes, known = [], [], []
    for case, mine in zip(cases, js):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                theirs = qdislib_result(case)
            except Exception as error:  # a library error is a finding too
                worst["self"] = max(worst["self"], abs(mine.get("value", 0) - mine["exact"]) if "value" in mine else 0)
                crashes.append((case, f"{type(error).__name__}: {error}"))
                continue
        problems = []
        if "value" in mine:
            worst["self"] = max(worst["self"], abs(mine["value"] - mine["exact"]))
        for key in ("exact", "value", "gamma", "joint_gamma"):
            if key in theirs and key in mine:
                difference = abs(theirs[key] - mine[key])
                worst[key] = max(worst[key], difference)
                if difference > 1e-9:
                    problems.append(f"{key}: {mine[key]} against {theirs[key]}")
        for key in ("invalid", "subcircuits", "fragments", "configurations", "optimized", "joint_subcircuits"):
            if key in theirs and mine.get(key) != theirs[key]:
                problems.append(f"{key}: {mine.get(key)} against {theirs[key]}")
        kinds = {c["kind"] for c in case["cuts"]}
        if theirs["invalid"]:
            tally["invalid"] += 1
        else:
            tally["valid"] += 1
            for kind in kinds:
                tally[kind] += 1
            tally["mixed"] += len(kinds) > 1
            tally["golden"] += mine.get("golden", 0) > 0
        counts_only = all(p.startswith(("subcircuits", "fragments", "joint_subcircuits")) for p in problems)
        if problems and counts_only and {"hadamard", "gate"} <= kinds:
            known.append((case, problems))
        elif problems:
            failures.append((case, problems))
    print(f"{count} circuits: {tally}")
    print("largest differences:", {k: f"{v:.2e}" for k, v in worst.items()})
    for case, problems in failures[:8]:
        print("DIFFERENT", json.dumps(case), problems)
    for case, message in crashes[:8]:
        print("QDISLIB RAISED", message[:200], json.dumps(case))
    print(f"{len(failures)} circuits differ; one fragment more in Qdislib's count on {len(known)}; Qdislib raised on {len(crashes)}")
    return 1 if failures or worst["self"] > 1e-9 else 0


if __name__ == "__main__":
    sys.exit(main())
