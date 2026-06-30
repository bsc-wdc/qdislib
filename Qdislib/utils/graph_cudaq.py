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

"""
Qdislib CUDA-Q graph utils.

This file contains all auxiliary CUDA-Q graph classes and functions.
"""

import logging

logger = logging.getLogger(__name__)

__all__ = [
    "circuit_cudaq_to_dag",
    "dag_to_circuit_cudaq",
    "dag_to_circuit_cudaq_kernel",
]

import ast
import json
import math
import networkx as nx
import numpy as np

try:
    import cudaq

    cudaq_available = True
except ImportError:
    cudaq = None
    cudaq_available = False

try:
    from pycompss.api.task import task

    pycompss_available = True

except ImportError:
    logger.warning("PyCOMPSs is not installed; running serially (no parallelism).")

    def task(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    pycompss_available = False


def _cudaq_gate_name(call):
    """Upper-cased gate name for a CUDA-Q gate call AST node, or ``None``.

    Handles plain calls (``h(...)`` -> ``"H"``) and the controlled-gate attribute
    syntax (``x.ctrl(c, t)`` -> ``"CX"``, ``z.ctrl(c, t)`` -> ``"CZ"``). Returns
    ``None`` for non-gate calls such as ``cudaq.qvector(...)``.
    """
    func = call.func
    if isinstance(func, ast.Name):
        return func.id.upper()
    if (
        isinstance(func, ast.Attribute)
        and func.attr == "ctrl"
        and isinstance(func.value, ast.Name)
    ):
        base = func.value.id.upper()
        return {"X": "CX", "Z": "CZ"}.get(base, "C" + base)
    return None


def _eval_cudaq_param(node):
    """Best-effort evaluation of a numeric gate-parameter expression node.

    Supports literals and expressions such as ``np.pi / 2`` / ``math.pi``;
    returns ``None`` if the value cannot be resolved to a float.
    """
    try:
        value = eval(  # noqa: S307 - trusted kernel source, restricted namespace
            compile(ast.Expression(node), "<cudaq-param>", "eval"),
            {"np": np, "numpy": np, "math": math, "pi": math.pi},
        )
        return float(value)
    except Exception:
        return None


def circuit_cudaq_to_dag(kernel, obs_I=None) -> nx.DiGraph:
    """Convert a CUDA-Q kernel into a Qdislib DAG.

    The kernel's Python source (``kernel.to_json()["funcSrc"]``) is parsed with
    the :mod:`ast` module rather than scanned with regular expressions, so it is
    robust to formatting/comments, the qubit-register variable name, expression
    parameters (e.g. ``np.pi / 2``) and the ``gate.ctrl(...)`` controlled syntax.

    :param kernel: A CUDA-Q kernel (``PyKernelDecorator``).
    :param obs_I: Optional observable string; each ``"I"`` appends an identity
        observable sink node on that qubit (internal use when cutting).
    :return: A :class:`networkx.DiGraph` whose nodes carry ``gate`` (upper-cased
        name), ``qubits`` (tuple of ints) and ``parameters`` (list of floats),
        with blue (nearest-predecessor) and red (other-predecessor) edges.
    """
    dag = nx.DiGraph()

    payload = kernel.to_json()
    if isinstance(payload, str):
        payload = json.loads(payload)
    tree = ast.parse(payload["funcSrc"])
    func = next(
        (node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)), None
    )
    if func is None:
        return dag

    gate_idx = 0
    for statement in func.body:
        # Gate operations are bare calls (``h(q[0])``) or measurement
        # assignments (``m = mz(q[0])``); everything else is skipped.
        if isinstance(statement, ast.Expr) and isinstance(statement.value, ast.Call):
            call = statement.value
        elif isinstance(statement, ast.Assign) and isinstance(
            statement.value, ast.Call
        ):
            call = statement.value
        else:
            continue

        name = _cudaq_gate_name(call)
        if name is None:
            continue

        # Subscript arguments (``q[i]``) are qubits; the rest are parameters.
        qubits, parameters = [], []
        for arg in call.args:
            if isinstance(arg, ast.Subscript):
                index = arg.slice.value if isinstance(arg.slice, ast.Constant) else None
                if isinstance(index, int):
                    qubits.append(index)
            else:
                value = _eval_cudaq_param(arg)
                if value is not None:
                    parameters.append(value)

        if not qubits:  # e.g. cudaq.qvector(...) or a whole-register op -> skip
            continue

        gate_idx += 1
        gate_name = f"{name}_{gate_idx}"
        dag.add_node(gate_name, gate=name, qubits=tuple(qubits), parameters=parameters)

        # Blue edges: nearest preceding gate on each shared qubit.
        for qubit in qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                ):
                    dag.add_edge(pred_gate, gate_name, color="blue")
                    break

        # Red edges: other earlier gates on a shared qubit, if not already linked.
        for qubit in qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                    and not dag.has_edge(pred_gate, gate_name)
                ):
                    dag.add_edge(pred_gate, gate_name, color="red")

    # Append identity-observable sink nodes when requested.
    if obs_I:
        sinks = [n for n in dag.nodes if dag.out_degree(n) == 0]
        for idx, observable in enumerate(obs_I):
            if observable == "I":
                for sink in sinks:
                    if dag.nodes[sink].get("qubits") == (idx,):
                        obs_node = f"OBSI_{idx}"
                        dag.add_node(
                            obs_node, gate="Observable I", qubits=(idx,), parameters=[]
                        )
                        dag.add_edge(sink, obs_node, color="blue")

    return dag


@task(returns=1)
def dag_to_circuit_cudaq(dag: nx.DiGraph, num_qubits: int):
    if dag is None:
        return None

    topo_order = list(nx.topological_sort(dag))

    # Extract "Observable I" nodes if present
    obs_I = []
    for node in topo_order:
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)

    if obs_I:
        observables = ["Z"] * num_qubits
        for element in obs_I:
            observables[element] = "I"

        observables = "".join(observables)
        return [(dag, num_qubits), observables]

    observables = ["Z"] * num_qubits
    observables = "".join(observables)
    return [(dag, num_qubits), observables]


# @task(returns=1)
def dag_to_circuit_cudaq_kernel(dag: nx.DiGraph, num_qubits: int):
    """
    Convert a DAG into a CUDA-Q kernel using the modern builder API.

    Produces code like::

        kernel = cudaq.make_kernel()
        q = [kernel.qalloc() for _ in range(num_qubits)]
        kernel.h(q[0])
        result = kernel.mz(q[0])
        kernel.if_result(result, lambda: kernel.x(q[1]))

    :param dag: networkx.DiGraph whose nodes carry ``gate`` (str), ``qubits``
        (tuple of indices) and ``parameters`` (list) attributes.
    :param num_qubits: Total number of qubits.
    :return: ``[kernel, observables]``.
    """
    if dag is None:
        return None

    # Initialize CUDA-Q kernel
    kernel = cudaq.make_kernel()

    # Allocate qubits explicitly
    q = [kernel.qalloc() for _ in range(num_qubits)]
    # Traverse DAG
    for node in nx.topological_sort(dag):
        node_data = dag.nodes[node]
        gate = node_data["gate"].upper()
        qubits = node_data["qubits"]
        params = node_data["parameters"]

        # --- Common 1-qubit gates ---
        if gate in ["H", "X", "Y", "Z", "T", "S", "SDG"]:
            for qb in qubits:
                getattr(kernel, gate.lower())(q[qb])

        # --- Rotation gates ---
        elif gate in ["RX", "RY", "RZ"]:
            angle = params[0] if params else 0.0
            for qb in qubits:
                getattr(kernel, gate.lower())(angle, q[qb])

        # --- Controlled gates ---
        elif gate in ["CNOT", "CX"]:
            ctrl, tgt = qubits
            kernel.cx(q[ctrl], q[tgt])

        elif gate == "CZ":
            ctrl, tgt = qubits
            kernel.cz(q[ctrl], q[tgt])

        # --- Measurement gates ---
        elif gate in ["MZ", "MX", "MY"]:
            for qb in qubits:
                getattr(kernel, gate.lower())(q[qb])

        # --- Gate-cut MEASURE: collapse the qubit, then force it to |1> so it
        # --- contributes -1, matching the Qiskit/Qibo gate-cut convention
        # --- (qiskit: measure + "if c==0: x"; qibo: M(collapse) + RX(pi*(1-m))).
        elif gate == "MEASURE":
            for qb in qubits:
                kernel.mz(q[qb])
                kernel.reset(q[qb])
                kernel.x(q[qb])

        else:
            raise ValueError(f"Unsupported gate type: {gate}")

    return kernel
