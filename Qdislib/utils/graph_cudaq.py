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
Qdislib qiskit graph utils.

This file contains all auxiliary qiskit graph classes and functions.
"""

import inspect
import networkx as nx
import cudaq
import re
import json
from Qdislib.utils.graph_qibo import _update_qubits_serie

try:
    from pycompss.api.task import task

    pycompss_available = True

except ImportError:
    print("PyCOMPSs is NOT installed. Proceeding with serial execution (no parallelism).")

    def task(*args, **kwargs):
        def decorator(func):
            return func

        return decorator

    pycompss_available = False





def circuit_cudaq_to_dag(kernel, obs_I=None) -> nx.DiGraph:
    """
    Convert a CUDA-Q kernel into a DAG by parsing its Python source.
    """
    dag = nx.DiGraph()

    # Grab source code from kernel JSON
    kjson = kernel.to_json()
    if isinstance(kjson, str):
        kjson = json.loads(kjson)
    src = kjson["funcSrc"].splitlines()

    gate_idx = 0
    for line in src:
        line = line.strip()
        if not line or line.startswith("@") or line.startswith("def") or line.startswith("q =") or line.startswith("b0"):
            continue

        # Match gates: e.g., h(q[0]), cx(q[0], q[1]), rx(1.57, q[0])
        m = re.match(r"([a-zA-Z_][a-zA-Z0-9_]*)\((.*)\)", line)
        if not m:
            continue

        op_name = m.group(1).upper()
        args = m.group(2)

        # Extract qubit indices: q[0], q[1], etc.
        qubits = tuple(int(idx) for idx in re.findall(r"q\[(\d+)\]", args))

        # Remove q[...] parts, keep only potential numeric literals
        args_clean = re.sub(r"q\[\d+\]", "", args)
        params = [float(x) for x in re.findall(r"[-+]?\d*\.\d+|\d+", args_clean)]

        gate_idx += 1
        gate_name = f"{op_name}_{gate_idx}"

        dag.add_node(
            gate_name,
            gate=op_name,
            qubits=qubits,
            parameters=params
        )

         # --- Add edges based on qubit dependencies ---
        for qubit in qubits:
            # Blue edges: connect to the nearest preceding gate on the same qubit
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                ):
                    dag.add_edge(pred_gate, gate_name, color="blue")
                    break

        for qubit in qubits:
            # Red edges: connect to all earlier gates on the same qubit that aren't already connected
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and pred_gate != gate_name
                    and not dag.has_edge(pred_gate, gate_name)
                ):
                    dag.add_edge(pred_gate, gate_name, color="red")

    # Handle Observable I insertion
    if obs_I:
        sinks = [n for n in dag.nodes if dag.out_degree(n) == 0]
        for idx, i in enumerate(obs_I):
            if i == "I":
                for elem in sinks:
                    if dag.nodes[elem].get("qubits") == (idx,):
                        obs_node = f"OBSI_{idx}"
                        dag.add_node(
                            obs_node,
                            gate="Observable I",
                            qubits=(idx,),
                            parameters=[]
                        )
                        dag.add_edge(elem, obs_node, color="blue")

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

#@task(returns=1)
def dag_to_circuit_cudaq_kernel(dag: nx.DiGraph, num_qubits: int):
    """
    Convert a DAG into a CUDA-Q kernel using the modern builder API.

    Produces code like:
        kernel = cudaq.make_kernel()
        q = [kernel.qalloc() for _ in range(num_qubits)]
        kernel.h(q[0])
        result = kernel.mz(q[0])
        kernel.if_result(result, lambda: kernel.x(q[1]))

    Args:
        dag: networkx.DiGraph with node attributes:
            - 'gate': str (gate name)
            - 'qubits': tuple of qubit indices
            - 'parameters': list of parameters (if any)
        num_qubits: int, total number of qubits

    Returns:
        [kernel, obs_I]
    """
    if dag is None:
        return None

    topo_order = list(nx.topological_sort(dag))

    # Extract "Observable I" nodes if present
    """obs_I = []
    for node in topo_order:
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)"""

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

        # --- Qiskit-style MEASURE → conditional ---
        elif gate == "MEASURE":
            for qb in qubits:
                result = kernel.mz(q[qb])
                if result:
                    kernel.h(q[qb])

        else:
            raise ValueError(f"Unsupported gate type: {gate}")
    
    """if obs_I:
        observables = ["Z"] * num_qubits
        for element in obs_I:
            observables[element] = "I"

        observables = "".join(observables)
        return [kernel, observables]

    observables = ["Z"] * num_qubits
    observables = "".join(observables)"""
    return kernel



@task(returns=1)
def _dag_to_circuit_cudaq_subcircuits(dag: nx.DiGraph, num_qubits: int):

    return 0
