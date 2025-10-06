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




def dag_to_circuit_cudaq(dag: nx.DiGraph, num_qubits: int):
    """
    Reconstruct a CUDA-Q kernel from a DAG representation.

    Args:
        dag (networkx.DiGraph): DAG with nodes containing:
            - 'gate': str, name of the gate
            - 'qubits': tuple of qubit indices
            - 'parameters': list of parameters (if any)
        num_qubits (int): number of qubits in the circuit.

    Returns:
        kernel: a @cudaq.kernel function implementing the circuit.
    """
    if dag is None:
        return None

    topo_order = list(nx.topological_sort(dag))

    # Track observable I (to return separately)
    obs_I = []
    for node in topo_order:
        if dag.nodes[node]["gate"].upper() == "OBSERVABLE I":
            obs_I.append(dag.nodes[node]["qubits"][0])
            dag.remove_node(node)
    
    dag, highest_qubit, smalles_qubit = _update_qubits_serie(dag)

    # Define the CUDA-Q kernel dynamically
    @cudaq.kernel
    def kernel():
        q = cudaq.qvector(num_qubits)

        for node in topo_order:
            node_data = dag.nodes[node]
            gate = node_data["gate"].upper()
            qubits = node_data["qubits"]
            params = node_data["parameters"]

            if gate == "OBSERVABLE I":
                # skip (we only track it separately)
                continue
            if gate in ["CNOT", "CX"]:
                ctrl, tgt = qubits
                x.ctrl(q[ctrl], q[tgt])
            elif gate == "CZ":
                ctrl, tgt = qubits
                z.ctrl(q[ctrl], q[tgt])
            elif gate in ["MZ"]:
                mz([q[i] for i in qubits])
            elif gate == "MX":
                mx([q[i] for i in qubits])
            elif gate == "MY":
                my([q[i] for i in qubits])
            elif gate =="MEASURE":
                # Execute measurement
                result = mz([q[i] for i in qubits])
                # Conditional X as in Qiskit if_test
                if result:
                    x([q[i] for i in qubits])

            else:
                # generic gate: e.g. h, rx, ry, rz, x, z, etc.
                gate_fn = globals().get(gate.lower())
                if gate_fn is None:
                    raise ValueError(f"Unsupported gate: {gate}")

                if params:
                    gate_fn(*params, *[q[i] for i in qubits])
                else:
                    gate_fn(*[q[i] for i in qubits])

    return [kernel, obs_I if obs_I else None]



@task(returns=1)
def _dag_to_circuit_cudaq_subcircuits(dag: nx.DiGraph, num_qubits: int):
    """
    Reconstruct a CUDA-Q kernel from a DAG (subcircuit version).
    
    Args:
        dag (networkx.DiGraph): DAG with nodes containing:
            - 'gate': str, name of the gate
            - 'qubits': tuple of qubit indices
            - 'parameters': list of parameters (if any)
        num_qubits (int): Number of qubits in the circuit.
    
    Returns:
        kernel: A @cudaq.kernel function implementing the circuit.
    """
    if dag is None:
        return None

    # Topological order for dependency resolution
    topo_order = list(nx.topological_sort(dag))

    # Remove Observable I nodes (we don't encode them in CUDA-Q subcircuits)
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data["gate"] == "Observable I":
            dag.remove_node(node)

    dag, highest_qubit, smalles_qubit = _update_qubits_serie(dag)
    
    topo_order = list(nx.topological_sort(dag))

    @cudaq.kernel
    def kernel():
        q = cudaq.qvector(num_qubits)

        for node in topo_order:
            node_data = dag.nodes[node]
            gate_name = node_data["gate"].upper()
            qubits = node_data["qubits"]
            params = node_data["parameters"]

            # Skip Observable I
            if gate_name == "OBSERVABLE I":
                continue

            # Special case: CNOT
            if gate_name in ["CNOT", "CX"]:
                ctrl, tgt = qubits
                x.ctrl(q[ctrl], q[tgt])
                continue

            # Special case: CZ
            if gate_name == "CZ":
                ctrl, tgt = qubits
                z.ctrl(q[ctrl], q[tgt])
                continue

            # Measurement handling
            if gate_name == "MZ":
                mz([q[i] for i in qubits])
                continue
            elif gate_name == "MX":
                mx([q[i] for i in qubits])
                continue
            elif gate_name == "MY":
                my([q[i] for i in qubits])
                continue
            elif gate_name == "MEASURE":
                # Execute measurement
                result = mz([q[i] for i in qubits])
                # Conditional X as in Qiskit if_test
                if result:
                    x([q[i] for i in qubits])
                continue

            # Generic single- or multi-qubit gates (rx, ry, rz, h, x, y, z, etc.)
            gate_fn = globals().get(gate_name.lower())
            if gate_fn is None:
                raise ValueError(f"Unsupported gate: {gate_name}")

            if params:
                gate_fn(*params, *[q[i] for i in qubits])
            else:
                gate_fn(*[q[i] for i in qubits])

    return kernel
