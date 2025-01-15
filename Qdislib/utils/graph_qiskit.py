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
import networkx

from qiskit import QuantumCircuit
from qiskit import QuantumRegister
from qiskit import ClassicalRegister


def circuit_qiskit_to_dag(circuit: QuantumCircuit) -> networkx.DiGraph:
    """Convert a Qibo circuit to a DAG where each node stores gate information.

    :param circuit: The Qibo circuit to transform.
    :param num_qubits: The number of qubits in the circuit.
    :return: A directed acyclic graph (DAG) with nodes containing gate information.
    """
    # Create a directed graph
    dag = networkx.DiGraph()

    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.data, start=1):

        if gate.operation.name != "barrier":

            # Unique identifier for each gate instance
            gate_name = f"{gate.operation.name}_{gate_idx}"

            qubits = ()
            for elem in gate.qubits:
                qubits = qubits + (elem._index,)
                print(qubits)

            # Add the gate to the DAG, including the gate type, qubits, and parameters
            dag.add_node(
                gate_name,
                gate=gate.operation.name,
                qubits=qubits,
                parameters=gate.operation.params,
            )

            # Connect gates based on qubit dependencies
            for qubit in qubits:
                # Skip the last node since it is the current gate being added
                for pred_gate in reversed(list(dag.nodes)):
                    # Check if the qubit is in the node's qubits
                    if (
                        dag.nodes[pred_gate].get("qubits")
                        and qubit in dag.nodes[pred_gate]["qubits"]
                        and gate_name != pred_gate
                    ):
                        dag.add_edge(pred_gate, gate_name, color="blue")
                        break

            for qubit in qubits:
                for pred_gate in reversed(list(dag.nodes)):
                    if (
                        dag.nodes[pred_gate].get("qubits")
                        and qubit in dag.nodes[pred_gate]["qubits"]
                        and gate_name != pred_gate
                        and not dag.has_edge(pred_gate, gate_name)
                    ):
                        dag.add_edge(pred_gate, gate_name, color="red")

    return dag


def dag_to_circuit_qiskit(dag: networkx.DiGraph, num_qubits: int) -> QuantumCircuit:
    """Reconstruct a Qibo circuit from a DAG.

    :param dag: A networkx DiGraph representing the circuit.
    :param num_qubits: The number of qubits in the original circuit.
    :return: A Qibo circuit reconstructed from the DAG.
    """
    # Create an empty Qibo circuit
    qreg_q = QuantumRegister(num_qubits, "q")
    creg_c = ClassicalRegister(num_qubits, "c")
    circuit = QuantumCircuit(qreg_q, creg_c)

    # Traverse the DAG in topological order
    topo_order = list(networkx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data["gate"]

        # Skip the measurement nodes (we'll handle them separately)
        if gate_name == "Observable I":
            continue

        # Get the qubits this gate acts on

        qubits = node_data["qubits"]
        parameters = node_data["parameters"]

        # Get the gate class from the qibo.gates module
        gate_class = gate_name.lower()
        print(gate_class)

        # Get the signature of the gate's __init__ method
        signature = inspect.signature(gate_class.__init__)

        # Count the number of required positional arguments (excluding 'self')
        param_count = len(signature.parameters) - 1  # exclude 'self'

        print(signature)
        print(param_count)
        print(parameters)

        # Check if parameters are provided and the gate requires them
        if parameters:
            # Pass qubits and parameters if the gate requires both
            # circuit.gate_class(parameters, *qubits)
            tmp = getattr(circuit, gate_class)
            tmp(*parameters, *qubits)
        else:
            if gate_class == "measure":
                clbit = qubits
                tmp = getattr(circuit, gate_class)
                tmp(*qubits, *clbit)

            else:
                # Otherwise, pass only the qubits
                # circuit.gate_class(*qubits)
                tmp = getattr(circuit, gate_class)
                tmp(*qubits)

    # TODO: check the obs_I initialization since it will never go into the next if
    # Optionally handle measurements, assuming all qubits are measured at the end
    obs_I = []
    # for node in topo_order:
    #     node_data = dag.nodes[node]
    #     if node_data['gate'] == "Observable I":
    #         obs_I.append(node_data['qubits'][0])
    #         dag.remove_node(node)
    #         #print(obs_I)
    #         #circuit.add(gates.M(node_data['qubit']))

    if obs_I:
        return [circuit, obs_I]

    return [circuit, None]
