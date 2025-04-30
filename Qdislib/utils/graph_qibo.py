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
Qdislib graph utils.

This file contains all auxiliary graph classes and functions.
"""

import inspect
import networkx
import matplotlib.pyplot as plt
import typing
from qibo import models
from qibo import gates

from pycompss.api.task import task
from pycompss.api.parameter import INOUT
from pycompss.api.constraint import *
import numpy as np


def circuit_qibo_to_dag(circuit: models.Circuit, obs_I=None) -> networkx.DiGraph:
    """Convert a Qibo circuit into a directed acyclic graph (DAG) where each node represents a gate.

    This function constructs a DAG from a Qibo circuit. Each node in the graph corresponds to a gate in the circuit,
    storing information such as gate type, applied qubits, and parameters. Edges represent qubit-based dependencies 
    between gates, and optionally, nodes representing identity observables ("I") can be appended to sink nodes in 
    the graph.

    :param circuit: Qibo :class:`Circuit` object to convert into a DAG.
    :param obs_I: Optional iterable indicating which qubits have identity ("I") observables appended. Used for internal workflow for cutting. Default None.

    :return: A :class:`networkx.DiGraph` object where nodes represent gates and edges indicate gate dependencies.

    Node Attributes:
        - ``gate``: Name of the gate (e.g., "H", "CNOT").
        - ``qubits``: Tuple of qubit indices the gate acts on.
        - ``parameters``: Tuple of parameters for the gate (e.g., rotation angles).

    Edge Attributes:
        - ``color``: Connection type (e.g., "blue" or "red") used to distinguish dependency types.

    Example:
        .. code-block:: python

            from qibo import models, gates
            import networkx as nx
            from Qdislib.core.utils.graph_qibo import circuit_qibo_to_dag

            # Define a Qibo circuit
            c = models.Circuit(2)
            c.add(gates.H(0))
            c.add(gates.CNOT(0, 1))

            # Convert it to a DAG
            dag = circuit_qibo_to_dag(c)

            # Visualize or analyze the DAG
            print(dag.nodes(data=True))
            print(dag.edges(data=True))
    """
    # Create a directed graph
    dag = networkx.DiGraph()

    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.queue, start=1):
        # Unique identifier for each gate instance
        gate_name = f"{gate.__class__.__name__}_{gate_idx}".upper()

        # Add the gate to the DAG, including the gate type, qubits, and parameters
        dag.add_node(
            gate_name,
            gate=gate.__class__.__name__,
            qubits=gate.qubits,
            parameters=gate.parameters,
        )

        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
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

        for qubit in gate.qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and gate_name != pred_gate
                    and not dag.has_edge(pred_gate, gate_name)
                ):
                    dag.add_edge(pred_gate, gate_name, color="red")

    if obs_I:
        sinks = [node for node in dag.nodes if dag.out_degree(node) == 0]
        #print(sinks)
        for idx, i in enumerate(obs_I):
            if i == "I":
                #print(obs_I)
                for elem in sinks:
                    if dag.nodes[elem].get("qubits") == (idx,):
                        #print(dag.nodes[elem])
                        #print(dag.nodes[elem].get("qubits"))
                        dag.add_node(f"OBSI_{idx}", gate="Observable I", qubits=(idx,),parameters=())
                        dag.add_edge(elem, f"OBSI_{idx}", color="blue")

    return dag


def plot_dag(dag: networkx.DiGraph) -> None:
    """
    Visualize a directed acyclic graph (DAG) representation of a quantum circuit.

    This function uses `matplotlib` and `networkx` to plot a DAG where each node represents 
    a quantum gate and edges indicate qubit dependencies between gates. Edges can be colored 
    differently based on the "color" attribute in the graph:
    
    - **Blue**: Standard gate dependencies.
    - **Red (dotted)**: Additional or alternative dependencies (e.g., same qubit paths).

    :param dag: A :class:`networkx.DiGraph` object representing a quantum circuit DAG.

    :return: None

    The plot includes:
        - Nodes representing gates, labeled with gate names.
        - Edges with color-coded dependencies.
        - A spring layout for visually appealing positioning.

    Example:
        .. code-block:: python

            import networkx as nx
            from Qdislib.core.utils.graph_qibo import circuit_qibo_to_dag, plot_dag

            # Define a Qibo circuit
            c = models.Circuit(2)
            c.add(gates.H(0))
            c.add(gates.CNOT(0, 1))

            # Convert it to a DAG
            dag = circuit_qibo_to_dag(c)

            plot_dag(dag)
    """
    # Set up graph layout
    pos = networkx.spring_layout(dag)

    # Draw edges for the first group with blue color
    edges_first_group = [
        (edge[0], edge[1]) for edge in dag.edges.data("color") if edge[2] == "blue"
    ]
    networkx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_first_group,
        edge_color="blue",
        width=2.0,
        alpha=0.7,
    )

    edges_second_group = [
        (edge[0], edge[1]) for edge in dag.edges.data("color") if edge[2] == "red"
    ]
    networkx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_second_group,
        edge_color="red",
        width=2.0,
        alpha=0.7,
        style="dotted",
    )

    # Draw nodes
    node_labels = networkx.get_node_attributes(dag, "gate")
    networkx.draw_networkx_nodes(dag, pos, node_color="skyblue", node_size=1000)
    networkx.draw_networkx_labels(dag, pos, labels=node_labels, font_size=10)

    # Show plot
    plt.title("Circuit DAG")
    plt.show()


@task(returns=1)
def dag_to_circuit_qibo(dag: networkx.DiGraph, num_qubits: int) -> models.Circuit:
    """
    Reconstruct a Qibo quantum circuit from a DAG representation.

    This function takes a directed acyclic graph (DAG) where each node represents a quantum gate,
    and reconstructs a Qibo circuit by adding gates in topological order. Optionally, it detects
    "Observable I" nodes and separates them from the circuit for measurement processing.

    :param dag: A :class:`networkx.DiGraph` representing the quantum circuit DAG. Each node must 
                contain `gate`, `qubits`, and `parameters` attributes.
    :param num_qubits: The number of qubits in the original circuit (used to initialize the circuit).

    :return: A list of two elements:
             - A :class:`qibo.models.Circuit` object reconstructed from the DAG.
             - A list of qubit indices where "Observable I" measurements are applied (or `None` if not present).

    Gate nodes must have the following format:
        - `gate`: Name of the gate (must match a class in `qibo.gates`).
        - `qubits`: A tuple of qubit indices the gate acts on.
        - `parameters`: A tuple of gate parameters (or `None` if not applicable).

    Notes:
        - Measurement gates (e.g., `"Observable I"`) are excluded from reconstruction but their qubit 
          indices are returned for postprocessing.
        - The reconstruction uses Python's `inspect` module to dynamically determine how to instantiate gates.

    Example:
        .. code-block:: python

            from mymodule import dag_to_circuit_qibo
            import networkx as nx

            dag = ...  # previously constructed DAG
            circuit, observables = dag_to_circuit_qibo(dag, num_qubits=5)
            print(circuit)
            print(observables)
    """
    # Traverse the DAG in topological order
    topo_order = list(networkx.topological_sort(dag))

    # Optionally handle measurements, assuming all qubits are measured at the end
    obs_I = []
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data["gate"] == "Observable I":
            #print(node)
            obs_I.append(node_data["qubits"][0])
            dag.remove_node(node)

    dag, highest_qubit, smalles_qubit = _update_qubits_serie(dag)

    # Create an empty Qibo circuit
    circuit = models.Circuit(highest_qubit)

    topo_order = list(networkx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data["gate"].upper()

        # Skip the measurement nodes (we'll handle them separately)
        if gate_name == "OBSERVABLE I":
            #print(gate_name)
            continue

        if gate_name == "MEASURE":
            #qubits = node_data["qubits"]
            #circuit.add(gates.M(*qubits, collapse=True))
            #circuit.add(gates.RX(*qubits, theta=np.pi * output.symbols[0] / 4))
            continue

        # Get the qubits this gate acts on
        qubits = node_data["qubits"]
        parameters = node_data["parameters"]

        # Get the gate class from the qibo.gates module
        gate_class = getattr(gates, gate_name)

        # Get the signature of the gate's __init__ method
        signature = inspect.signature(gate_class.__init__)

        # Count the number of required positional arguments (excluding 'self')
        param_count = len(signature.parameters) - 1  # exclude 'self'

        # Check if parameters are provided and the gate requires them
        if parameters is not None and param_count > len(qubits):
            # Pass qubits and parameters if the gate requires both
            circuit.add(gate_class(*qubits, parameters))
        else:
            # Otherwise, pass only the qubits
            circuit.add(gate_class(*qubits))



    if obs_I:
        return [circuit, obs_I]

    return [circuit, None]


def _max_qubit(graph: networkx.DiGraph) -> float:
    """Get the highest Qubit value.

    :param graph: Graph to explore.
    :return: The value of the highest Qubit.
    """
    # Initialize a variable to keep track of the highest Qubits value
    max_qubits = float("-inf")  # Start with the lowest possible number
    # # It is possible to store the node:
    # max_node = None  # Store the node with the highest qubit

    # Iterate over the nodes and check their 'Qubits' attribute
    for node, data in graph.nodes(data=True):
        qubits = data.get("qubits", 0)  # Default to 0 if 'Qubits' is not present
        for qubit in qubits:
            if qubit > max_qubits:
                max_qubits = qubit
                # # Keep track of the node with the highest Qubit value
                # max_node = node
    return max_qubits

#@constraint(processors=[{"processorType": "GPU", "computingUnits": "1"}])
@task(returns=2, s=INOUT)
def _update_qubits(
    s: typing.List[typing.Any],
) -> typing.Tuple[typing.List[typing.Any], float]:
    """Update qubits task.

    :param s: Graph to explore.
    :return: The updated graph and the highest qubit value.
    """
    if not s.nodes():
        return None, None
    my_set = set()
    for node, _ in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    for node, _ in s.nodes(data=True):
        new_tuple = ()
        for qubit in s.nodes[node]["qubits"]:
            len_missing = _count_missing_up_to(my_set, qubit)
            new_qubit = qubit - len_missing
            new_tuple = new_tuple + (new_qubit,)
        s.nodes[node]["qubits"] = new_tuple

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    return s, highest_qubit


def _remove_red_edges(graph: networkx.DiGraph) -> networkx.DiGraph:
    """Remove red edges from the given graph.

    :param graph: Graph to process.
    :return: Updated input graph without red edges.
    """
    copy_dag = graph.copy()
    red_edges = []

    for ed in copy_dag.edges:
        if copy_dag.get_edge_data(ed[0], ed[1])["color"] == "red":
            red_edges.append(ed)

    copy_dag.remove_edges_from(red_edges)
    return copy_dag


def _count_missing_up_to(nums: typing.Set[int], max_num: int) -> int:
    """Retrieve the amount of missing numbers in the nums set.

    # TODO: Maybe max_num is not necessary and can be taken from max(nums)

    :param nums: Set of numbers.
    :param max_num: Highest number in nums.
    :return: Amount of missing numbers.
    """
    # Create a set of all numbers from 0 to max_num
    full_set = set(range(max_num + 1))

    # Subtract the given set from the full set to get the missing numbers
    missing_numbers = full_set - nums

    # Return the count of missing numbers
    return len(missing_numbers)


def _update_qubits_serie(
    s: typing.List[typing.Any],
) -> typing.Tuple[typing.List[typing.Any], int, int]:
    """Update the given serie of qubits.

    :param s: Input qubit serie.
    :return: Updated qubit serie, highest qubit value and lowest qubit value.
    """
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    for node, data in s.nodes(data=True):
        new_tuple = ()
        for qubit in s.nodes[node]["qubits"]:
            len_missing = _count_missing_up_to(my_set, qubit)
            new_qubit = qubit - len_missing
            new_tuple = new_tuple + (new_qubit,)
        s.nodes[node]["qubits"] = new_tuple

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    smallest_qubit = min(my_set) - _count_missing_up_to(my_set, min(my_set))
    return s, highest_qubit, smallest_qubit


def _max_qubits_graph(
    s: typing.List[typing.Any],
) -> typing.Tuple[typing.List[typing.Any], int, int]:
    """Update the given serie of qubits.

    :param s: Input qubit serie.
    :return: Updated qubit serie, highest qubit value and lowest qubit value.
    """
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)

    highest_qubit = max(my_set) + 1 - _count_missing_up_to(my_set, max(my_set))
    smallest_qubit = min(my_set) - _count_missing_up_to(my_set, min(my_set))
    return highest_qubit

@task(returns=1)
def _dag_to_circuit_qibo_subcircuits(dag: networkx.DiGraph, num_qubits: int) -> models.Circuit:
    # Traverse the DAG in topological order
    topo_order = list(networkx.topological_sort(dag))

    # Optionally handle measurements, assuming all qubits are measured at the end
    obs_I = []
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data["gate"] == "Observable I":
            #print(node)
            obs_I.append(node_data["qubits"][0])
            dag.remove_node(node)

    dag, highest_qubit, smalles_qubit = _update_qubits_serie(dag)

    # Create an empty Qibo circuit
    circuit = models.Circuit(highest_qubit)

    topo_order = list(networkx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data["gate"].upper()

        # Skip the measurement nodes (we'll handle them separately)
        if gate_name == "OBSERVABLE I":
            #print(gate_name)
            continue

        if gate_name == "MEASURE":
            #qubits = node_data["qubits"]
            #circuit.add(gates.M(*qubits, collapse=True))
            #circuit.add(gates.RX(*qubits, theta=np.pi * output.symbols[0] / 4))
            continue

        # Get the qubits this gate acts on
        qubits = node_data["qubits"]
        parameters = node_data["parameters"]

        # Get the gate class from the qibo.gates module
        gate_class = getattr(gates, gate_name)

        # Get the signature of the gate's __init__ method
        signature = inspect.signature(gate_class.__init__)

        # Count the number of required positional arguments (excluding 'self')
        param_count = len(signature.parameters) - 1  # exclude 'self'

        # Check if parameters are provided and the gate requires them
        if parameters is not None and param_count > len(qubits):
            # Pass qubits and parameters if the gate requires both
            circuit.add(gate_class(*qubits, parameters))
        else:
            # Otherwise, pass only the qubits
            circuit.add(gate_class(*qubits))

    return circuit