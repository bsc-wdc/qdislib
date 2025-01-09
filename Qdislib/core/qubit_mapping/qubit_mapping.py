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

"""Qbit mapping algorithms."""

# TODO: Complete file commented out. Remove?
# import networkx as nx
# import matplotlib.pyplot as plt
# from qibo import models, gates
# from Qdislib.classes.circuit_classes import _NewCircuit
#
#
# def architecture_x():
#     """
#     Description
#     -----------
#     Generate a graph representing the architecture of a quantum device with an 'X' shape topology. The graph consists of nodes representing qubits and edges representing connections between qubits.
#
#     Returns
#     -------
#     graph: NetworkX graph
#         A graph representing the 'X' shape topology architecture.
#
#     Example
#     -------
#     >>> from networkx import draw
#     >>> G = architecture_X()
#     >>> draw(G, with_labels=True)
#     """
#     x_graph = nx.Graph()
#     x_graph.add_nodes_from(["A", "B", "C", "D", "E"])
#     x_graph.add_edges_from([("A", "B"), ("B", "C"), ("B", "D"), ("B", "E")])
#     return x_graph
#
#
# def _qubit_arch(circuit, draw=False):
#     """
#     Extract qubit connections from a quantum circuit and generates
#     a graph representing the architecture based on these connections. Each edge
#     in the graph represents a connection between two qubits.
#
#     :param circuit: Circuit.
#     :param drwa: bool.
#     :return: graph
#     """
#     graph = nx.Graph()
#     lst = []
#     for gate in circuit.queue:
#         if len(gate.qubits) > 1:
#             lst.append(gate)
#             graph.add_edge(gate.qubits[0], gate.qubits[1])
#
#     pos = nx.spring_layout(graph)  # positions for all nodes
#     if draw:
#         nx.draw(
#             graph,
#             pos,
#             with_labels=True,
#             font_weight="bold",
#             node_size=700,
#             node_color="skyblue",
#             font_color="black",
#             font_size=10,
#         )
#         plt.show()
#     return graph
#
#
# def subgraph_matcher(architecture, circuit, draw=False):
#     """
#
#     Description
#     -----------
#     Check if a given circuit graph (pattern graph) is a subgraph of another graph representing the architecture (target graph). It uses NetworkX's GraphMatcher to perform subgraph isomorphism checks.
#
#     Parameters
#     ----------
#     architecture: NetworkX graph
#         The target graph representing the architecture.
#     circuit_graph: NetworkX graph
#         The pattern graph representing the circuit.
#     draw: bool, optional
#         If True, the function will visualize the target and pattern graphs. Default is False.
#
#     Returns
#     -------
#     all_subgraph_isomorphisms: list
#         A list of all subgraph isomorphisms found between the pattern graph and the target graph. Each isomorphism is represented as a dictionary mapping nodes from the pattern graph to nodes in the target graph.
#
#     Example
#     -------
#     >>> target_graph = nx.Graph()
#     >>> target_graph.add_edges_from([(1, 2), (1, 3), (2, 3)])
#     >>> pattern_graph = nx.Graph()
#     >>> pattern_graph.add_edges_from([(1, 2), (2, 3)])
#     >>> all_subgraph_isomorphisms = subgraph_matcher(target_graph, pattern_graph, draw=True)
#     """
#     # Create an example target graph
#     circuit_graph = _qubit_arch(circuit)
#
#     target_graph = architecture
#
#     pos = nx.spring_layout(target_graph)  # positions for all nodes
#     if draw:
#         nx.draw(
#             target_graph,
#             pos,
#             with_labels=True,
#             font_weight="bold",
#             node_size=700,
#             node_color="skyblue",
#             font_color="black",
#             font_size=10,
#         )
#         plt.show()
#
#     # Create an example pattern graph
#     pattern_graph = circuit_graph
#
#     pos = nx.spring_layout(pattern_graph)  # positions for all nodes
#     if draw:
#         nx.draw(
#             pattern_graph,
#             pos,
#             with_labels=True,
#             font_weight="bold",
#             node_size=700,
#             node_color="skyblue",
#             font_color="black",
#             font_size=10,
#         )
#         plt.show()
#
#     # Initialize GraphMatcher with the target and pattern graphs
#     matcher = nx.algorithms.isomorphism.GraphMatcher(
#         target_graph, pattern_graph
#     )
#
#     # Check if the pattern graph is a subgraph of the target graph
#     is_subgraph = matcher.subgraph_is_isomorphic()
#
#     # If there is a subgraph isomorphism, print the mapping
#     if is_subgraph:
#         mapping = matcher.mapping
#         print("Subgraph found! Node mapping:", mapping)
#     else:
#         print("No subgraph isomorphism found.")
#
#     # You can also get all subgraph isomorphisms
#     all_subgraph_isomorphisms = list(matcher.subgraph_isomorphisms_iter())
#     print("All subgraph isomorphisms:", all_subgraph_isomorphisms)
#     return all_subgraph_isomorphisms
#
#
# def mapping_order(target, pattern, order):
#     """
#     Description
#     -----------
#     Find the best mapping between nodes in a pattern graph and nodes in a target graph based on a desired order of nodes. It uses NetworkX's GraphMatcher to perform subgraph isomorphism checks and calculates the best mapping based on the specified order.
#
#     Parameters
#     ----------
#     target: NetworkX graph
#         The target graph.
#     pattern: NetworkX graph
#         The pattern graph.
#     order: list
#         A list specifying the desired order of nodes. Nodes present in the order list will be prioritized in the mapping process.
#
#     Returns
#     -------
#     best_architecture: dict
#         The best mapping found between nodes in the pattern graph and nodes in the target graph, based on the desired order.
#
#     Example
#     -------
#     >>> target_graph = nx.Graph()
#     >>> target_graph.add_edges_from([(1, 2), (1, 3), (2, 3)])
#     >>> pattern_graph = nx.Graph()
#     >>> pattern_graph.add_edges_from([(1, 2), (2, 3)])
#     >>> order = [3, 1, 2]
#     >>> best_architecture = mapping_order(target_graph, pattern_graph, order)
#     """
#
#     # Create an example target graph
#     target_graph = target
#
#     # Create an example pattern graph
#     pattern_graph = pattern
#
#     # Specify the desired order of nodes
#     desired_order = order
#
#     # Initialize GraphMatcher with the target and pattern graphs
#     matcher = nx.algorithms.isomorphism.GraphMatcher(
#         target_graph, pattern_graph
#     )
#
#     # Get all subgraph isomorphisms
#     all_subgraph_isomorphisms = list(matcher.subgraph_isomorphisms_iter())
#
#     results = []
#     for element in all_subgraph_isomorphisms:
#         result = 0
#         for index, x in enumerate(desired_order):
#             if x in element:
#                 result = result + (index + 1)
#         results.append([element, result])
#     print(results)
#     min_second_element = min(results, key=lambda x: x[1])
#     best_architecture = min_second_element[0]
#     print("Best architecture mapping: ", best_architecture)
#     return best_architecture
#
#
# def rename_qubits(subcirc, qubit_middle, best_arch, middle_arch_qubit):
#     """
#     Description
#     -----------
#     Rename the qubits in a given circuit according to a specified mapping provided by the ``best_arch`` dictionary. It updates the qubit labels based on the provided mapping, where the ``middle_arch_qubit`` represents the qubit to be renamed, and the ``qubit_middle`` represents its new label.
#
#     Parameters
#     ----------
#     subcirc: Circuit.
#         The circuit whose qubits need to be renamed.
#     qubit_middle: int
#         The new label for the qubit to be renamed.
#     best_arch: dict
#         A dictionary representing the best mapping between qubits in the original circuit and the target circuit.
#     middle_arch_qubit: int
#         The qubit in the target architecture that corresponds to the qubit to be renamed.
#
#     Returns
#     -------
#     new_circuit: Circuit.
#         The modified circuit with renamed qubits.
#
#     Example
#     -------
#     >>> subcirc = models.Circuit(3)
#     >>> subcirc.add(models.H(0))
#     >>> subcirc.add(models.CNOT(0, 1))
#     >>> best_arch = {0: 2, 1: 1, 2: 0}
#     >>> new_circuit = rename_qubits(subcirc, 1, best_arch, 0)
#     """
#     if subcirc.nqubits <= 2:
#         new_subcirc = models.Circuit(subcirc.nqubits + 1)
#         new_subcirc.queue = subcirc.queue
#         subcirc = new_subcirc
#     print(subcirc.draw())
#     target_qubit = best_arch[middle_arch_qubit]
#     print(type(subcirc))
#     if isinstance(subcirc, _NewCircuit):
#         _NewCircuit_class = True
#         new_circuit = models.Circuit(subcirc.circuit.nqubits)
#         new_circuit.queue = subcirc.circuit.queue
#         subcirc = subcirc.circuit
#     else:
#         _NewCircuit_class = False
#         new_circuit = models.Circuit(subcirc.nqubits)
#         new_circuit.queue = subcirc.queue
#     for element in subcirc.queue:
#         qubit = element.qubits
#         if len(qubit) < 2:
#             qubit = qubit[0]
#             if qubit == target_qubit:
#                 element._set_target_qubits((qubit_middle,))
#             elif qubit == qubit_middle:
#                 element._set_target_qubits((target_qubit,))
#         else:
#             qubit_0 = element.qubits[0]
#             qubit_1 = element.qubits[1]
#             if qubit_0 == target_qubit:
#                 element._set_control_qubits((qubit_middle,))
#             elif qubit_1 == target_qubit:
#                 element._set_target_qubits((qubit_middle,))
#             if qubit_0 == qubit_middle:
#                 element._set_control_qubits((target_qubit,))
#             elif qubit_1 == qubit_middle:
#                 element._set_target_qubits((target_qubit,))
#     if _NewCircuit_class:
#         new_circuit = _NewCircuit(new_circuit)
#     return new_circuit
