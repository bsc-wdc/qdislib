#!/usr/bin/python
#
#  Copyright 2002-2024 Barcelona Supercomputing Center (www.bsc.es)
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

import networkx as nx
import matplotlib.pyplot as plt

from pycompss.api.task import task
from qibo import models

class DAGgraph:
    """Direct Aciclyc Graph class.

    Representation of a direct aciclyc graph.
    """

    def __init__(self):
        self.nodes = []
        self.edges = []
        self.edges2 = []

    def add_node(self, gate):
        self.nodes.append(gate)

    def add_edge(self, gate1, gate2):
        self.edges.append((gate1, gate2))

    def add_edge2(self, gate1, gate2):
        self.edges2.append((gate1, gate2))

    def print_nodes(self):
        print("Nodes :", self.nodes)

    def print_edges(self):
        print("Edges: ", self.edges)


def build_dag(circuit, dag):
    for gate in circuit.queue:
        dag.add_node(gate)
        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
            for node in reversed(
                dag.nodes[:-1]
            ):  # Skip the last node since it is the current gate being added
                if (
                    qubit in node.qubits
                ):  # Check if the qubit is in the node's qubits
                    dag.add_edge(node, gate)
                    break
    for gate in circuit.queue:
        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
            for node in reversed(
                dag.nodes[:-1]
            ):  # Skip the last node since it is the current gate being added
                if (
                    qubit in node.qubits
                ):  # Check if the qubit is in the node's qubits
                    if (node, gate) in dag.edges or (gate, node) in dag.edges:
                        pass
                    else:
                        if node != gate:
                            dag.add_edge2(node, gate)
    return dag


def print_graph(graph):
    pos = nx.spring_layout(graph)  # Define layout for the nodes

    # Draw edges for the first group with blue color
    edges_first_group = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "blue"
    ]
    nx.draw_networkx_edges(
        graph,
        pos,
        edgelist=edges_first_group,
        edge_color="blue",
        width=2.0,
        alpha=0.7,
    )

    edges_second_group = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "red"
    ]
    nx.draw_networkx_edges(
        graph,
        pos,
        edgelist=edges_second_group,
        edge_color="red",
        width=2.0,
        alpha=0.7,
        style="dotted",
    )

    # Draw nodes and labels
    nx.draw_networkx_nodes(graph, pos, node_color="skyblue", node_size=500)
    nx.draw_networkx_labels(graph, pos, font_weight="bold", font_size=12)

    # plt.title(
    #    "Directed Graph with Two Edge Groups (Red edges in dotted line)"
    # )
    plt.show()


def create_graph(dag, digraph):
    new_nodes = [index + 1 for index, i in enumerate(dag.nodes)]
    labels = {
        i: index + 1 for index, i in enumerate(dag.nodes)
    }  # Start numbering from 0
    new_edges = [(labels[gate1], labels[gate2]) for gate1, gate2 in dag.edges]
    new_edges2 = [
        (labels[gate1], labels[gate2]) for gate1, gate2 in dag.edges2
    ]

    digraph.add_nodes_from(new_nodes)
    digraph.add_edges_from(new_edges, color="blue")
    digraph.add_edges_from(new_edges2, color="red")
    return digraph


def del_empty_qubits(circuit):
    empty_qubits = []
    for gate in circuit.queue:
        for i in gate.qubits:
            if i not in empty_qubits:
                empty_qubits.append(i)
        gate.qubits
    return empty_qubits

def print_graph(graph):
    pos = nx.spring_layout(graph)  # Define layout for the nodes

    # Draw edges for the first group with blue color
    edges_first_group = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "blue"
    ]
    nx.draw_networkx_edges(
        graph,
        pos,
        edgelist=edges_first_group,
        edge_color="blue",
        width=2.0,
        alpha=0.7,
    )

    edges_second_group = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "red"
    ]
    nx.draw_networkx_edges(
        graph,
        pos,
        edgelist=edges_second_group,
        edge_color="red",
        width=2.0,
        alpha=0.7,
        style="dotted",
    )

    # Draw nodes and labels
    nx.draw_networkx_nodes(graph, pos, node_color="skyblue", node_size=500)
    nx.draw_networkx_labels(graph, pos, font_weight="bold", font_size=12)

    # plt.title(
    #    "Directed Graph with Two Edge Groups (Red edges in dotted line)"
    # )
    plt.show()

@task(returns=list)
def gen_graph_circuit(new_circuit, observable_dict=None, verbose=False):
    list_subcircuits = []
    # convert to DAG and DIGRPAPH
    digraph = nx.Graph()
    dag = DAGgraph()

    build_dag(new_circuit, dag)
    create_graph(dag, digraph)

    subgraphs = list(nx.connected_components(digraph))
    if verbose:
        print(subgraphs)
    diff_list = []
    for subgraph in subgraphs:
        subgraph = sorted(subgraph)
        selected_elements = [dag.nodes[i - 1] for i in subgraph]
        # circuit_copy = copy.deepcopy(new_circuit)

        # remove specific qubit

        circuit_copy = models.Circuit(new_circuit.nqubits)
        circuit_copy.add(selected_elements)

        non_empty_qubits = del_empty_qubits(circuit_copy)
        non_empty_qubits.sort()
        # print(non_empty_qubits)
        difference_list = [
            value - index for index, value in enumerate(non_empty_qubits)
        ]
        # print("Non empty qubit ",difference_list)
        subtracted_list = [
            x - y for x, y in zip(non_empty_qubits, difference_list)
        ]
        if verbose:
            print("Substracted list: ", subtracted_list)
        diff_list.append(non_empty_qubits)
        for gate in circuit_copy.queue:
            if len(gate.qubits) > 1:
                control = subtracted_list[
                    non_empty_qubits.index(gate.qubits[0])
                ]
                gate._set_control_qubits((control,))
                target = subtracted_list[
                    non_empty_qubits.index(gate.qubits[1])
                ]
                gate._set_target_qubits((target,))
            else:
                target = subtracted_list[
                    non_empty_qubits.index(gate.qubits[0])
                ]
                gate._set_target_qubits((target,))
        circuit_copy.nqubits = len(non_empty_qubits)
        circuit_copy.queue.nmeasurements = 0
        list_subcircuits.append(circuit_copy)

    if verbose:
        print(non_empty_qubits)
    if verbose:
        print(diff_list)
    if observable_dict is not None:
        list_obs = []
        for p in diff_list:
            new_obs = {}
            for index, x in enumerate(p):
                if verbose:
                    print(observable_dict)

                new_obs[index] = observable_dict[x]
                if verbose:
                    print(new_obs)
            list_obs.append(new_obs)
        if verbose:
            print(list_obs)
        list_subcircuits_obs = [list_subcircuits, list_obs]
    else:
        list_subcircuits_obs = list_subcircuits

    return list_subcircuits_obs