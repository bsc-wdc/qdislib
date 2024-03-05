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

import copy
import networkx as nx
import matplotlib.pyplot as plt


class DAGgraph:
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

    # print("Numeric nodes ",new_nodes)
    # print("Numeric edges ",digraph.edges(data=True))

    articulation_points = list(nx.articulation_points(digraph))
    # print("Articulation points ", articulation_points)

    # print_graph(digraph)
    return digraph


def del_empty_qubits(circuit):
    empty_qubits = []
    for gate in circuit.queue:
        for i in gate.qubits:
            if i not in empty_qubits:
                empty_qubits.append(i)
        gate.qubits
    return empty_qubits
