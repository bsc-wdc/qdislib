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

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *

import networkx as nx
import matplotlib.pyplot as plt
import copy


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
    #     "Directed Graph with Two Edge Groups (Red edges in dotted line)"
    # )
    plt.show()


def gates_dict(circuit):
    double_gates = {}
    for index, gate in enumerate(circuit.queue):
        if len(gate.qubits) > 1:
            start, end = gate.qubits
            new_tuple = [(i, i + 1) for i in range(start, end)]
            for tupl in new_tuple:
                if tupl not in double_gates:
                    double_gates[tupl] = []
                    double_gates[tupl].append(index + 1)
                else:
                    double_gates[tupl].append(index + 1)
    return double_gates


def has_number(lst, number):
    for num in lst:
        if num == number:
            return True
    return False


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

    # articulation_points = list(nx.articulation_points(digraph))
    # print("Articulation points ", articulation_points)
    # print("CIRCUIT")
    # print_graph(digraph)
    return digraph


def del_empty_qubits(circuit):
    empty_qubits = []
    for gate in circuit.queue:
        for i in gate.qubits:
            if i not in empty_qubits:
                empty_qubits.append(i)
    return empty_qubits


@task(returns=list)
def gate_selector(
    digraph, circuit, max_qubits=None, num_subcirucits=None, draw=False
):
    # ----------------------------------------------------
    # LOOP THROUGH DOUBLE GATES POINTS
    # ----------------------------------------------------
    # num_nodes_in_subgraphs = []
    temp = {}
    right_subgrafs = []
    computational_cost = []
    double_gates = gates_dict(circuit)
    for key, array in double_gates.items():
        num_nodes = []

        copy_dag = copy.deepcopy(digraph)
        for i in array:
            copy_dag.remove_node(i)

        # Check graph connectivity
        num_components = nx.number_connected_components(copy_dag)

        # Calculate connected components (subgraphs) after removing the point
        subgraphs = list(nx.connected_components(copy_dag))

        # Store number of nodes in each subgraph
        result_list = []
        for subgraph in subgraphs:
            subgraph = sorted(subgraph)
            # if draw: print("Subgraph: ",subgraph)
            # --------------------------------------------------------------
            selected_elements = [circuit.queue[i - 1] for i in subgraph]
            # print(selected_elements)
            circuit_copy = copy.deepcopy(circuit)
            circuit_copy.queue = selected_elements
            # print(circuit_copy.draw())
            non_empty_qubits = del_empty_qubits(circuit_copy)
            non_empty_qubits.sort()
            # print(non_empty_qubits)
            # difference_list = [
            #     value - index for index, value in enumerate(non_empty_qubits)
            # ]
            result_list.append(len(non_empty_qubits))
            # --------------------------------------------------------------
            num_nodes.append(len(subgraph))

        # result_list = []
        # previous_sum = 0
        # for element in max_qubit:
        #    result_list.append(abs(element - previous_sum))
        #    previous_sum += element

        if max_qubits is None or has_number(result_list, max_qubits):
            if num_subcirucits is None or num_components == num_subcirucits:
                right_subgrafs.append(array)
                str_array = str(array)
                temp[str_array] = {}
                # num_nodes_in_subgraphs.append(num_nodes)
                comp_cost = (
                    2 ** len(array)
                    * num_components
                    * max(result_list)
                    * (max(num_nodes) - min(num_nodes))
                )
                computational_cost.append(comp_cost)
                temp[str_array]["computacional_cost"] = comp_cost
                temp[str_array]["num_components"] = num_components
                temp[str_array]["num_qubits"] = result_list
                temp[str_array]["num_gates"] = num_nodes
                if draw:
                    print("Removing gate ", array)
                if draw:
                    print("Number of connected components:", num_components)
                if draw:
                    print("Number of qubits per subcircuit: ", result_list)
                if draw:
                    print("Number of gates per subcircuit: ", num_nodes)
                if draw:
                    print("Computational Cost: ", comp_cost)
                if draw:
                    print_graph(copy_dag)

    # -------------------------------------------------------
    # DECIDE WHAT GATE TO CUT
    # -------------------------------------------------------
    print(temp)
    gate_cut = right_subgrafs[0]
    min_gate = min(computational_cost)
    if len(right_subgrafs) > 1:
        print("Options where to cut: ", right_subgrafs)
        # print("Num nodes each subgraph: ",num_nodes_in_subgraphs)

        # results = []
        # for option in num_nodes_in_subgraphs:
        # result = option[0]
        # Subtract each subsequent element from the result
        # for num in option[1:]:
        #    result -= num
        # results.append(max(option)-min(option))
        # print("Node difference between subgraphs: ",results)
        print("Costs: ", computational_cost)

        # min_gate = min(results)
        min_gate = min(computational_cost)
        index_of_smallest = computational_cost.index(min_gate)
        gate_cut = right_subgrafs[index_of_smallest]
        print("Choosing the smallest number")

    print("Gate where to cut for balanced subgraphs: ", gate_cut)
    print("Computational cost of cutting these gates: ", min_gate)
    return_list = [gate_cut, min_gate]
    return return_list


def find_edge_to_split_graph(graph):
    edges_cut = []
    blue_edges = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "blue"
    ]
    red_edges = [
        (edge[0], edge[1])
        for edge in graph.edges.data("color")
        if edge[2] == "red"
    ]
    for edge in blue_edges:
        # Remove an edge and check if it splits the graph
        # into disconnected components
        temp_graph = graph.copy()
        temp_graph.remove_edges_from(red_edges)
        temp_graph.remove_edge(*edge)
        if nx.number_connected_components(temp_graph) > 1:
            edges_cut.append(edge)
    return edges_cut


@task(returns=list)
def wire_selector(digraph, circuit, max_qubits, draw=False):
    # Finding an edge that splits the graph
    edge_to_remove = find_edge_to_split_graph(digraph)
    # print(edge_to_remove)
    print("\n")
    red_edges = [
        (edge[0], edge[1])
        for edge in digraph.edges.data("color")
        if edge[2] == "red"
    ]
    # num_nodes_in_subgraphs = []
    right_subgraph = []
    computational_cost = []
    for edge in edge_to_remove:
        num_nodes = []
        temp_graph = digraph.copy()
        temp_graph.remove_edges_from(red_edges)

        temp_graph.remove_edge(edge[0], edge[1])

        # Calculate connected components (subgraphs) after removing the point
        subgraphs = list(nx.connected_components(temp_graph))

        result_list = []
        for subgraph in subgraphs:
            subgraph = sorted(subgraph)
            # if draw: print("Subgraph: ",subgraph)
            # --------------------------------------------------------------
            selected_elements = [circuit.queue[i - 1] for i in subgraph]
            # print(selected_elements)
            circuit_copy = copy.deepcopy(circuit)
            circuit_copy.queue = selected_elements
            # print(circuit_copy.draw())
            non_empty_qubits = del_empty_qubits(circuit_copy)
            non_empty_qubits.sort()
            # print(non_empty_qubits)
            # difference_list = [
            #     value - index for index, value in enumerate(non_empty_qubits)
            # ]
            result_list.append(len(non_empty_qubits))
            # ---------------------------------------------------------------
            num_nodes.append(len(subgraph))

        # result_list = []
        # previous_sum = 0
        # for element in max_qubit:
        #    result_list.append(abs(element - previous_sum))
        #    previous_sum += element

        if max_qubits is None or has_number(result_list, max_qubits):
            # num_nodes_in_subgraphs.append(num_nodes)
            right_subgraph.append(edge)
            comp_cost = (
                10 * max(result_list) * (max(num_nodes) - min(num_nodes))
            )
            computational_cost.append(comp_cost)
            if draw:
                print(
                    f"Edge {edge} splits the graph into two \
                    disconnected components."
                )
            if draw:
                print("Number of qubits per subcircuit: ", result_list)
            if draw:
                print("Number of gates per subcircuit: ", num_nodes)
            if draw:
                print("Computational Cost: ", comp_cost)
            if draw:
                print_graph(temp_graph)

    wire_cut = right_subgraph
    min_wire = min(computational_cost)

    if len(right_subgraph) > 1:
        print("Options where to cut: ", right_subgraph)
        # print("Num nodes each subgraph: ",num_nodes_in_subgraphs)
        print("Costs: ", computational_cost)

        # min_wire = min(results)
        min_wire = min(computational_cost)
        index_of_smallest_wire = computational_cost.index(min_wire)
        wire_cut = [right_subgraph[index_of_smallest_wire]]
        wire_cost = computational_cost[index_of_smallest_wire]
        print("Choosing the smallest number")

    print("Wire where to cut for balanced subgraphs: ", wire_cut)
    print("Computational cost of cutting this wire: ", wire_cost)
    return_list = [wire_cut, min_wire]
    return return_list


def optimal_cut(
    circuit,
    max_qubits=None,
    num_subcircuits=None,
    max_cuts=None,
    gate_cut=True,
    wire_cut=True,
    draw=False,
    # max_subcircuit_size=None
):
    digraph2 = nx.Graph()
    dag2 = DAGgraph()

    build_dag(circuit, dag2)
    create_graph(dag2, digraph2)

    if gate_cut:
        return_list = gate_selector(
            digraph2, circuit, max_qubits, num_subcircuits, draw
        )
    if wire_cut:
        return_list2 = wire_selector(digraph2, circuit, max_qubits, draw)

    if gate_cut:
        return_list = compss_wait_on(return_list)
    if wire_cut:
        return_list2 = compss_wait_on(return_list2)

    if gate_cut:
        gate_cut2, min_gate2 = return_list
    if wire_cut:
        wire_cut2, min_wire2 = return_list2

    print("\n")
    if gate_cut and wire_cut:
        print("Minimum computational cost of gate cutting: ", min_gate2)
        print("Minimum computational cost of wire cutting: ", min_wire2)
        if min_gate2 > min_wire2:
            print("BEST OPTION WIRE: ", wire_cut2)
            return wire_cut2
        elif min_gate2 < min_wire2:
            print("BEST OPTION GATE: ", gate_cut2)
            return gate_cut2
        else:
            print("WIRE AND GATE SAME OPTION")
            return gate_cut2
    elif gate_cut:
        print("BEST GATE: ", gate_cut2)
        return gate_cut2
    elif wire_cut:
        print("BEST WIRE: ", wire_cut2)
        return wire_cut2
