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
from pycompss.api.parameter import COLLECTION_IN

import numpy as np
import qibo
from qibo import models, gates, hamiltonians, callbacks
from qibo.models import Circuit
from qibo.symbols import X, Y, Z, I

# for connecting with the quantum computer
from qiboconnection.connection import ConnectionConfiguration
from qiboconnection.api import API

from collections import Counter
from functools import reduce
from itertools import product

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
        print("Nodes: ", self.nodes)

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


def has_number_or_less(lst, number):
    for num in lst:
        if num <= number:
            return True
    return False


# ---------------------------------------------------
# Create and set DAG and nx.Graph()
# ---------------------------------------------------


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


def gates_dict(circuit):
    double_gates = {}
    for index, gate in enumerate(circuit.queue):
        if len(gate.qubits) > 1:
            start, end = gate.qubits
            new_tuple = [(i, i + 1) for i in range(start, end)]
            # print(new_tuple)
            for tuple in new_tuple:
                if tuple not in double_gates:
                    double_gates[tuple] = []
                    double_gates[tuple].append(index + 1)
                else:
                    double_gates[tuple].append(index + 1)

    # print(double_gates)
    return double_gates


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


def del_empty_qubits(circuit):
    empty_qubits = []
    for gate in circuit.queue:
        for i in gate.qubits:
            if i not in empty_qubits:
                empty_qubits.append(i)
        gate.qubits
    return empty_qubits


@task(returns=list)
def gen_graph_circuit(new_circuit):
    list_subcircuits = []
    # convert to DAG and DIGRPAPH
    digraph = nx.Graph()
    dag = DAGgraph()

    build_dag(new_circuit, dag)
    create_graph(dag, digraph)

    subgraphs = list(nx.connected_components(digraph))
    print(subgraphs)

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
        # print("Substracted list: ", subtracted_list)

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
    return list_subcircuits


def split_gates(gates_cut, circuit, draw=False):
    # ------------------------------------
    # SPLIT IN 4 SUBCIRCUITS
    # ------------------------------------
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    combinations_list = generate_combinations(len(gates_cut), type_gates)
    generated_circuits = []
    for index2, combination in enumerate(combinations_list):
        circuit1 = circuit.copy(True)
        target_gates = []
        for gate in gates_cut:
            target_gates.append(circuit1.queue[gate - 1])

        if any(
            type(element) != type(target_gates[0]) for element in target_gates
        ):
            print("All the gates to cut have to be the same type")
            return

        target_qubit = target_gates[0].control_qubits[0]

        for index, gate in enumerate(circuit1.queue):
            if gate in target_gates:
                if type(gate) == gates.CNOT:
                    # print("CNOT")
                    idx = target_gates.index(gate)
                    control_qubit = gate.control_qubits[0]
                    target_qubit = gate.target_qubits[0]
                    if len(combination[idx]) > 2:
                        circuit1.queue[index] = combination[idx][0](
                            control_qubit, np.pi / 2
                        )
                        circuit1.queue.insert(
                            index + 1, combination[idx][1](control_qubit)
                        )
                        circuit1.queue.insert(
                            index + 2,
                            combination[idx][2](target_qubit, np.pi / 2),
                        )
                        circuit1.queue.insert(
                            index + 3, combination[idx][3](target_qubit)
                        )
                    else:
                        circuit1.queue[index] = combination[idx][0](
                            control_qubit, np.pi / 2
                        )
                        circuit1.queue.insert(
                            index + 1,
                            combination[idx][1](target_qubit, np.pi / 2),
                        )
                elif type(gate) == gates.CZ:
                    # print("CZ")
                    idx = target_gates.index(gate)
                    control_qubit = gate.control_qubits[0]
                    target_qubit = gate.target_qubits[0]
                    circuit1.queue[index] = combination[idx][0](control_qubit)
                    circuit1.queue.insert(
                        index + 1, combination[idx][1](target_qubit)
                    )
                    # changed_gates.append(
                    #     (circuit1.queue[index],circuit1.queue[index+1])
                    # )

        if draw:
            print("\n Circuit " + str(index2 + 1))
            print(circuit1.draw())
        generated_circuits.append(circuit1)

    list_subcircuits = []
    for new_circuit in generated_circuits:
        new_list = gen_graph_circuit(new_circuit)
        list_subcircuits.append(new_list)

    list_subcircuits = compss_wait_on(list_subcircuits)
    list_subcircuits = concatenate_lists(list_subcircuits)

    if draw:
        for index, sub in enumerate(list_subcircuits):
            print("\n Subcircuit " + str(index + 1))
            print(sub.draw())
    return list_subcircuits


def concatenate_lists(lst):
    conc_list = []
    for el in lst:
        conc_list = conc_list + el
    return conc_list


@task(returns=qibo.states.CircuitResult)
def gate_simulation(i, shots=30000):
    # ------------------------------------
    # SIMULATION
    # ------------------------------------
    result = i(nshots=shots)
    return result


@task(returns=list)
def gate_frequencies(result):
    # ------------------------------------
    # FREQUENCIES
    # ------------------------------------
    freq = dict(result.frequencies(binary=True))
    return freq


@task(returns=float)
def gate_expectation_value(freq, shots=30000):
    # ------------------------------------
    # EXPECTATION VALUE
    # ------------------------------------
    expec = 0
    for key, value in freq.items():
        ones = key.count("1")
        if ones % 2 == 0:
            expec += float(value) / shots
        else:
            expec -= float(value) / shots
    return expec


def gate_reconstruction(type_gates, gates_cut, exp_values):
    # --------------------------------------
    # RECONSTRUCTION
    # --------------------------------------
    num_generated = int(len(exp_values) / 2 ** len(gates_cut))
    result = [
        eval("*".join(map(str, exp_values[i : i + num_generated])))
        for i in range(0, len(exp_values), num_generated)
    ]

    # result = [exp_values[i] * exp_values[i + 1]
    #               for i in range(0, len(exp_values), 2)]
    result1 = [x * 1j if i % 2 == 0 else x for i, x in enumerate(result)]
    result2 = [x * 1j if i % 2 != 0 else x for i, x in enumerate(result)]
    if type_gates == gates.CZ:
        reconstruction = (
            (1 / (1j + 1))
            * (sum(result1) + sum(result2))
            / 2 ** len(gates_cut)
        )
    elif type_gates == gates.CNOT:
        reconstruction = (
            1j
            * np.exp(-1j * np.pi / 4)
            / np.sqrt(2)
            * (sum(result1) + sum(result2))
            / 2 ** len(gates_cut)
        )
    print("\n")

    print("Reconstructed expected value: ", reconstruction)
    print("Absolute value of reconstruction ", np.absolute(reconstruction))
    print("Reconstruction value: ", reconstruction.real)
    # print("Anlytical expected value: ",exp_full_circuit)
    return reconstruction


@task(returns=dict, dicts=COLLECTION_IN)
def sum_dicts(dicts):
    summed_dict = reduce(lambda a, b: a + Counter(b), dicts, Counter())
    return dict(summed_dict)


def generate_combinations(n, gate_type):
    objects = []
    if gate_type == gates.CZ:
        objects = [(gates.S, gates.S), (gates.SDG, gates.SDG)]
    elif gate_type == gates.CNOT:
        objects = [
            (gates.RZ, gates.RX),
            (gates.RZ, gates.Z, gates.RX, gates.X),
        ]
    all_combinations = list(product(objects, repeat=n))
    return all_combinations


def gate_cutting(gates_cut, circuit, shots=30000, chunk=1, draw=False):
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    subcircuits = split_gates(gates_cut, circuit, draw)
    exp_value = []
    for i in subcircuits:
        list_freq = []
        i.add(gates.M(*range(i.nqubits)))
        for p in range(0, chunk):
            result = gate_simulation(i, int(shots / chunk))
            freq = gate_frequencies(result)
            # frq in a list COLLECTION
            list_freq.append(freq)
        # task per sumar dicts COLLECTIONS
        total_freq = sum_dicts(list_freq)
        exp_value.append(gate_expectation_value(total_freq, shots))

    exp_value = compss_wait_on(exp_value)
    # print("Expected values list: ",exp_value)
    result = gate_reconstruction(type_gates, gates_cut, exp_value)
    return result
