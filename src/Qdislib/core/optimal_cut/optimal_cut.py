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
import copy
import qibo

'''from Qdislib.core.cutting_algorithms.wire_cutting import (
    wire_cutting,
    wire_cutting_QC,
)'''
from Qdislib.core.cutting_algorithms.gate_cutting import (
    gate_cutting,
    _gates_dict,
    _has_number,
    gate_cutting_QC,
)
from Qdislib.utils.graph import (
    _build_dag,
    _create_graph,
    _del_empty_qubits,
    print_graph,
    _DAGgraph,
)


def _double_gates(circuit, digraph, max_qubits, num_subcirucits, draw):
    """
    Analyze a circuit and its corresponding DAG to identify points
    where double gates can be applied. It calculates the computational cost and
    identifies suitable subcircuits based on the specified criteria.

    :param circuit: Circuit.
    :param digraph: Object DAG.
    :param max_qubits: int.
    :param num_subcircuits: int.
    :param draw: bool.
    :return: right_subgraphs, computational_cost

    """
    temp = {}
    right_subgrafs = []
    computational_cost = []
    _double_gates = _gates_dict(circuit)
    for _, array in _double_gates.items():
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
        if num_components > 1:
            for subgraph in subgraphs:
                subgraph = sorted(subgraph)
                if draw:
                    print("Subgraph: ", subgraph)
                # --------------------------------------------------------------
                selected_elements = [circuit.queue[i - 1] for i in subgraph]
                circuit_copy = copy.deepcopy(circuit)
                circuit_copy.queue = selected_elements
                non_empty_qubits = _del_empty_qubits(circuit_copy)
                non_empty_qubits.sort()
                result_list.append(len(non_empty_qubits))
                # --------------------------------------------------------------
                num_nodes.append(len(subgraph))
            if max_qubits is None or _has_number(result_list, max_qubits):
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
    return right_subgrafs, computational_cost


@task(returns=list)
def _gate_selector(
    digraph, circuit, max_qubits=None, num_subcirucits=None, draw=False
):
    """
    Select the appropriate gates to cut in order to split
    the circuit into balanced subcircuits based on the provided DAG and
    circuit information. It utilizes the `_double_gates` function to
    identify candidate points for gate cutting and then selects the
    optimal gate based on computational cost.

    :param digraph: Object DAG.
    :param circuit: Circuit.
    :param max_qubits: int.
    :param num_subcircuits: int.
    :param draw: bool.
    :return: gate_list
    """
    right_subgrafs, computational_cost = _double_gates(
        circuit, digraph, max_qubits, num_subcirucits, draw
    )
    gate_cut = right_subgrafs[0]
    min_gate = min(computational_cost)
    if len(right_subgrafs) > 1:
        if draw:
            print("Options where to cut: ", right_subgrafs)
        if draw:
            print("Costs: ", computational_cost)

        # min_gate = min(results)
        min_gate = min(computational_cost)
        index_of_smallest = computational_cost.index(min_gate)
        gate_cut = right_subgrafs[index_of_smallest]
        if draw:
            print("Choosing the smallest number")

    if draw:
        print("Gate where to cut for balanced subgraphs: ", gate_cut)
    if draw:
        print("Computational cost of cutting these gates: ", min_gate)
    gate_list = [gate_cut, min_gate]
    return gate_list


def _find_edge_to_split_graph(graph):
    """
    Identify edges in a directed acyclic graph (DAG) that can be
    cut to split the graph into multiple connected components. It specifically
    looks for edges marked as "blue" (interconected gates) and ensures that
    cutting them would result in more than one connected component.

    :param graph: Object DAG.
    :return: edges_cut.
    """
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
        temp_graph = graph.copy()
        temp_graph.remove_edges_from(red_edges)
        temp_graph.remove_edge(*edge)
        if nx.number_connected_components(temp_graph) > 1:
            edges_cut.append(edge)
    return edges_cut


@task(returns=list)
def _wire_selector(digraph, circuit, max_qubits, draw=False):
    """
    Select wires in a directed acyclic graph (DAG) to split the
    graph into multiple connected components. It identifies edges that, when
    removed, result in the graph being split into disconnected components, and
    selects the wire that minimizes computational cost.

    :param digraph: Object DAG.
    :param circuit: Circuit.
    :param max_qubits: int.
    :param draw: bool.
    :return: wire_list
    """
    edge_to_remove = _find_edge_to_split_graph(digraph)
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
            # --------------------------------------------------------------
            selected_elements = [circuit.queue[i - 1] for i in subgraph]
            circuit_copy = copy.deepcopy(circuit)
            circuit_copy.queue = selected_elements
            non_empty_qubits = _del_empty_qubits(circuit_copy)
            non_empty_qubits.sort()
            result_list.append(len(non_empty_qubits))
            # ---------------------------------------------------------------
            num_nodes.append(len(subgraph))
        if max_qubits is None or _has_number(result_list, max_qubits):
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
        if draw:
            print("Options where to cut: ", right_subgraph)
        # print("Num nodes each subgraph: ",num_nodes_in_subgraphs)
        if draw:
            print("Costs: ", computational_cost)

        # min_wire = min(results)
        min_wire = min(computational_cost)
        index_of_smallest_wire = computational_cost.index(min_wire)
        wire_cut = [right_subgraph[index_of_smallest_wire]]
        min_wire = computational_cost[index_of_smallest_wire]
        if draw:
            print("Choosing the smallest number")

    if draw:
        print("Wire where to cut for balanced subgraphs: ", wire_cut)
    if draw:
        print("Computational cost of cutting this wire: ", min_wire)
    wire_list = [wire_cut, min_wire]
    return wire_list


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
    """

    Description
    -----------
    Determine the optimal cut strategy for a given circuit, considering gate cutting and wire cutting options. It evaluates the computational cost of each cutting strategy and returns the one with the lowest cost.

    Parameters
    ----------
    circuit: Circuit.
        The circuit object to be analyzed.
    max_qubits: int, optional.
        The maximum number of qubits allowed in each subcircuit. Defaults to None.
    num_subcircuits: int, optional.
        The desired number of subcircuits after cutting. Defaults to None.
    max_cuts: int, optional.
        The maximum number of cuts allowed. Defaults to None.
    gate_cut: bool, optional.
        Whether to consider gate cutting as a strategy. Defaults to True.
    wire_cut: bool, optional
        Whether to consider wire cutting as a strategy. Defaults to True.
    draw: bool, optional.
        Whether to visualize the cutting process. Defaults to False.

    Returns
    -------
    best_cut: list
        The optimal cutting strategy, either gate cut or wire cut, with the lowest computational cost.

    Raises
    ------
    ValueError: Raised when there is an error in determining the best cutting strategy.

    Example
    -------
    >>> best_cut = optimal_cut(circuit, max_qubits=5, num_subcircuits=3,
    >>>                        gate_cut=True, wire_cut=True, draw=True)
    """
    digraph2 = nx.Graph()
    dag2 = _DAGgraph()

    _build_dag(circuit, dag2)
    _create_graph(dag2, digraph2)

    if gate_cut:
        return_list = _gate_selector(
            digraph2, circuit, max_qubits, num_subcircuits, draw
        )
    if wire_cut:
        return_list2 = _wire_selector(digraph2, circuit, max_qubits, draw)

    if gate_cut:
        return_list = compss_wait_on(return_list)
    if wire_cut:
        return_list2 = compss_wait_on(return_list2)

    if gate_cut:
        gate_cut2, min_gate2 = return_list
    if wire_cut:
        wire_cut2, min_wire2 = return_list2

    if gate_cut and wire_cut:
        if draw:
            print("Minimum computational cost of gate cutting: ", min_gate2)
        if draw:
            print("Minimum computational cost of wire cutting: ", min_wire2)
        if min_gate2 > min_wire2:
            if draw:
                print("BEST OPTION WIRE: ", wire_cut2)
            return wire_cut2
        elif min_gate2 < min_wire2:
            if draw:
                print("BEST OPTION GATE: ", gate_cut2)
            return gate_cut2
        else:
            if draw:
                print("WIRE AND GATE SAME OPTION")
            return gate_cut2
    elif gate_cut:
        if draw:
            print("BEST GATE: ", gate_cut2)
        return gate_cut2
    elif wire_cut:
        if draw:
            print("BEST WIRE: ", wire_cut2)
        return wire_cut2
    else:
        ValueError("Error on returning best option")


'''def execute_optimal_cut(
    observables, circuit, cut, shots=30000, chunk=1, verbose=False, sync=True, gpu=False, gpu_counter=0
):
    """

    Description
    -----------
    Execute the optimal cut strategy determined by the `optimal_cut` function. It applies either gate cutting or wire cutting to the circuit based on the specified cutting strategy.

    Parameters
    ----------
    circuit: Circuit.
        The circuit object to be modified.
    cut: list.
        The optimal cutting strategy determined by the `optimal_cut` function. It can be a list of gates or a tuple of gates for cutting the wire in the middle.
    shots: int.
        Number of shots that will be executed the circuit for the simulation.
    verbose: bool, optional.
        Whether to print additional information during execution. Defaults to False.

    Returns
    -------
    reconstruction: float
        The reconstructed circuit after applying the optimal cutting strategy.

    Example
    -------
    >>> reconstruction = execute_optimal_cut(circuit, cut=[2,14], verbose=True)
    """

    if isinstance(cut[0], tuple):
        reconstruction = wire_cutting(
            observables, circuit, cut, shots, verbose=verbose, sync=sync
        )
        if verbose:
            print("RECONSTRUCTION WIRE CUTTING: ", reconstruction)
    else:
        reconstruction = gate_cutting(
            observables,
            circuit,
            cut,
            shots,
            chunk=chunk,
            verbose=verbose,
            sync=sync,
            gpu=gpu,
            gpu_counter=gpu_counter
        )
        if verbose:
            print("RECONSTRUCTION GATE CUTTING: ", reconstruction)
    return reconstruction'''


'''def execute_optimal_cut_QC(
    connection, observables, circuit, cut, shots=30000, verbose=False
):
    """

    Description
    -----------
    Execute the optimal cut strategy determined by the `optimal_cut` function. It applies either gate cutting or wire cutting to the circuit based on the specified cutting strategy.

    Parameters
    ----------
    circuit: Circuit.
        The circuit object to be modified.
    cut: list.
        The optimal cutting strategy determined by the `optimal_cut` function. It can be a list of gates or a tuple of gates for cutting the wire in the middle.
    shots: int.
        Number of shots that will be executed the circuit for the simulation.
    verbose: bool, optional.
        Whether to print additional information during execution. Defaults to False.

    Returns
    -------
    reconstruction: float
        The reconstructed circuit after applying the optimal cutting strategy.

    Example
    -------
    >>> reconstruction = execute_optimal_cut(circuit, cut=[2,14], verbose=True)
    """

    if isinstance(cut[0], tuple):
        reconstruction = wire_cutting_QC(
            connection, observables, circuit, cut, shots, verbose=verbose
        )
        if verbose:
            print("RECONSTRUCTION WIRE CUTTING: ", reconstruction)
    else:
        reconstruction = gate_cutting_QC(
            connection, observables, circuit, cut, shots, verbose=verbose
        )
        if verbose:
            print("RECONSTRUCTION GATE CUTTING: ", reconstruction)
    return reconstruction'''
