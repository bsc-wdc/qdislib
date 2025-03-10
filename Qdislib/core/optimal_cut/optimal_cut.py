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

"""Optimal cut algorithms."""

import itertools
import networkx
import typing

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *

from Qdislib.utils.graph_qibo import update_qubits
from Qdislib.utils.graph_qibo import update_qubits_serie
from Qdislib.utils.graph_qibo import remove_red_edges

import qibo
import qiskit

from Qdislib.utils.graph_qibo import circuit_to_dag
from Qdislib.utils.graph_qiskit import circuit_qiskit_to_dag, dag_to_circuit_qiskit

import qiskit.qasm2

from qiskit_addon_cutting.automated_cut_finding import (
    find_cuts,
    OptimizationParameters,
    DeviceConstraints,
    )


def find_nodes_with_qubit(
    graph: networkx.Graph,
    node: typing.Any,
    qubit: typing.Any,
    direction: str = "predecessor",
) -> typing.List[typing.Any]:
    """Find predecessor or successor nodes with specific qubit.

    :param graph: Graph.
    :param node: Node.
    :param qubit: Qbit.
    :param direction: Direction, defaults to "predecessor".
    :raises ValueError: Unexpected direction error. Must be "predecessor" or "successor".
    :return: List of nodes with qubit.
    """
    if direction == "predecessor":
        neighbors = graph.predecessors(node)
    elif direction == "successor":
        neighbors = graph.successors(node)
    else:
        raise ValueError("Direction must be either 'predecessor' or 'successor'")

    # Filter neighbors based on the qubit data
    nodes_with_qubit = [n for n in neighbors if qubit in graph.nodes[n]["qubits"]]
    return nodes_with_qubit

@task(returns=2)
def evaluate_cut(
    graph: networkx.Graph,
    cut_edges: typing.List[typing.Any],
    cut_nodes: typing.List[typing.Any],
    threshold: int,
) -> typing.Tuple[bool, float]:
    """Evaluate if the cut is feasible and compute the score.

    :param graph: NetworkX graph object.
    :param cut_edges: List of edges to remove for this cut.
    :param cut_nodes: List of nodes to remove for this cut (treated as edges).
    :param threshold: The max allowable difference between min and max qubit in each subgraph.
    :return: (valid_cut: bool, score: float)
    """
    # Create a copy of the graph and remove the edges
    graph_copy = graph.copy()

    # Process cut_nodes for multi-qubit gates only and convert them to edge cuts
    for elem in cut_nodes:
        target_qubits = graph_copy.nodes[elem]["qubits"]

        # Only consider multi-qubit gates (nodes with more than 1 qubit)
        if len(target_qubits) > 1:
            pred_0 = find_nodes_with_qubit(
                graph_copy, elem, qubit=target_qubits[0], direction="predecessor"
            )
            pred_1 = find_nodes_with_qubit(
                graph_copy, elem, qubit=target_qubits[1], direction="predecessor"
            )
            succ_0 = find_nodes_with_qubit(
                graph_copy, elem, qubit=target_qubits[0], direction="successor"
            )
            succ_1 = find_nodes_with_qubit(
                graph_copy, elem, qubit=target_qubits[1], direction="successor"
            )

            # Collect edges connected to the multi-qubit nodes and treat them as wire cuts
            if pred_0:
                cut_edges.append((pred_0[0], elem))
            if pred_1:
                cut_edges.append((pred_1[0], elem))
            if succ_0:
                cut_edges.append((elem, succ_0[0]))
            if succ_1:
                cut_edges.append((elem, succ_1[0]))

    # Remove edges from the graph copy
    graph_copy.remove_edges_from(cut_edges)

    # Find all connected components after the cut
    components = [
        graph_copy.subgraph(c).copy()
        for c in networkx.connected_components(graph_copy.to_undirected())
    ]
    if len(components) < 2:
        return False, float("inf")

    num_nodes = []
    for component in components:
        component, _, _ = update_qubits_serie(component)
        highest_qubit = float("-inf")
        smallest_qubit = float("inf")
        for node in component:
            qubits = component.nodes[node]["qubits"]
            if max(qubits) > highest_qubit:
                highest_qubit = max(qubits)
            if min(qubits) < smallest_qubit:
                smallest_qubit = min(qubits)
        if highest_qubit - smallest_qubit > (threshold - 1):
            return False, float(
                "inf"
            )  # Invalid cut if any component exceeds the threshold
        num_nodes.append(len(component))

    cut_size = len(cut_edges)
    num_components = len(components)
    diff_num_nodes = abs(max(num_nodes) - min(num_nodes))

    # Weights for the objective function
    w1 = 2  # Weight for cut size (minimize this)
    w2 = -1  # Weight for number of components (maximize this)
    w3 = 1  # Weight for balanced graphs (minimize this)

    # Calculate the score
    score = w1 * cut_size + w2 * num_components + w3 * diff_num_nodes
    return True, score

#@task(returns=3)
def optimal_cut_wire(
    graph: networkx.Graph, threshold=None, verbose=False
) -> typing.Tuple[typing.List[typing.Any], typing.List[typing.Any], int]:
    """Find the best cut based on the given constraints.

    :param graph: NetworkX graph object.
    :param threshold: The max allowable difference between min and max qubit in each subgraph.
    :return: Best cut as a list of edges and its score.
    """
    best_cut_edges = None
    best_score_components = []
    best_cut_components = []

    graph = remove_red_edges(graph)
    components = [
        graph.subgraph(c).copy()
        for c in networkx.connected_components(graph.to_undirected())
    ]
    
    for idx, component in enumerate(components):
        best_score = []
        best_cut_edges = []
        component, highest_qubit, smallest_qubit = update_qubits_serie(component)
        highest_qubit = highest_qubit - 1
        if verbose:
            print(f"highest_qubit: {highest_qubit}")
            print(f"smallest_qubit: {smallest_qubit}")

        if threshold is None:
            threshold = 1000

        if highest_qubit - smallest_qubit > (threshold - 1):
            if verbose:
                print(f"Component {idx} out of {len(components)}")
                print("Cut required due to threshold violation")

            centrality = networkx.edge_betweenness_centrality(component)
            sorted_edges = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
            edges = [edge for edge, _ in sorted_edges]

            if len(edges) > 25:
                edges = edges[:25]

            # Filter articulation points that are multi-qubit gates only
            articulation_points = [
                node
                for node in networkx.articulation_points(component.to_undirected())
                if len(component.nodes[node]["qubits"]) > 1  # Only multi-qubit gates
            ]

            # TODO: This variable is not used. Remove?
            #flag_best_score = False

            for r in range(1, 8 + 1):  # Number of edges to remove
                '''if flag_best_score:
                    break'''

                if r <= 3:  # Only edges
                    for cut_edges in itertools.combinations(edges, r):
                        valid, score = evaluate_cut(
                            component, list(cut_edges), [], threshold
                        )
                        #if valid:  # and abs(score) < best_score:
                        best_cut_edges.append(cut_edges)
                        best_score.append(score)
                        # Loop break condition
                        '''if score < 2:
                            flag_best_score = True
                            break'''

                elif (
                    r >= 4
                ):  # Combination of node and edge removal (nodes converted to edges)
                    num_nodes_to_remove = r // 4
                    num_edges_to_remove = r % 4
                    if len(articulation_points) >= num_nodes_to_remove:
                        for cut_nodes in itertools.combinations(
                            articulation_points, num_nodes_to_remove
                        ):
                            for cut_edges in itertools.combinations(
                                edges, num_edges_to_remove
                            ):
                                all_cut_edges = list(cut_edges)  # Start with the edges
                                valid, score = evaluate_cut(
                                    component, all_cut_edges, cut_nodes, threshold
                                )
                                #if valid:  # and abs(score) < best_score:
                                best_cut_edges.append(cut_edges)
                                best_score.append(score)
                                # Loop break condition
                                '''if score < 2:
                                    flag_best_score = True
                                    break'''

            best_score_components.append(best_score)
            best_cut_components.append(best_cut_edges)

        else:
            if verbose:
                print(f"Component {idx} out of {len(components)}")
                print("No cut required")

    return best_score_components, best_cut_components
    

#@task(returns=2)
def optimal_cut_gate(dag, max_qubits=None, max_components=None, max_cuts=None, verbose=False):
        double_gates = []
        for node, data in dag.nodes(data=True):
            if len(data["qubits"]) > 1:
                double_gates.append(node)
        
        if verbose:
            print(double_gates)

        print(double_gates)

        centrality = networkx.betweenness_centrality(dag)
        sorted_nodes = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
        double_gates = [nodes for nodes, _ in sorted_nodes if nodes in double_gates]

        print(sorted_nodes)

        print(double_gates)

        if len(double_gates) > 10:
            #(len(double_gates)//2)
            double_gates = double_gates[:(len(double_gates)//4)]


        if max_cuts is None:
            max_cuts = 8
        
        '''if max_components is None:
            max_components = float('inf')'''


        results = []
        final_cut = []


        for r in range(1,max_cuts):
            for cut in itertools.combinations(double_gates, r):
                score, cut = evaluate_cut_gate(dag,cut,max_components, max_qubits, verbose)
                results.append(score)
                final_cut.append(cut)

        return results, final_cut

        

@task(returns=2)
def evaluate_cut_gate(dag,cut,max_components, max_qubits, verbose=False):
    flag = True
    if verbose:
        print("TRYING CUT ", cut)
    copy_dag = dag.copy()
    copy_dag.remove_nodes_from(cut)

    S = [copy_dag.subgraph(c).copy() for c in networkx.connected_components(copy_dag.to_undirected())]
    num_components = len(S)

    if num_components <= 1:
        return float('inf'), float('inf')

    max_num_qubits = []
    for s in S:
        s_new, highest_qubit, smallest_qubit = update_qubits_serie(s)
        max_num_qubits.append(highest_qubit - smallest_qubit)


    # Weights for the objective function
    w1 = 3  # Weight for cut size (minimize this)
    w2 = -1  # Weight for number of components (maximize this)
    w3 = 1  # Weight for balanced graphs (minimize this)

    # Calculate the score
    score = w1 * len(cut) + w2 * num_components + w3 * max(max_num_qubits)
    if verbose:
        print("NUM COMPONENTS ", num_components)
        print("NUM MAX QUBITS ", max(max_num_qubits))
        print("SCORE ", score)



    if max_components is not None and max_components < num_components:
        flag = False


    if max_qubits is not None and max(max_num_qubits) > max_qubits:
        flag = False

    #print(max(max_num_qubits) > max_qubits)
    if flag:
        #results[counter] = score
        #final_cut[counter] = cut
        return score, cut
    else:
        return float('inf'), float('inf')

def optimal_cut(circuit, max_qubits=None, max_components=None, max_cuts=None, wire_cut=True, gate_cut=True, implementation='qdislib', verbose=False):
    if type(circuit) == qiskit.circuit.quantumcircuit.QuantumCircuit:
        dag = circuit_qiskit_to_dag(circuit)
        num_qubits = circuit.num_qubits

    elif type(circuit) == qibo.models.Circuit:
        dag = circuit_to_dag(circuit)
        num_qubits = circuit.nqubits

    elif type(circuit) == networkx.DiGraph:
        dag = circuit
    
    else:
        Exception("Type circuit not suported")

    if implementation == 'qdislib':
        if wire_cut:
            if max_qubits is None:
                max_qubits_wire_cut = num_qubits //2 +1
            else:
                max_qubits_wire_cut = max_qubits
            best_score_components, best_cut_components = optimal_cut_wire(dag, max_qubits_wire_cut, verbose)

            
        if gate_cut:
            '''if max_qubits is None:
                max_qubits = num_qubits'''
            results, final_cut = optimal_cut_gate(dag, max_qubits, max_components, max_cuts, verbose)


        #if gate_cut and wire_cut:
        best_score_components = compss_wait_on(best_score_components)
        best_cut_components = compss_wait_on(best_cut_components)
        results = compss_wait_on(results)
        final_cut = compss_wait_on(final_cut)

        

        if not wire_cut:
            cuts, scores, max_len_cut = None, [float('inf')], float('inf')
        else:
            if verbose:
                print(best_score_components)
                print(best_cut_components)

            cuts = []
            scores = []
            max_len_cut = float("-inf")

            if best_score_components != [[]]:
                for idx, best_score_comp in enumerate(best_score_components):
                    best_score = [abs(ele) for ele in best_score_comp]
                    index_min = best_score.index(min(best_score))
                    best_cut_edges = best_cut_components[idx][index_min]
                    max_len_cut = len(best_cut_components[idx])
                    best_score_min = min(best_score)
                    cuts = cuts + [*best_cut_edges]
                    scores.append(best_score_min)

            if verbose:
                print(scores)
                print(cuts)
        
        if not gate_cut:
            best_gate_score, cut_gate = float('inf'), None
        else:
            if verbose:
                print(results)

            if results == []:
                best_gate_score, cut_gate = float('inf'), []
            else:
                best_score_gate = min(results)
                min_index = results.index(best_score_gate)
                min_cut = final_cut[min_index]

                best_gate_score, cut_gate = best_score_gate, [min_cut]


        #cuts = compss_wait_on(cuts)
        #scores = compss_wait_on(scores)
        #max_len_cut = compss_wait_on(max_len_cut)
        #cut_gate = compss_wait_on(cut_gate)  
        #best_gate_score = compss_wait_on(best_gate_score)  

        scores = sum(scores)

        if verbose:
            print(cuts, scores, max_len_cut)
            print(best_gate_score, cut_gate)

        if not cuts:
            scores = float('inf')

        if scores < best_gate_score:
            print("Wire Cutting best cut: ", cuts)
            return cuts
        if best_gate_score < scores:
            print("Gate Cutting best cut: ", cut_gate)
            return cut_gate
        else:
            print("Same score both cuts, choosing gate cutting:  ", cut_gate)
            return cut_gate
    
    elif implementation == 'ibm-ckt':
        if type(circuit) == qibo.models.Circuit:
            qasm_str = circuit.to_qasm()
            circuit = qiskit.qasm2.loads(qasm_str)
        elif type(circuit) == networkx.DiGraph:
            circuit = dag_to_circuit_qiskit(circuit)

        # Specify settings for the cut-finding optimizer
        optimization_settings = OptimizationParameters(seed=50)

        # Specify the size of the QPUs available
        if max_qubits is None:
            max_qubits_ckt = num_qubits //2 +1
        else:
            max_qubits_ckt = max_qubits
        device_constraints = DeviceConstraints(qubits_per_subcircuit=max_qubits_ckt)

        cut_circuit, metadata = find_cuts(circuit, optimization_settings, device_constraints)
        if verbose:
            print(
                f'Found solution using {len(metadata["cuts"])} cuts with a sampling '
                f'overhead of {metadata["sampling_overhead"]}.\n'
                f'Lowest cost solution found: {metadata["minimum_reached"]}.'
            )
            for cut in metadata["cuts"]:
                print(f"{cut[0]} at circuit instruction index {cut[1]}")
            cut_circuit.draw("mpl", scale=0.8, fold=-1)
        
        list_gates_cut=[]
        for cut in metadata["cuts"]:
            if cut[0] == 'Gate Cut':
                gate_name = circuit.data[cut[1]].operation.name.upper()
                list_gates_cut.append(f"{gate_name}_{cut[1]+1}")
        print(list_gates_cut)
        return list_gates_cut



