#!/usr/bin/env python3
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
import itertools as it

from Qdislib.utils.graph import update_qubits, update_qubits_serie, remove_red_edges


# Function to find predecessor or successor nodes with specific qubit
def find_nodes_with_qubit(G, node, qubit, direction='predecessor'):
    if direction == 'predecessor':
        neighbors = G.predecessors(node)
    elif direction == 'successor':
        neighbors = G.successors(node)
    else:
        raise ValueError("Direction must be either 'predecessor' or 'successor'")
    
    # Filter neighbors based on the qubit data
    nodes_with_qubit = [n for n in neighbors if qubit in G.nodes[n]['qubits']]
    return nodes_with_qubit

@task(returns=2)
def evaluate_cut(graph, cut_edges, cut_nodes, threshold):
    """
    Evaluate if the cut is feasible and compute the score.
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
        target_qubits = graph_copy.nodes[elem]['qubits']

        # Only consider multi-qubit gates (nodes with more than 1 qubit)
        if len(target_qubits) > 1:
            pred_0 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[0], direction='predecessor')
            pred_1 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[1], direction='predecessor')
            succ_0 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[0], direction='successor')
            succ_1 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[1], direction='successor')

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
    components = [graph_copy.subgraph(c).copy() for c in nx.connected_components(graph_copy.to_undirected())]
    if len(components) < 2:
        return False, float('inf')

    num_nodes = []
    for component in components:
        component, _ , _= update_qubits_serie(component)
        highest_qubit = float('-inf')
        smallest_qubit = float('inf')
        for node in component:
            qubits = component.nodes[node]['qubits']
            if max(qubits) > highest_qubit:
                highest_qubit = max(qubits)
            if min(qubits) < smallest_qubit:
                smallest_qubit = min(qubits)
        if highest_qubit - smallest_qubit > (threshold - 1):
            return False, float('inf')  # Invalid cut if any component exceeds the threshold
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


def optimal_cut(graph, threshold):
    """
    Find the best cut based on the given constraints.
    :param graph: NetworkX graph object.
    :param threshold: The max allowable difference between min and max qubit in each subgraph.
    :return: Best cut as a list of edges and its score.
    """
    best_cut_edges = None
    best_score_components = []
    best_cut_components = []

    graph = remove_red_edges(graph)
    components = [graph.subgraph(c).copy() for c in nx.connected_components(graph.to_undirected())]
    cuts = []
    scores = []
    max_len_cut = float("-inf")
    for idx, component in enumerate(components):
        best_score = []
        best_cut_edges = []
        component, highest_qubit, smallest_qubit = update_qubits_serie(component)
        highest_qubit = highest_qubit - 1
        print(highest_qubit)
        print(smallest_qubit)
        '''highest_qubit = float('-inf')
        smallest_qubit = float('inf')
        for node in component:
            qubits = component.nodes[node]['qubits']
            if max(qubits) > highest_qubit:
                highest_qubit = max(qubits)
            if min(qubits) < smallest_qubit:
                smallest_qubit = min(qubits)'''

        if highest_qubit - smallest_qubit > (threshold - 1):
            print(f"Component {idx} out of {len(components)}")
            print("Cut required due to threshold violation")

            centrality = nx.edge_betweenness_centrality(component)
            sorted_edges = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
            edges = [edge for edge, _ in sorted_edges]

            # Filter articulation points that are multi-qubit gates only
            articulation_points = [
                node for node in nx.articulation_points(component.to_undirected())
                if len(component.nodes[node]['qubits']) > 1  # Only multi-qubit gates
            ]

            flag_best_score = False

            for r in range(1, 8 + 1):  # Number of edges to remove
                if flag_best_score:
                    break

                if r <= 3:  # Only edges
                    for cut_edges in it.combinations(edges, r):
                        valid, score = evaluate_cut(component, list(cut_edges), [], threshold)
                        if valid: #and abs(score) < best_score:
                            best_cut_edges.append(cut_edges)
                            best_score.append(score)
                            '''if best_score < 2:
                            flag_best_score = True
                            break'''

                elif r >= 4:  # Combination of node and edge removal (nodes converted to edges)
                    num_nodes_to_remove = r // 4
                    num_edges_to_remove = r % 4
                    if len(articulation_points) >= num_nodes_to_remove:
                        for cut_nodes in it.combinations(articulation_points, num_nodes_to_remove):
                            for cut_edges in it.combinations(edges, num_edges_to_remove):
                                all_cut_edges = list(cut_edges)  # Start with the edges
                                valid, score = evaluate_cut(component, all_cut_edges, cut_nodes, threshold)
                                if valid: #and abs(score) < best_score:
                                    best_cut_edges.append(cut_edges)
                                    best_score.append(score)
                                '''if best_score < 2:
                                    flag_best_score = True
                                    break'''

            best_score_components.append(best_score)
            best_cut_components.append(best_cut_edges)

            
        else:
            print(f"Component {idx} out of {len(components)}")
            print("No cut required")

    best_score_components = compss_wait_on(best_score_components)
    best_cut_components = compss_wait_on(best_cut_components)
    print(best_score_components)
    print(best_cut_components)
    
    for idx,best_score in enumerate(best_score_components):
        best_score = [abs(ele) for ele in best_score]
        index_min = best_score.index(min(best_score))
        best_cut_edges = best_cut_components[idx][index_min]
        max_len_cut = len(best_cut_components[idx])
        best_score = min(best_score)
        cuts = cuts + [*best_cut_edges]
        scores.append(best_score)

    print(scores)
    print(cuts)
    return cuts, scores, max_len_cut
