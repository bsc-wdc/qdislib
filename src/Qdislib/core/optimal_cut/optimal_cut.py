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
import itertools as it

from Qdislib.utils.graph import update_qubits, remove_red_edges


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

def evaluate_cut(graph, cut_edges,cut_nodes, threshold):
    """
    Evaluate if the cut is feasible and compute the score.
    :param graph: NetworkX graph object.
    :param cut_edges: List of edges to remove for this cut.
    :param threshold: The max allowable difference between min and max qubit in each subgraph.
    :return: (valid_cut: bool, score: float, num_components: int)
    """
    # Create a copy of the graph and remove the edges
    graph_copy = graph.copy()

    cut_nodes_edges = []

    for elem in cut_nodes:
        target_qubits = graph_copy.nodes[elem]['qubits']

        # Find predecessor node with qubit 1
        pred_0 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[0], direction='predecessor')

        print(elem ,target_qubits )
        # Find predecessor node with qubit 2
        pred_1 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[1], direction='predecessor')

        # Find successor node with qubit 1
        succ_0 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[0], direction='successor')

        # Find successor node with qubit 2
        succ_1 = find_nodes_with_qubit(graph_copy, elem, qubit=target_qubits[1], direction='successor')

        # Output the results
        #print(f"Predecessor nodes with qubit {target_qubits[0]}: {pred_0}")
        #print(f"Predecessor nodes with qubit {target_qubits[1]}: {pred_1}")
        #print(f"Successor nodes with qubit {target_qubits[0]}: {succ_0}")
        #print(f"Successor nodes with qubit {target_qubits[1]}: {succ_1}")
        if pred_0:
            cut_nodes_edges.append((pred_0[0], elem))
        if pred_1:
            cut_nodes_edges.append((pred_1[0], elem))
        if succ_0:
            cut_nodes_edges.append((elem, succ_0[0]))
        if succ_1:
            cut_nodes_edges.append((elem, succ_1[0]))


    #plot_dag(graph_copy)
    graph_copy.remove_edges_from(cut_edges)
    graph_copy.remove_nodes_from(cut_nodes)
    #plot_dag(graph_copy)
    # Find all connected components after the cut
    components = [graph_copy.subgraph(c).copy() for c in nx.connected_components(graph_copy.to_undirected())]
    if len(components) < 2:
        return False, float('inf'), 0
    #print(components)
    # Check the qubit difference constraint in each component



    num_nodes = []
    for component in components:
        #plot_dag(component)
        component, highest_qubit = update_qubits(component)
        #print("HIGHEST: ", highest_qubit)
        #print(dag_to_circuit(component,5)[0].draw())
        highest_qubit = float('-inf')
        smallest_qubit = float('inf')
        for node in component:
            qubits = component.nodes[node]['qubits']
            #print(qubits)
            if max(qubits) > highest_qubit:
                highest_qubit = max(qubits)
            if min(qubits) < smallest_qubit:
                smallest_qubit = min(qubits)
        #print(highest_qubit)
        #print(smallest_qubit)
        if highest_qubit - smallest_qubit > (threshold-1):
            #print(highest_qubit)
            #print(smallest_qubit)
            #print(threshold)
            return False, float('inf'), 0  # Invalid cut if any component exceeds the threshold
        num_nodes.append(len(component))
    # Score calculation (minimize edges removed, maximize components)

    print("LEN CUT NODES: ", len(cut_nodes_edges))
    print("LEN CUT NODES: ", cut_nodes_edges)
    cut_size = len(cut_edges) + (len(cut_nodes_edges))
    num_components = len(components)
    diff_num_nodes = abs(max(num_nodes) - min(num_nodes))
    
    # Weights for the objective function
    w1 = 1  # Weight for cut size (we want to minimize this)
    w2 = -1  # Weight for number of components (we want to maximize this)
    w3 = 1 #Weight for balanced graphs with num of nodes (we want to minimize this)
    
    # Calculate the score
    score = w1 * cut_size + w2 * num_components + w3 * diff_num_nodes
    print("SCORE ", score)
    return True, score, cut_nodes_edges

def find_best_cut(graph, threshold):
    """
    Find the best cut based on the given constraints.
    :param graph: NetworkX graph object.
    :param threshold: The max allowable difference between min and max qubit in each subgraph.
    :return: Best cut as a list of edges and its score.
    """
    best_cut_edges = None
    best_cut_nodes = None
    best_score = float('inf')


    graph = remove_red_edges(graph)
    # Find all connected components after the cut
    components = [graph.subgraph(c).copy() for c in nx.connected_components(graph.to_undirected())]
    cuts=[]
    scores=[]
    for component in components:
        highest_qubit = float('-inf')
        smallest_qubit = float('inf')
        for node in component:
            qubits = component.nodes[node]['qubits']
            if max(qubits) > highest_qubit:
                highest_qubit = max(qubits)
            if min(qubits) < smallest_qubit:
                smallest_qubit = min(qubits)    
        if highest_qubit - smallest_qubit > (threshold-1):
            print("Yes we need cut")
            # Try all combinations of edge removals

            #EDGES BETWENESS CENTRALITY
            # Calculate edge betweenness centrality for this component
            centrality = nx.edge_betweenness_centrality(component)
            
            # Sort edges by betweenness centrality in descending order (most central edges first)
            sorted_edges = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
            
            # Try combinations of edge removals, prioritizing high centrality edges
            edges = [edge for edge, _ in sorted_edges]

            initial_articulation_points = list(nx.articulation_points(component.to_undirected()))

            articulation_points = []
            for elem in initial_articulation_points:
                print(component.nodes[elem]['qubits'])
                if len(component.nodes[elem]['qubits'])> 1:
                    articulation_points.append(elem)

            print("ARTICULATION POINTS: ", articulation_points)
            #edges = list(component.edges())
            flag_best_score = False
            counter = 0
            for r in range(1, 8+1):  # Number of edges to remove
                '''for cut_edges in it.combinations(edges, r):
                    print(cut_edges)
                    #plot_dag(component)
                    valid, score, _ = evaluate_cut(component, cut_edges, threshold)
                    print(valid, score)
                    if valid and abs(score) < best_score:
                        best_cut = cut_edges
                        print(best_cut)
                        best_score = abs(score) 
                        print("HEY ", best_score)  '''
                if flag_best_score:
                    break

                if r <= 3:  # Only edges
                    for cut_edges in it.combinations(edges, r):
                        valid, score, _ = evaluate_cut(component, cut_edges, [], threshold)
                        
                        print("HERE")
                        if valid and abs(score) < best_score:
                            best_cut_edges = cut_edges
                            best_cut_nodes = []
                            best_score = abs(score)
                            counter += 1
                            print("BEST SCORE ", best_score)
                            if best_score < 2:
                                print("COUNTER: ",counter +1)
                                flag_best_score = True
                                break
                    if flag_best_score:
                        break

                elif r >= 4:  # Combination of node and edge removal
                    # Remove nodes for r >= 4
                    num_nodes_to_remove = r // 4  # Each 4 cuts = 1 node removed
                    num_edges_to_remove = r % 4  # Remaining cuts are edges
                    print("HERE 2")
                    if len(articulation_points) >= num_nodes_to_remove:  # Check if enough nodes to remove
                        for cut_nodes in it.combinations(articulation_points, num_nodes_to_remove):
                            for cut_edges in it.combinations(edges, num_edges_to_remove):
                                valid, score, cut_nodes_edges = evaluate_cut(component, cut_edges, cut_nodes, threshold)
                                
                                if valid and abs(score) < best_score:
                                    best_cut_edges = cut_edges
                                    #best_cut_nodes = cut_nodes
                                    
                                    best_cut_nodes = cut_nodes_edges
                                    print(best_cut_nodes)
                                    best_score = abs(score)
                                    counter += 1
                                    if best_score < 2:
                                        print("COUNTER: ",counter +1)
                                        flag_best_score = True
                                        break
                            if flag_best_score:
                                break
                        if flag_best_score:
                            break

                else:
                    print("Something went wrong")
                    raise ValueError


            cuts = cuts + [*best_cut_edges] + best_cut_nodes
            scores.append(best_score)
        else:
            print("NO we dont need cut")

    return cuts, scores