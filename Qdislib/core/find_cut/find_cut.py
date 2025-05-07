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

"""Find cut algorithms."""

import networkx
import typing

try:
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
    pycompss_available = True
except ImportError:
    print("NO PYCOMPSS AVAILABLE")
    # Define dummy decorators and functions to avoid breaking the code
    def task(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

    def compss_wait_on(obj):
        return obj

    pycompss_available = False

from Qdislib.utils.graph_qibo import _update_qubits
from Qdislib.utils.graph_qibo import _update_qubits_serie
from Qdislib.utils.graph_qibo import _remove_red_edges

import qibo
import qiskit
import pymetis

from Qdislib.utils.graph_qibo import circuit_qibo_to_dag
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

def circuit_to_qubit_graph(circuit):
    dag = networkx.MultiGraph()
    for i in range(0,circuit.nqubits):
        dag.add_node(i,gate=i)
    
    for idx, gate in enumerate(circuit.queue,start=1):
        if len(gate.qubits) > 1:
            dag.add_edge(*gate.qubits, gate=f'{gate.__class__.__name__}_{idx}')
    
    return dag

def circuit_to_gate_graph(circuit):
    """Convert a Qibo circuit to a DAG where each node stores gate information.

    :param circuit: The Qibo circuit to transform.
    :return: A directed acyclic graph (DAG) with nodes containing gate information.
    """
    # Create a directed graph
    dag = networkx.MultiGraph()
    
    counter = 0
    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.queue, start=1):

        if len(gate.qubits) < 2:
                continue

        # Unique identifier for each gate instance
        gate_name = f"{gate.__class__.__name__}_{gate_idx}".upper()

        # Add the gate to the DAG, including the gate type, qubits, and parameters
        dag.add_node(
            counter,gate=gate_name,qubits=gate.qubits
        )

        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
            for pred_gate in reversed(list(dag.nodes)):
                # Check if the qubit is in the node's qubits
                pred_name = dag.nodes[pred_gate]["gate"]
                if (
                    dag.nodes[pred_gate].get("qubits")
                    and qubit in dag.nodes[pred_gate]["qubits"]
                    and counter != pred_gate
                ):
                    dag.add_edge(pred_gate, counter, gate=(f'{pred_name}',f'{gate_name}'))
                    break
        counter += 1

    return dag


def find_best_kernighan_lin_partition(graph, max_qubits, max_cuts=8, num_runs=10):
    best_cutset_size = float('inf')
    best_partition = []
    best_cutset = []
    
    for _ in range(num_runs):
        #print(len(graph.nodes))
        # Step 2: Perform balanced partition using Kernighan-Lin
        partition = networkx.algorithms.community.kernighan_lin_bisection(graph.to_undirected())
        
        # Step 3: Extract the edges in the cut-set
        cutset = list(networkx.edge_boundary(graph, partition[0], partition[1], data=True))
        cutset_size = len(cutset)
        #cutset = [value[2]['gate'] for value in cutset]
        flag = False
        for ele in sorted(list(c) for c in partition):
            if max_qubits and len(ele) > max_qubits:
                flag=True
        if max_cuts and cutset_size > max_cuts:
            flag=True
        if flag:
            break

        # If this cut-set is smaller, update the best partition and cut-set
        if cutset_size < best_cutset_size:
            best_cutset_size = cutset_size
            best_partition = partition
            best_cutset = cutset
        best_cutset

    if best_cutset is not []:
        best_cutset = [el[2]['gate'] for el in best_cutset]
        #best_cutset = [(u,v) for u,v,c in best_cutset]
    else:
        return [],[]
    return best_partition, best_cutset

def find_best_girvan_newman_cutset(graph, max_qubits, max_cuts):
    from networkx.algorithms.community import girvan_newman
    # Perform Girvan-Newman community detection
    communities = girvan_newman(graph)

    # Get multiple partitions at different levels of edge removal
    for i, partition in enumerate(communities):
        partition_list = list(partition)  # Convert partition to list of sets
        #print(f"Partition {i+1}: {tuple(sorted(list(c) for c in partition_list))}")

        # Find edges that are between different partitions
        cut_edges = []
        for u, v in graph.edges():
            # Check if nodes u and v belong to different components
            component_u = [c for c in partition_list if u in c][0]
            component_v = [c for c in partition_list if v in c][0]
            if component_u != component_v:
                data = graph.get_edge_data(u,v)
                gates_list = [value['gate'] for value in data.values()]
                
                
                #gates_list = [(data[0],data[1])]

                cut_edges = cut_edges + gates_list
                #cut_edges = cut_edges + [(u,v)]

        #print(f"Edges that need to be cut for Partition {i+1}: {cut_edges}")
        flag = True
        for ele in sorted(list(c) for c in partition_list):
            if max_qubits and len(ele) > max_qubits:
                flag=False
        if flag:
            break

    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []

    return partition_list, cut_edges

def find_best_spectral_clustering(graph,clusters,max_qubits,max_cuts):
    from sklearn.cluster import SpectralClustering
    # Convert the graph to adjacency matrix form
    adjacency_matrix = networkx.to_numpy_array(graph)
    # Perform spectral clustering
    num_clusters = clusters  # Example: Find 3 components
    sc = SpectralClustering(n_clusters=num_clusters, affinity='precomputed', assign_labels='kmeans')
    labels = sc.fit_predict(adjacency_matrix)

    # Group nodes based on clusters
    clusters = {}
    for node, label in enumerate(labels):
        clusters.setdefault(label, []).append(node)

    flag = False
    for cluster_id, nodes in clusters.items():
        #print(f"Cluster {cluster_id}: {nodes}")
        if max_qubits and len(nodes) > max_qubits:
            flag=True
    
    if flag:
        return [], []

    clusters = list(clusters.values())
    # Output clusters (partitions)
    #print("Partitions (clusters):")
    #for cluster_id, nodes in clusters.items():
        #print(f"Cluster {cluster_id}: {nodes}")

    # Find the cut edges: edges between nodes in different clusters
    cut_edges = []
    for u, v in graph.edges():
        if labels[u] != labels[v]:
            data = graph.get_edge_data(u,v)
            gates_list = [value['gate'] for value in data.values()]
            cut_edges = cut_edges + gates_list

    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    
    #print(f"Cut edges: {cut_edges}")

    '''import matplotlib.pyplot as plt
    # Visualize the graph with clusters and cut edges
    colors = [f'C{label}' for label in labels]
    plt.figure(figsize=(8, 6))

    # Draw the graph with node coloring
    networkx.draw(graph, with_labels=True, node_color=colors, node_size=700, font_size=10, font_weight='bold', edge_color='gray')

    # Highlight the cut edges in red
    cut_edge_list = [(u, v) for u, v in cut_edges]
    networkx.draw_networkx_edges(graph, pos=networkx.spring_layout(graph), edgelist=cut_edge_list, edge_color='red', width=2.5)

    plt.title("Graph Partitioning and Cut Edges")
    plt.show()'''
    return clusters, cut_edges

def metis_partition(graph, num_clusters, max_nodes_per_cluster, max_cuts):
    # Convert the NetworkX graph to a format suitable for METIS
    print(num_clusters)
    adjacency_list = []
    
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        adjacency_list.append(neighbors)
    _, parts = pymetis.part_graph(adjacency=adjacency_list, nparts=num_clusters)
    # Group nodes based on the partition labels
    clusters = {}
    for node, part in enumerate(parts):
        clusters.setdefault(part, []).append(node)

    # Check for cluster size constraint violations
    for cluster_id, nodes in clusters.items():
        if max_nodes_per_cluster and len(nodes) > max_nodes_per_cluster:
            print(f"Cluster {cluster_id} exceeds the max size.")
            return [], []
            # You can apply further refinements here
    
    clusters = list(clusters.values())
    cut_edges = []
    cut_edges = []
    for u, v in graph.edges():
        if parts[u] != parts[v]:
            data = graph.get_edge_data(u,v)
            gates_list = [value['gate'] for value in data.values()]
            #gates_list = [value[2]['gate'] for value in data.values()]
            cut_edges = cut_edges + gates_list
    if max_cuts and len(cut_edges) > max_cuts:
        return [], []
    elif cut_edges == []:
        return [], []
    return clusters, cut_edges

@task(returns=2)
def find_cut_gate(circuit,max_qubits=None,max_cuts=None, max_components=None, verbose=False):
    #clusters = circuit.nqubits // max_qubits 
    #Gate cutting
    dag = circuit_to_qubit_graph(circuit)
    best_kl_partition, best_kl_cutset = find_best_kernighan_lin_partition(graph=dag, max_qubits=max_qubits, max_cuts=max_cuts)
    best_gn_partition, best_gn_cutset = find_best_girvan_newman_cutset(dag, max_qubits, max_cuts=max_cuts)
    best_sc_partition, best_sc_cutset = find_best_spectral_clustering(dag,clusters = max_components, max_qubits=max_qubits, max_cuts=max_cuts)
    best_met_partition, best_met_cutset = metis_partition(dag, num_clusters=max_components, max_nodes_per_cluster=max_qubits, max_cuts=max_cuts)
    
    # Weights for the objective function
    w1 = 100  # Weight for cut size (minimize this)
    w2 = -1  # Weight for number of components (maximize this)
    w3 = 1  # Weight for balanced graphs (minimize this)
    
    results = []
    final = []

    print(best_kl_cutset)
    print("NUM COMPONENTS: ", len(best_kl_partition))
    qubit_diff = []
    for el in best_kl_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_kl_cutset) + w2 *  len(best_kl_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_kl_cutset)
    
    
    print(best_gn_cutset)
    print("NUM COMPONENTS: ", len(best_gn_partition))
    qubit_diff = []
    for el in best_gn_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_gn_cutset) + w2 *  len(best_gn_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_gn_cutset)
    
    
    print(best_sc_cutset)
    print("NUM COMPONENTS: ", len(best_sc_partition))
    qubit_diff = []
    for el in best_sc_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_sc_cutset) + w2 *  len(best_sc_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_sc_cutset)

    # Perform partitioning with METIS
    print(best_met_cutset)
    print("NUM COMPONENTS: ", len(best_met_partition))
    qubit_diff = []
    for el in best_met_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_met_cutset) + w2 *  len(best_met_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_met_cutset)

    print(results)
    best_score = min(results)
    best_cut = final[results.index(best_score)]
    print(best_cut)

    if verbose:
        import matplotlib.pyplot as plt
        # Plot the graph with partition
        pos = networkx.spring_layout(dag)  # Positions for all nodes
        plt.figure(figsize=(12, 12))
        
        # Get the best partition (assuming the best partition is the one corresponding to the best cut)
        best_partition = None
        if results.index(best_score) == 0:
            best_partition = best_kl_partition
        elif results.index(best_score) == 1:
            best_partition = best_gn_partition
        elif results.index(best_score) == 2:
            best_partition = best_sc_partition
        else:
            best_partition = best_met_partition

        # Assign colors to partitions
        colors = ['orange', 'purple', 'green', 'red', 'blue']
        for i, partition in enumerate(best_partition):
            networkx.draw_networkx_nodes(dag, pos, nodelist=partition, node_color=colors[i % len(colors)], node_size=200, alpha=0.5)
        
        # Draw the edges and labels
        networkx.draw_networkx_edges(dag, pos, alpha=0.5)
        networkx.draw_networkx_labels(dag, pos)

        # Display the plot
        plt.title("Graph Partition Visualization")
        plt.show()
    return best_cut, best_score

@task(returns=2)
def find_cut_wire(circuit,max_qubits=None,max_cuts=None, max_components=None, verbose=False):
    #clusters = circuit.nqubits // max_qubits 
    #Gate cutting
    dag = circuit_to_gate_graph(circuit)
    import matplotlib.pyplot as plt
    best_kl_partition, best_kl_cutset = find_best_kernighan_lin_partition(graph=dag, max_qubits=max_qubits, max_cuts=max_cuts)
    print("KERNIGHAN DONE...")
    best_gn_partition, best_gn_cutset = find_best_girvan_newman_cutset(dag, max_qubits=max_qubits, max_cuts=max_cuts)
    print("GIRVAN NEWMAN DONE...")
    best_sc_partition, best_sc_cutset = find_best_spectral_clustering(dag,clusters = max_components, max_qubits=max_qubits, max_cuts=max_cuts)
    print("SPECTRAL CLUSTERING DONE...")
    best_met_partition, best_met_cutset = metis_partition(dag, num_clusters=max_components, max_nodes_per_cluster=max_qubits, max_cuts=max_cuts)
    print("METIS DONE...")

    # Weights for the objective function
    w1 = 100  # Weight for cut size (minimize this)
    w2 = -1  # Weight for number of components (maximize this)
    w3 = 1  # Weight for balanced graphs (minimize this)
    
    results = []
    final = []

    print(best_kl_cutset)
    print("NUM COMPONENTS: ", len(best_kl_partition))
    qubit_diff = []
    for el in best_kl_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_kl_cutset) + w2 *  len(best_kl_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_kl_cutset)
    
    
    print(best_gn_cutset)
    print("NUM COMPONENTS: ", len(best_gn_partition))
    qubit_diff = []
    for el in best_gn_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_gn_cutset) + w2 *  len(best_gn_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_gn_cutset)
    
    
    print(best_sc_cutset)
    print("NUM COMPONENTS: ", len(best_sc_partition))
    qubit_diff = []
    for el in best_sc_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_sc_cutset) + w2 *  len(best_sc_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_sc_cutset)

    # Perform partitioning with METIS
    print(best_met_cutset)
    print("NUM COMPONENTS: ", len(best_met_partition))
    qubit_diff = []
    for el in best_met_partition:
        qubit_diff.append(len(el))
    qubit_diff.append(float('-inf'))
    print("MAX QUBITS: ", max(qubit_diff))
    score = w1 * len(best_met_cutset) + w2 *  len(best_met_partition) + w3 * max(qubit_diff)
    print("SCORE: ", score)
    results.append(abs(score))
    final.append(best_met_cutset)

    tmp = min(results)
    tmp2= results.index(tmp)
    cut = final[tmp2]
    new_cut = []
    for pair in cut:
        u,v = pair
        _, u_idx = u.split('_')
        u_idx = int(u_idx)
        _, v_idx = v.split('_')
        v_idx = int(v_idx)
        initial_gate = circuit.queue[u_idx-1]
        final_gate = circuit.queue[v_idx-1]
        flag = True
        for idx, gate in enumerate(circuit.queue[u_idx:v_idx-1],start=u_idx):
            if len(gate.qubits) > 2:
                pass
            if set(gate.qubits).intersection(initial_gate.qubits, final_gate.qubits):
                if initial_gate.name.upper() == 'CX':
                    initial_name = 'CNOT'
                else:
                    initial_name = initial_gate.name.upper()
                if gate.name.upper() == 'CX':
                    gate_name = 'CNOT'
                else:
                    gate_name = gate.name.upper()
                new_cut.append((f'{initial_name}_{u_idx}',f'{gate_name}_{idx+1}'))
                flag = False
        if flag:
            new_cut.append(pair)
    print(new_cut)
    return new_cut, tmp


def find_cut(circuit,max_qubits=None,max_cuts=None, max_components=None, wire_cut=True, gate_cut=True, implementation='qdislib', verbose=False):
    if implementation == 'qdislib':
        if max_qubits and not max_components:
            max_components = circuit.nqubits // max_qubits + 1
        elif not max_qubits and not max_components:
            max_components =  2

        if gate_cut:
            cut_gate, score_gate = find_cut_gate(circuit,max_qubits=max_qubits,max_cuts=max_cuts, max_components=max_components,verbose=verbose)
        
        if wire_cut:
            cut_wire,  score_wire = find_cut_wire(circuit,max_qubits=max_qubits,max_cuts=max_cuts, max_components=max_components,verbose=verbose)

        if gate_cut and wire_cut:
            cut_gate = compss_wait_on(cut_gate)
            score_gate = compss_wait_on(score_gate)
            cut_wire = compss_wait_on(cut_wire)
            score_wire = compss_wait_on(score_wire)
        elif gate_cut:
            cut_gate = compss_wait_on(cut_gate)
            score_gate = compss_wait_on(score_gate)
        elif wire_cut:
            cut_wire = compss_wait_on(cut_wire)
            score_wire = compss_wait_on(score_wire)
        else:
            raise NameError

        if wire_cut and gate_cut:
            if score_gate < score_wire:
                return cut_gate
            elif score_gate > score_wire:
                return cut_wire
            else:
                return cut_gate
        elif gate_cut:
            return cut_gate
        elif wire_cut:
            return cut_wire
        else:
            raise TypeError
        
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
            max_qubits_ckt = circuit.num_qubits //2 +1
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

        

