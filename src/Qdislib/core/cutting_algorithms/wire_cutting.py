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
from pycompss.api.parameter import COLLECTION_IN, DICTIONARY_IN
from pycompss.api.parameter import *

import numpy as np
import qibo
from qibo import models, gates, hamiltonians  # , callbacks
import networkx as nx

from Qdislib.api import *
from qiboconnection.connection import ConnectionConfiguration
from qiboconnection.api import API
from Qdislib.api import *
from scipy.optimize import minimize
import time
import itertools as it

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.parameter import *

import igraph as ig
import matplotlib.pyplot as plt
import random

import math

import inspect


def analytical_solution(observables, circuit, verbose=False):
    """Description
    -----------
    Calculate the analytical expected value of a whole circuit.

    Parameters
    ----------
    observables: string
        String containing observables.
    circuit: Circuit
        Circuit object.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    analytical expected value: float
        Analytical expected value.

    Example
    -------
    >>> analytical_value = analytical_solution(observables="ZZZZ", circuit=circuit, verbose=True)
    """

    state = circuit()
    counter = 0
    final = []
    for i in observables:
        if i == "Z":
            final.append(qibo.symbols.Z(counter))
        if i == "X":
            final.append(qibo.symbols.X(counter))
        if i == "Y":
            final.append(qibo.symbols.Y(counter))
        if i == "I":
            final.append(qibo.symbols.I(counter))
        counter = counter + 1

    expectation_value = np.prod(final)
    if verbose:
        print(expectation_value)

    # We convert the expectation value in a symbolic Hamiltonian
    new_expectation_value = hamiltonians.SymbolicHamiltonian(expectation_value)
    # Finally we compute the expectation value
    exp_full_circuit = float(
        new_expectation_value.expectation(
            state.state(numpy=True), normalize=False
        )
    )

    if verbose:
        print(
            "The expectation value of",
            expectation_value,
            "in the entire circuit is ",
            exp_full_circuit,
        )
    return exp_full_circuit



def circuit_to_dag(circuit):
    """
    Convert a Qibo circuit to a DAG where each node stores gate information.
    
    Args:
    - circuit: The Qibo circuit to transform.
    - num_qubits: The number of qubits in the circuit.
    
    Returns:
    - dag: A directed acyclic graph (DAG) with nodes containing gate information.
    """
    # Create a directed graph
    dag = nx.DiGraph()

    # Add gates to the DAG as nodes with unique identifiers
    for gate_idx, gate in enumerate(circuit.queue, start=1):
        # Unique identifier for each gate instance
        gate_name = f"{gate.__class__.__name__}_{gate_idx}"

        # Add the gate to the DAG, including the gate type, qubits, and parameters
        dag.add_node(gate_name, gate=gate.__class__.__name__, qubits=gate.qubits, parameters=gate.parameters)
            
    
        # Connect gates based on qubit dependencies
        for qubit in gate.qubits:
            for pred_gate in reversed(list(dag.nodes)): # Skip the last node since it is the current gate being added
                if dag.nodes[pred_gate].get('qubits') and qubit in dag.nodes[pred_gate]['qubits']:  # Check if the qubit is in the node's qubits
                    if gate_name != pred_gate:
                        dag.add_edge(pred_gate, gate_name, color="blue")
                        break
    
        for qubit in gate.qubits:
            for pred_gate in reversed(list(dag.nodes)):
                if dag.nodes[pred_gate].get('qubits') and qubit in dag.nodes[pred_gate]['qubits']:
                    if gate_name != pred_gate:
                        if not dag.has_edge(pred_gate, gate_name):
                            dag.add_edge(pred_gate,gate_name, color="red")

    return dag

def plot_dag(dag):
    """
    Plot the DAG graph using matplotlib and networkx.
    
    Args:
    - dag: A networkx DiGraph representing the circuit.
    """
    # Set up graph layout
    pos = nx.spring_layout(dag)

    # Draw edges for the first group with blue color
    edges_first_group = [
        (edge[0], edge[1])
        for edge in dag.edges.data("color")
        if edge[2] == "blue"
    ]
    nx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_first_group,
        edge_color="blue",
        width=2.0,
        alpha=0.7,
    )
    
    edges_second_group = [
        (edge[0], edge[1])
        for edge in dag.edges.data("color")
        if edge[2] == "red"
    ]
    nx.draw_networkx_edges(
        dag,
        pos,
        edgelist=edges_second_group,
        edge_color="red",
        width=2.0,
        alpha=0.7,
        style="dotted",
    )

    # Draw nodes
    node_labels = nx.get_node_attributes(dag, 'gate')
    nx.draw_networkx_nodes(dag, pos, node_color="skyblue", node_size=1000)
    nx.draw_networkx_labels(dag, pos, labels=node_labels, font_size=10)

    # Show plot
    plt.title("Circuit DAG")
    plt.show()





@task(returns=1)
def dag_to_circuit(dag, num_qubits):
    """
    Reconstruct a Qibo circuit from a DAG.
    
    Args:
    - dag: A networkx DiGraph representing the circuit.
    - num_qubits: The number of qubits in the original circuit.
    
    Returns:
    - circuit: A Qibo circuit reconstructed from the DAG.
    """
    
    # Create an empty Qibo circuit
    circuit = models.Circuit(num_qubits)
    
    # Traverse the DAG in topological order
    topo_order = list(nx.topological_sort(dag))

    for node in topo_order:
        node_data = dag.nodes[node]
        gate_name = node_data['gate']
        
        # Skip the measurement nodes (we'll handle them separately)
        if gate_name == "Observable I":
            continue
        
        # Get the qubits this gate acts on
        
        qubits = node_data['qubits']
        parameters = node_data['parameters']

        # Get the gate class from the qibo.gates module
        gate_class = getattr(gates, gate_name)

        # Get the signature of the gate's __init__ method
        signature = inspect.signature(gate_class.__init__)

        # Count the number of required positional arguments (excluding 'self')
        param_count = len(signature.parameters) - 1  # exclude 'self'

        # Check if parameters are provided and the gate requires them
        if parameters is not None and param_count > len(qubits):
            # Pass qubits and parameters if the gate requires both
            circuit.add(gate_class(*qubits, parameters))
        else:
            # Otherwise, pass only the qubits
            circuit.add(gate_class(*qubits))
        
        '''if gate_name == 'H':
            circuit.add(gates.H(*qubits))
        elif gate_name == 'CNOT':
            circuit.add(gates.CNOT(*qubits))
        elif gate_name == 'CZ':
            circuit.add(gates.CZ(*qubits))
        elif gate_name == 'X':
            circuit.add(gates.X(*qubits))
        elif gate_name == 'M':
            circuit.add(gates.M(*qubits))
        elif gate_name == 'S':
            circuit.add(gates.S(*qubits))
        elif gate_name == 'T':
            circuit.add(gates.T(*qubits))
        elif gate_name == 'SDG':
            circuit.add(gates.SDG(*qubits))
        elif gate_name == 'RZ':
            circuit.add(gates.RZ(*qubits, parameters))
        elif gate_name == 'RY':
            circuit.add(gates.RY(*qubits, parameters))
        elif gate_name == 'RX':
            circuit.add(gates.RX(*qubits, parameters))
        elif gate_name == 'Z':
            circuit.add(gates.Z(*qubits))

        else:
            raise ValueError(f"Unsupported gate type: {gate_name}")'''

    # Optionally handle measurements, assuming all qubits are measured at the end
    obs_I = []
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data['gate'] == "Observable I":
            obs_I.append(node_data['qubits'][0])
            dag.remove_node(node)
            #print(obs_I)
            #circuit.add(gates.M(node_data['qubit']))

    if obs_I:
        return [circuit, obs_I]

    return [circuit, None]


def wire_cutting(rand_qc,cut,sync=True,gate_cutting=False):
    if type(rand_qc) == models.Circuit:
        dag = circuit_to_dag(rand_qc)
         
    else:
        dag = rand_qc
    
    
    if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [dag.subgraph(c).copy() for c in nx.connected_components(dag.to_undirected())]
        results = []
        for s in S:
            num_qubits = max_qubit(s)
            tmp_cuts = []
            for c in cut:
                if s.has_edge(*c):
                    tmp_cuts.append(c)
            if tmp_cuts:
                #print(dag_to_circuit(s,5)[0].draw())
                graphs = generate_wire_cutting(s, tmp_cuts , num_qubits=num_qubits)
                print("GRAPHS ",graphs)
                results.append(sum(graphs)) #1/(2**len(tmp_cuts))*sum(graphs)
            else:
                s_new, highest_qubit = update_qubits(s)
                subcirc = dag_to_circuit(s_new,highest_qubit)
                #print(subcirc.draw())
                expected_value = execute_subcircuits(subcirc)
                results.append(expected_value)       
                print("EV ",expected_value)
        
        if sync:
            results = compss_wait_on(results)

        '''if gate_cutting:
            final_recons = 1/2*sum(results)
        else:'''
        print(results)
        final_recons = 1/(2**len(cut))*math.prod(results)
        return final_recons
    
    else:
        num_qubits = max_qubit(dag)
        if cut:
            results = generate_wire_cutting(dag, cut , num_qubits=num_qubits)

            if sync:
                results = compss_wait_on(results)
            #print(results)
            final_recons = 1/(2**len(cut))*sum(results)
        else:
            s_new, highest_qubit = update_qubits(dag)
            subcirc = dag_to_circuit(s_new,highest_qubit)
            final_recons = execute_subcircuits(subcirc)
        return final_recons

def generate_wire_cutting(dag, edges_to_replace, num_qubits):
    """
    Replace a specific edge in the DAG with a source and end node.
    
    Args:
    - dag: The directed acyclic graph (DAG) to modify.
    - edge_to_replace: The edge to remove (tuple of nodes).
    - num_qubits: The current number of qubits in the circuit.
    
    Returns:
    - updated_dag: The modified DAG with new source and end nodes.
    """
    
    reconstruction = []

    for index, edge_to_replace in enumerate(edges_to_replace, start=1):
    
        # Extract the nodes of the edge to be replaced
        source, target = edge_to_replace
        
        # Remove the original edge
        dag.remove_edge(source, target)

        source_gate_info = dag.nodes[source]

        target_gate_info = dag.nodes[target]

        common_qubit = list(set(target_gate_info.get('qubits')).intersection(set(source_gate_info.get('qubits'))))

        successors = []
        # Iterate over all nodes in the graph
        for node in dag.nodes:
            if dag.has_edge(target, node):
                successors.append(node)
    
        # Include the target node itself
        nodes = [target] + successors

        for successor in nodes:
            qubits = dag.nodes[successor].get('qubits')
            for qubit in qubits:

                if common_qubit[0] is qubit:
                    temp_list = list(dag.nodes[successor].get('qubits'))
            
                    # Replace the common element with the new value
                    for i in range(len(temp_list)):

                        if temp_list[i] == common_qubit[0]:
                            temp_list[i] = num_qubits+index
            
                    updated_tuple = tuple(temp_list)
                    dag.nodes[successor]['qubits'] = updated_tuple

        dag.add_node(f"O_{index}", gate='S', qubits=common_qubit, parameters=())
        
        # Add the new end node with the same properties as the target node
        dag.add_node(f"PS_{index}", gate='T', qubits=(num_qubits+index,), parameters=())

        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")
    
        copy_dag = dag.copy()
        red_edges = []
        for ed in dag.edges:
            if dag.get_edge_data(ed[0],ed[1])["color"] == "red":
                red_edges.append(ed)
        
        copy_dag.remove_edges_from(red_edges)

    #print(dag_to_circuit(dag,6)[0].draw())

    graphs = []
    for i in range(8**len(edges_to_replace)):
        graphs.append(dag.copy())

    for index, graph in enumerate(graphs, start=0):
        copy_graph = graph.copy()

        copy_graph = remove_red_edges(copy_graph)

        num_components = nx.number_connected_components(copy_graph.to_undirected())

        graph_components = []
        for i in range(num_components):
            graph_components.append(nx.DiGraph().copy()) 

        graph = generate_subcircuits_wire_cutting(copy_graph, num_qubits+len(edges_to_replace),index, edges_to_replace, graph_components)

        exp_value = []
        for s in graph_components:
            s_new, highest_qubit = update_qubits(s)
            #print(s_new.nodes(data=True))
            subcirc = dag_to_circuit(s_new,highest_qubit)
            #print(subcirc.draw())
            expected_value = execute_subcircuits(subcirc)
            exp_value.append(expected_value)

        #print(exp_value)
        exp_value = change_sign(exp_value, index)
        #print(exp_value)
        reconstruction.append(exp_value)
    return reconstruction

def max_qubit(graph):
    # Initialize a variable to keep track of the highest Qubits value
    max_qubits = float('-inf')  # Start with the lowest possible number
    max_node = None  # Store the node with the highest qubit

    # Iterate over the nodes and check their 'Qubits' attribute
    for node, data in graph.nodes(data=True):
        qubits = data.get('qubits', 0)  # Default to 0 if 'Qubits' is not present
        for qubit in qubits:
            if qubit > max_qubits:
                max_qubits = qubit
                max_node = node
    return max_qubits


def connected_components(G):
    """Generate connected components.

    Parameters
    ----------
    G : NetworkX graph
       An undirected graph

    Returns
    -------
    comp : generator of sets
       A generator of sets of nodes, one for each component of G.

    Raises
    ------
    NetworkXNotImplemented
        If G is directed.

    Examples
    --------
    Generate a sorted list of connected components, largest first.

    >>> G = nx.path_graph(4)
    >>> nx.add_path(G, [10, 11, 12])
    >>> [len(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    [4, 3]

    If you only want the largest connected component, it's more
    efficient to use max instead of sort.

    >>> largest_cc = max(nx.connected_components(G), key=len)

    To create the induced subgraph of each component use:

    >>> S = [G.subgraph(c).copy() for c in nx.connected_components(G)]

    See Also
    --------
    strongly_connected_components
    weakly_connected_components

    Notes
    -----
    For undirected graphs only.

    """
    seen = set()
    for v in G:
        if v not in seen:
            c = plain_bfs(G, v)
            seen.update(c)
            yield c

@task(returns=1)
def plain_bfs(G, source):
    """A fast BFS node generator"""
    adj = G._adj
    n = len(adj)
    seen = {source}
    nextlevel = [source]
    while nextlevel:
        thislevel = nextlevel
        nextlevel = []
        for v in thislevel:
            for w in adj[v]:
                if w not in seen:
                    seen.add(w)
                    nextlevel.append(w)
            if len(seen) == n:
                return seen
    return seen


@task(returns=1)
def compute_exp_value(exp_value, expected_value):
    return exp_value * expected_value


@task(returns=1)
def get_subgraph(graph, c):
    return graph.subgraph(c).copy()

@task(returns=2, s=INOUT)
def update_qubits(s):
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)


    for node, data in s.nodes(data=True):
        new_tuple = ()
        for qubit in s.nodes[node]["qubits"]:
            len_missing = count_missing_up_to(my_set, qubit)
            new_qubit = qubit - len_missing
            new_tuple = new_tuple + (new_qubit,)
        s.nodes[node]["qubits"] = new_tuple

    highest_qubit = max(my_set)+1 - count_missing_up_to(my_set, max(my_set))
    return s, highest_qubit

def remove_red_edges(graph):
    copy_dag = graph.copy()
    red_edges = []
    
    for ed in copy_dag.edges:
        if copy_dag.get_edge_data(ed[0],ed[1])["color"] == "red":
            red_edges.append(ed)

    copy_dag.remove_edges_from(red_edges)
    return copy_dag


@task(returns=1)
def generate_set(s):
    my_set = set()
    for node, data in s.nodes(data=True):
        for qubit in s.nodes[node]["qubits"]:
            my_set.add(qubit)
    return my_set


def count_missing_up_to(nums, max_num):
    # Create a set of all numbers from 0 to max_num
    full_set = set(range(max_num + 1))
    
    # Subtract the given set from the full set to get the missing numbers
    missing_numbers = full_set - nums
    
    # Return the count of missing numbers
    return len(missing_numbers)


@task(returns=1, graph_components=COLLECTION_OUT)
def generate_subcircuits_wire_cutting(updated_dag, num_qubits, idx, edges_to_replace, graph_components):

    base8_rep = oct(idx)[2:]
    base8_rep = base8_rep.zfill(len(edges_to_replace))
    list_substitutions = list(map(int, base8_rep))

    for idx2, index in enumerate(list_substitutions, start=0):
        idx2 = idx2+1
        
        # I 0 
        if index == 0:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'Observable I'
            #updated_dag.remove_node(f'O_{idx2}')
            updated_dag.remove_node(f'PS_{idx2}')

        # I 1
        elif index == 1:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'Observable I'
            #updated_dag.remove_node(f'O_{idx2}')
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'X'

        # X +
        elif index == 2:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'H'
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'H'
        
        # X -
        elif index == 3:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'H'
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'H'
            updated_dag.add_node(f'PS2_{idx2}', gate='X', qubits=updated_dag.nodes[f'PS_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'PS2_{idx2}', f'PS_{idx2}', color="blue")

        # Y +i
        elif index == 4:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'SDG'
            updated_dag.add_node(f"O2_{idx2}", gate='H', qubits=updated_dag.nodes[f'O_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'O_{idx2}', f'O2_{idx2}', color="blue")
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'S'
            updated_dag.add_node(f"PS2_{idx2}", gate='H', qubits=updated_dag.nodes[f'PS_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'PS2_{idx2}', f'PS_{idx2}', color="blue")

        # Y -i
        elif index == 5:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'SDG'
            updated_dag.add_node(f"O2_{idx2}", gate='H', qubits=updated_dag.nodes[f'O_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'O_{idx2}', f'O2_{idx2}', color="blue")
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'S'
            updated_dag.add_node(f"PS2_{idx2}", gate='H', qubits=updated_dag.nodes[f'PS_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'PS2_{idx2}', f'PS_{idx2}',color="blue")
            updated_dag.add_node(f"PS3_{idx2}", gate='X', qubits=updated_dag.nodes[f'PS_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'PS3_{idx2}', f'PS2_{idx2}', color="blue")

        # Z 0
        elif index == 6:
            updated_dag.remove_node(f'O_{idx2}')
            updated_dag.remove_node(f'PS_{idx2}')

        # Z 1
        elif index == 7:
            updated_dag.remove_node(f'O_{idx2}')
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'X'

        else:
            print("ERROR")
            raise ValueError

    updated_dag = remove_red_edges(updated_dag)
    for i, c in enumerate(nx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")
    
    return updated_dag

@task(returns=1)
def execute_subcircuits(subcirc):
    tmp = subcirc[1]
    subcirc = subcirc[0]
    if tmp:
        obs_I=tmp
    else:
        obs_I = None
    observables = ['Z']*subcirc.nqubits

    if obs_I:
        for element in obs_I:
            observables[element] = 'I'

    observables = ''.join(observables)
    print(observables)
    
    qibo.set_backend("numpy")
    shots=50000
    subcirc.add(gates.M(*range(subcirc.nqubits))) #
    result = subcirc(nshots=shots)
    freq = dict(result.frequencies(binary=True))

    expectation_value = 0
    for key, value in freq.items():
        contribution = 1
        for bit, obs in zip(key, observables):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")
        
        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / shots)
    #print(expectation_value)
    return expectation_value

@task(returns=1, expectation_value=COLLECTION_IN)
def change_sign(expectation_value, index):
    expectation_value = math.prod(expectation_value)
    number = index

    change_sign = False

    while number != 0:
        digit = number % 8  # Get the last digit
        if digit in {3, 5, 7}:  # Check if the digit is 3, 5, or 7
            change_sign = not change_sign  # Flip the sign change flag
        number //= 8 # Move to the next digit
    #print(change_sign)

    # If change_sign is True, we flip the sign of the original number
    if change_sign:
        return -expectation_value
    else:
        return expectation_value

def random_circuit(qubits, gate_max, num_cz, p):
    if p == None:
        graph = ig.Graph.Erdos_Renyi(n=qubits, m=num_cz, directed=False, loops=False)
    elif num_cz == None:
        graph = ig.Graph.Erdos_Renyi(n=qubits, p=p, directed=False, loops=False)
    else:
        print(
            "Error: only the number of edges or the probability of adding an edge must be specified"
        )

    # adj_mat=graph.get_adjacency()
    # fig, ax = plt.subplots()
    # ig.plot(graph, target=ax, vertex_label=range(qubits))
    # graph.degree()

    edge_list = graph.get_edgelist()

    gates_pull = [gates.X, gates.H, gates.S, gates.T]  # pull of single-qubit gates
    circuit = models.Circuit(qubits)
    for i in range(len(edge_list)):
        rand_tmp = random.randint(
            0, gate_max
        )  # number of single-qubit gates between CZ
        for j in range(rand_tmp):
            sel_gate = random.choice(gates_pull)  # gate selected from the pull
            sel_qubit = random.randint(
                0, qubits - 1
            )  # qubit selected to apply the gate
            circuit.add(sel_gate(sel_qubit))

        circuit.add(
            gates.CZ(edge_list[i][0], edge_list[i][1])
        )  # 2-qubit gate from graph

    return circuit

def generate_combinations(n, gate_type):
    """
    Generate combinations for gate cutting
    depending on the type of the gate and the
    number of gates being cut.

    :param n: int.
    :param gate_type: string
    :return: all_combinations.
    """
    objects = []
    if gate_type == 'CZ':
        objects = [('S', 'S'), ('SDG', 'SDG')]
    elif gate_type == 'CNOT':
        objects = [
            ('RZ','RX'),
            ('RZ','Z','RX', 'X'),
        ]
    all_combinations = list(it.product(objects, repeat=n))
    return all_combinations

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

def generate_gate_cutting(dag, target_node, index, cuts_wc):
    dag_wc = dag.copy()

    target_qubits = dag.nodes[target_node]['qubits']

    # Find predecessor node with qubit 1
    pred_0 = find_nodes_with_qubit(dag, target_node, qubit=target_qubits[0], direction='predecessor')

    # Find predecessor node with qubit 2
    pred_1 = find_nodes_with_qubit(dag, target_node, qubit=target_qubits[1], direction='predecessor')

    # Find successor node with qubit 1
    succ_0 = find_nodes_with_qubit(dag, target_node, qubit=target_qubits[0], direction='successor')

    # Find successor node with qubit 2
    succ_1 = find_nodes_with_qubit(dag, target_node, qubit=target_qubits[1], direction='successor')

    # Output the results
    print(f"Predecessor nodes with qubit {target_qubits[0]}: {pred_0}")
    print(f"Predecessor nodes with qubit {target_qubits[1]}: {pred_1}")
    print(f"Successor nodes with qubit {target_qubits[0]}: {succ_0}")
    print(f"Successor nodes with qubit {target_qubits[1]}: {succ_1}")

    dag.remove_node(target_node)
    dag_wc.remove_node(target_node)

    dag.add_node(f"GCA_{index}", gate='S', qubits=(target_qubits[0],), parameters=())
    dag.add_node(f"GCB_{index}", gate='T', qubits=(target_qubits[1],), parameters=())

    if pred_0:
        dag.add_edge(pred_0[0],f"GCA_{index}", color="blue")

    if succ_0:
        dag.add_edge(f"GCA_{index}", succ_0[0], color="blue")

    if pred_1:
        dag.add_edge(pred_1[0],f"GCB_{index}", color="blue")
    
    if succ_1:
        dag.add_edge(f"GCB_{index}", succ_1[0], color="blue")

    if pred_0 and succ_0:
        dag_wc.add_edge(pred_0[0],succ_0[0], color="blue")
        cuts_wc.append((pred_0[0],succ_0[0]))

    if pred_1 and succ_1:
        dag_wc.add_edge(pred_1[0],succ_1[0], color="blue")
        cuts_wc.append((pred_1[0],succ_1[0]))

    return dag, dag_wc, cuts_wc

@task(returns=1, graph_components=COLLECTION_OUT)
def generate_subcircuits_gate_cutting(index, target_nodes,graphs, all_combinations, graph_components):
    idx = index
    combination = all_combinations[index]
    for idx2, nodes in enumerate(target_nodes):
        graphs[idx].nodes[f'GCA_{idx2+1}']['gate'] = combination[idx2][0]
        graphs[idx].nodes[f'GCA_{idx2+1}']['parameters'] = np.pi/2
        graphs[idx].nodes[f'GCB_{idx2+1}']['gate'] = combination[idx2][1]
        graphs[idx].nodes[f'GCB_{idx2+1}']['parameters'] = np.pi/2
    
    updated_dag = remove_red_edges(graphs[idx])
    for i, c in enumerate(nx.connected_components(updated_dag.to_undirected())):
        new_subgraph = updated_dag.subgraph(c).copy()
        graph_components[i].add_nodes_from(new_subgraph.nodes(data=True))
        graph_components[i].add_edges_from(new_subgraph.edges(data=True), color="blue")
    return updated_dag

@task(returns=1, exp_value=COLLECTION_IN)
def math_prod(exp_value):
    expectation_value = math.prod(exp_value)
    return expectation_value

from collections import defaultdict
def draw_to_circuit(text_draw, parameters=None):
    split = text_draw.splitlines()
    #print(split) 
    print(split)
    qubits_lst = []
    split = [element for element in split if element.strip()]
    split = [element for element in split if element != ""]
    print(split)
    for line in split:
        index = line.index('─')
        qubits_lst.append(line[index:])
        

    list_multiple_gates = defaultdict(list)
    # Now we will process each line to identify multi-qubit gates
    for idx, qubit_line in enumerate(qubits_lst):
        qubit_number = idx  # Line number corresponds to the qubit (q0 is index 0)
        qubit_state = list(qubit_line)
        
        # Boolean to track if we are inside a multi-qubit gate
        for i, symbol in enumerate(qubit_state):
            if symbol == 'o':
                index = i
                for idx2, qubit in enumerate(qubits_lst[idx+1:]):
                    if list(qubit)[index] != '|':
                        name = list(qubit)[index]
                        if name == 'Z':
                            name = 'CZ'
                        elif name == 'X':
                            name = 'CNOT'
                        list_multiple_gates[idx].append((name, (idx,idx2+idx+1)))
                        qubits_lst[idx2+idx+1] = qubits_lst[idx2+idx+1][:index] + '─' + qubits_lst[idx2+idx+1][index+1:] 
                        break

    circuit = models.Circuit(len(qubits_lst))
    
    num_steps = len(list(qubits_lst[0]))  # Total number of time steps (columns)
    
    for step in range(num_steps):
        print(qubits_lst)
        saved_qubit = []
        for idx, qubit_line in enumerate(qubits_lst):
            qubit_state = list(qubit_line)
            parameter_tracker = 0

            char = qubit_state[step]
            #for idx2, char in enumerate(qubit_state):
            if char != '─' and char != '|':
                if char != 'o':
                    if qubit_state[step+1] == '─' and qubit_state[step-1] == '─':
                        tmp = char
                        #print("Add gate: ", tmp, " qubit ", (idx,))
                        #circuit.add(getattr(gates, tmp)(idx))
                        print(tmp)
                        gate_name = tmp
                        qubits = idx
                        
                        # Get the gate class from the qibo.gates module
                        gate_class = getattr(gates, gate_name)

                        # Get the signature of the gate's __init__ method
                        signature = inspect.signature(gate_class.__init__)

                        # Count the number of required positional arguments (excluding 'self')
                        param_count = len(signature.parameters) - 1  # exclude 'self'

                        # Check if parameters are provided and the gate requires them
                        if parameters is not None and param_count > 1:
                            param = parameters[idx][parameter_tracker][1]
                            # Pass qubits and parameters if the gate requires both
                            circuit.add(gate_class(qubits, param))
                            parameter_tracker += 1
                        else:
                            # Otherwise, pass only the qubits
                            circuit.add(gate_class(qubits))


                        
                    elif qubit_state[step-1] == '─' and qubit_state[step+1] != '─':
                        tmp = ''
                        print(qubit_state)
                        print(qubit_state[step+1])
                        print(range(step,num_steps))
                        for i in range(step,num_steps):
                            print(qubit_state[i])
                            if qubit_state[i+1] == '─':
                                print("HEY")
                                tmp = tmp + qubit_state[i]

                                gate_name = tmp
                                qubits = idx

                                print(tmp)

                                # Get the gate class from the qibo.gates module
                                gate_class = getattr(gates, gate_name)

                                # Get the signature of the gate's __init__ method
                                signature = inspect.signature(gate_class.__init__)

                                # Count the number of required positional arguments (excluding 'self')
                                param_count = len(signature.parameters) - 1  # exclude 'self'

                                # Check if parameters are provided and the gate requires them
                                if parameters is not None and param_count > 1:
                                    print("HEY2")
                                    param = parameters[idx][parameter_tracker][1]
                                    # Pass qubits and parameters if the gate requires both
                                    circuit.add(gate_class(qubits, param))
                                    parameter_tracker += 1
                                    break
                                else:
                                    print("HEY3")
                                    # Otherwise, pass only the qubits
                                    print(gate_class)

                                    circuit.add(gate_class(qubits))
                                    break

                            else:
                                tmp = tmp + qubit_state[i]
                        

                
                elif char == 'o':
                    saved_qubit.append(idx)


        for idx in saved_qubit:
        #if list_multiple_gates[idx]:
            print("Add gate: ", list_multiple_gates[idx][0][0] ," qubit ", list_multiple_gates[idx][0][1])
            circuit.add(getattr(gates, list_multiple_gates[idx][0][0])(*list_multiple_gates[idx][0][1]))        
            list_multiple_gates[idx].remove(list_multiple_gates[idx][0])

    print(circuit.draw())
    return circuit


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