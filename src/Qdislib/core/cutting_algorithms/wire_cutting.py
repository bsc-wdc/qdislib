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

import qibo
from qibo import models, gates
import networkx as nx

from Qdislib.api import *

import math

from Qdislib.utils.graph import circuit_to_dag, dag_to_circuit, max_qubit, update_qubits, remove_red_edges

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
                graphs = sum_results(graphs)
                results.append(graphs) #1/(2**len(tmp_cuts))*sum(graphs)
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

@task(returns=1, lst=COLLECTION_IN)
def sum_results(lst):
    return sum(lst)

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
