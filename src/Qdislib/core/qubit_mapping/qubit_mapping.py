import networkx as nx
import matplotlib.pyplot as plt
from qibo import models

def qubit_arch(circuit):
    G1 = nx.Graph()
    lst = []
    for gate in circuit.queue:
        if len(gate.qubits) > 1:
            lst.append(gate)
            G1.add_edge(gate.qubits[0], gate.qubits[1])
            
    pos = nx.spring_layout(G1)  # positions for all nodes
    nx.draw(G1, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', font_color='black', font_size=10)
    plt.show()
    return G1


def subgraph_matcher(architecture,circuit_graph):
    # Create an example target graph
    target_graph = architecture
    
    pos = nx.spring_layout(target_graph)  # positions for all nodes
    nx.draw(target_graph, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', font_color='black', font_size=10)
    plt.show()
    
    # Create an example pattern graph
    pattern_graph = circuit_graph
    
    pos = nx.spring_layout(pattern_graph)  # positions for all nodes
    nx.draw(pattern_graph, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', font_color='black', font_size=10)
    plt.show()
    
    # Initialize GraphMatcher with the target and pattern graphs
    matcher = nx.algorithms.isomorphism.GraphMatcher(target_graph, pattern_graph)
    
    # Check if the pattern graph is a subgraph of the target graph
    is_subgraph = matcher.subgraph_is_isomorphic()
    
    # If there is a subgraph isomorphism, print the mapping
    if is_subgraph:
        mapping = matcher.mapping
        print("Subgraph found! Node mapping:", mapping)
    else:
        print("No subgraph isomorphism found.")
    
    # You can also get all subgraph isomorphisms
    all_subgraph_isomorphisms = list(matcher.subgraph_isomorphisms_iter())
    print("All subgraph isomorphisms:", all_subgraph_isomorphisms)
    return all_subgraph_isomorphisms


def mapping_order(target, pattern, order):
    # Create an example target graph
    target_graph = target
    
    # Create an example pattern graph
    pattern_graph = pattern
    
    # Specify the desired order of nodes
    desired_order = order
    
    # Initialize GraphMatcher with the target and pattern graphs
    matcher = nx.algorithms.isomorphism.GraphMatcher(target_graph, pattern_graph)
    
    # Get all subgraph isomorphisms
    all_subgraph_isomorphisms = list(matcher.subgraph_isomorphisms_iter())
    
    results = []
    for element in all_subgraph_isomorphisms:
        keys_list = list(element.keys())
        result = 0
        for index, x in enumerate(desired_order):
            if x in element:
                result = result + (index+1)
        results.append([element, result])
    print(results)
    min_second_element = min(results, key=lambda x: x[1])
    best_architecture = min_second_element[0]
    print('Best architecture mapping: ',best_architecture)
    return best_architecture

def rename_qubits(subcirc, qubit_middle, best_arch):
    target_qubit = best_arch['B']
    new_circuit = models.Circuit(subcirc.nqubits)
    new_circuit.queue = subcirc.queue
    for index, element in enumerate(subcirc.queue):
        qubit = element.qubits
        if len(qubit) < 2:
            qubit = qubit[0]
            if qubit == target_qubit:
                element._set_target_qubits((qubit_middle,))
            elif qubit == qubit_middle:
                element._set_target_qubits((target_qubit,))
        else:
            qubit_0 = element.qubits[0]
            qubit_1 = element.qubits[1]
            if qubit_0 == target_qubit:
                element._set_control_qubits((qubit_middle,))
            elif qubit_1 == target_qubit:
                element._set_target_qubits((qubit_middle,))
            if qubit_0 == qubit_middle:
                element._set_control_qubits((target_qubit,))
            elif qubit_1 == qubit_middle:   
                element._set_target_qubits((target_qubit,))
    return new_circuit