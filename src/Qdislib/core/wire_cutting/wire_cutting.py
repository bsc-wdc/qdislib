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

import numpy as np
import qibo
from qibo import models, gates, hamiltonians  # , callbacks
from qibo.symbols import Z, X, Y, I

import copy
import networkx as nx
import matplotlib.pyplot as plt

# for connecting with the quantum computer
# from qiboconnection.connection import ConnectionConfiguration
# from qiboconnection.api import API

from Qdislib.core.wire_cutting.pycompss_functions import *
from Qdislib.utils.graph import *

from pycompss.api.api import compss_wait_on


def wire_cutting(
    observables, circuit, gate_tuple, shots=30000, draw=False, verbose=False
):
    """Implements the algorithm of circuit cutting in one function.
    Cuts the circuit in 2 after the specified gate, and calculates the
    expected value of the reconstruction.

     :param ciruit: Circuit.
     :param gate_tuple: int tuple.
     :param observable: str.
     :param shots: int.
     :param verbose: bool.
     :return: reconstruction value.
    """

    # Show qibo version
    if verbose:
        print(f"qibo version: {qibo.__version__}")

    qubit, list_subcircuits, lst_observables = split(observables, circuit, gate_tuple,draw,verbose)      
    reconstruction = simulation(lst_observables,qubit,list_subcircuits[0], list_subcircuits[1],shots)
    return reconstruction


def split(observables, circuit, edge_remove, draw=False, verbose=False):

    """Splits a circuit in two subcircuits. Cuts 
       after the gate we pass as the parameter.

        :param ciruit: Circuit.
        :param numgates: int.
        :param draw: bool.
        :param verbose: bool.
        :return: circuit1, circuit2.
        """

    circuit = circuit.copy()
    qubits_first_gate = circuit.queue[edge_remove[0]-1].qubits
    qubits_second_gate = circuit.queue[edge_remove[1]-1].qubits
    first_gate = circuit.queue[edge_remove[0]-1]
    second_gate = circuit.queue[edge_remove[1]-1]
    if verbose: print(first_gate, second_gate)
    #print(circuit.queue[edge_remove[0]-1].qubits > circuit.queue[edge_remove[1]-1].qubits)
    qubit = list(set(circuit.queue[edge_remove[0]-1].qubits).intersection(set(circuit.queue[edge_remove[1]-1].qubits)))
    qubit = qubit[0]
    pos1 = qubits_first_gate.index(qubit)
    pos2 = qubits_second_gate.index(qubit)

    observable_dict = {}
    for num_qubit in range(0, circuit.nqubits):
        observable_dict[num_qubit] = observables[num_qubit]
    if verbose: print(observable_dict)

    #convert to DAG and DIGRPAPH 
    digraph = nx.Graph()
    dag = DAGgraph()

    build_dag(circuit,dag)
    create_graph(dag,digraph)
    
    red_edges = [(edge[0], edge[1]) for edge in digraph.edges.data('color') if edge[2] == 'red']
    digraph.remove_edges_from(red_edges)
    
    digraph.remove_edge(edge_remove[0],edge_remove[1])

    subgraphs = list(nx.connected_components(digraph))
    if len(subgraphs) > 2:
        print("MORE THAN 2 SUBGRAPH")
        return None, None
    #print(subgraphs)

    #print_graph(digraph)
    list_subcircuits = []
    non_empty_list = []
    for subgraph in subgraphs:
        subgraph = sorted(subgraph)
        selected_elements = [dag.nodes[i-1] for i in subgraph]
        #circuit_copy = copy.deepcopy(new_circuit)
        
    #remove specific qubit
        

        circuit_copy = models.Circuit(circuit.nqubits)
        circuit_copy.add(selected_elements)

        non_empty_qubits = del_empty_qubits(circuit_copy)
        non_empty_qubits.sort()
        #print("Non empty qubits!: ",non_empty_qubits)
        difference_list = [value - index for index, value in enumerate(non_empty_qubits)]
        #print("Different list: ",difference_list)
        non_empty_list.append(difference_list)
        subtracted_list = [x - y for x, y in zip(non_empty_qubits, difference_list)]
        #print("Substracted list: ", subtracted_list)

        for gate in circuit_copy.queue:
            if len(gate.qubits) > 1:
                control = subtracted_list[non_empty_qubits.index(gate.qubits[0])]
                gate._set_control_qubits((control,))
                target = subtracted_list[non_empty_qubits.index(gate.qubits[1])]
                gate._set_target_qubits((target,))
            else:
                target = subtracted_list[non_empty_qubits.index(gate.qubits[0])]
                gate._set_target_qubits((target,))
        circuit_copy.nqubits = len(non_empty_qubits)
        circuit_copy.queue.nmeasurements = 0
        list_subcircuits.append(circuit_copy)
    
    if verbose:print(qubits_first_gate, qubits_second_gate)
    if verbose:print(non_empty_list)
    if verbose:print(non_empty_qubits)
    if verbose:print(difference_list)
    if verbose:print(subtracted_list)
    new_qubit = [first_gate.qubits[pos1], second_gate.qubits[pos2]]
    if verbose:print(new_qubit)
    if verbose:print(first_gate.qubits, second_gate.qubits)

    obs1 = {}
    obs2 = {}

    for index, x in enumerate(non_empty_list[0]):
        if first_gate in list_subcircuits[0].queue:
            if index == new_qubit[0]:
                obs1[index] = '-'
            else:
                obs1[index] = observable_dict[index + x]
        else:
            obs1[index] = observable_dict[index + x]

    for index, x in enumerate(non_empty_list[1]):
        if first_gate not in list_subcircuits[0].queue:
            if index == new_qubit[0]:
                obs2[index] = '-'
            else:
                obs2[index] = observable_dict[index + x]
        else:
            obs2[index] = observable_dict[index + x]

    if verbose:print(obs1)
    if verbose:print(obs2)
    lst_observables = [obs1,obs2]

    if first_gate not in list_subcircuits[0].queue:
        list_subcircuits.reverse()
        lst_observables.reverse()
    if draw:
        for circ in list_subcircuits:
            print(circ.draw())
            print("\n")
    return new_qubit, list_subcircuits, lst_observables


def simulation(lst_observables,qubit, circuit_1, circuit_2=None, shots=30000, verbose=False):

    """Performs the execution of a cirucuit to calculate the expected value.
       It accepts one or two circuits. With 1 circuit it calculates the expected
       value straight forward, with 2 it performs a reeconstruction in order 
       to provie the expected value.

        :param observables: string.
        :param ciruit1: Circuit.
        :param ciruit2: Circuit.
        :param shots: int.
        :param verbose: bool.
        :return: reeconstruction value.
        """
    
    if verbose: print(f"qibo version: {qibo.__version__}")

    if circuit_2 != None:
        # SIMULATION

        basis = ["X", "Y", "Z", "I"]
        states = ["0", "1", "+", "+i"]

        observables_1 = lst_observables[0]
        observables_2 = lst_observables[1]
        #num_z_sub1 = observables[:(sub_circuit_1_dimension - 1)]
        #num_z_sub2 = observables[(sub_circuit_1_dimension - 1):]

        # first subcircuit:
        exp_value_1 = {}
        for b in basis:
            if verbose: print("Basis: ", b)
            copy_circuit1 = models.Circuit(circuit_1.nqubits)
            copy_circuit = circuit_1.copy(True)
            copy_circuit1.queue = copy_circuit.queue
            circuit1 = first_subcircuit_basis(copy_circuit1, b,qubit[0])
            result = circuit1.execute_COMPSs(nshots=shots)
            # obtain probability distribution
            freq = result.frequencies_COMPSs(binary=True)
            # we call the function that computes the e
            #new_obs = num_z_sub1[:qubit[0]] + b + num_z_sub1[qubit[0]:]
            #print(new_obs)
            for key, value in observables_1.items():
                if value == '-':
                    observables_1[key] = b
            new_obs = ''.join([value for key, value in sorted(observables_1.items())])
            if verbose: print("OBSERVABLES: ", new_obs)
            exp_value_1[b] = compute_expectation_value(freq, new_obs, shots=shots)

        # second subcircuit:
        exp_value_2 = {}
        if verbose: print("SECOND ONE")
        for s in states:
            if verbose:print("States: ", s)
            copy_circuit2 = models.Circuit(circuit_2.nqubits)
            copy_circuit = circuit_2.copy(True)
            copy_circuit2.queue = copy_circuit.queue
            circuit2 = second_subcircuit_states(copy_circuit2, s, qubit[1])
            result = circuit2.execute_COMPSs(nshots=shots)
            # obtain probability distribution
            freq = result.frequencies_COMPSs(binary=True)
            new_obs2 = ''.join([value for key, value in sorted(observables_2.items())])
            # we call the function that computes the expectation value
            if verbose: print("OBSERVABLES: ", new_obs2)
            exp_value_2[s] = compute_expectation_value(freq, new_obs2, shots=shots)

        exp_value_1 = compss_wait_on(exp_value_1)
        exp_value_2 = compss_wait_on(exp_value_2)

        reconstruction = (
            1
            / 2
            * (
                (exp_value_1["I"] + exp_value_1["Z"]) * exp_value_2["0"]
                + (exp_value_1["I"] - exp_value_1["Z"]) * exp_value_2["1"]
                + exp_value_1["X"] * (2 * exp_value_2["+"] - exp_value_2["0"] - exp_value_2["1"])
                + exp_value_1["Y"] * (2 * exp_value_2["+i"] - exp_value_2["0"] - exp_value_2["1"])
            )
        )

        print("Expectation value after circuit cutting and reconstruction:", reconstruction)
        return reconstruction
    else:
        circuit_1.add(gates.M(*range(circuit_1.nqubits)))
        result = circuit_1(nshots=shots)
        freq = dict(result.frequencies(binary=True))

        expec = 0
        for key, value in freq.items():
            ones = key.count('1')
            if ones%2 == 0: 
                expec += float(value)/shots
            else:
                expec -= float(value)/shots
        print("Expectation value for the circuit: ", expec)
        return expec


def quantum_computer(
    observables, circuit1, circuit2, connection, shots=30000, verbose=False
):
    """Sends the execution to the quantum computer to calculate the expected
    value, instead of performing a simulation. (in process).

     :param observables: string.
     :param ciruit1: Circuit.
     :param ciruit2: Circuit.
     :param shots: int.
     :param verbose: bool.
     :return: reeconstruction value.
    """

    if verbose:
        print(f"qibo version: {qibo.__version__}")

    basis = ["X", "Y", "Z", "I"]
    states = ["0", "1", "+", "-", "+i", "-i"]

    sub_circuit_1_dimension = circuit1.nqubits
    if verbose:
        print("Sub circuit 1 dimension: ", sub_circuit_1_dimension)

    num_z_sub1 = observables[: (sub_circuit_1_dimension - 1)]
    num_z_sub2 = observables[(sub_circuit_1_dimension - 1) :]

    # first subcircuit:
    exp_value_1 = {}
    for b in basis:
        if verbose:
            print("Basis: ", b)
        # circuit_copy = circuit_1.copy(True)
        circuit1 = first_subcircuit_basis(circuit1, b)
        result = circuit1.execute_qc_compss(connection, nshots=shots)
        # obtain probability distribution
        freq = result.frequencies_compss(binary=True)
        # we call the function that computes the expectation value
        exp_value_1[b] = compute_expectation_value(
            freq, num_z_sub1 + b, shots=shots
        )
        if verbose:
            print(exp_value_1[b])

    # second subcircuit:
    exp_value_2 = {}

    for s in states:
        if verbose:
            print("States: ", s)
        circuit2 = second_subcircuit_states(circuit2, s)
        result = circuit2.execute_qc_compss(connection, nshots=shots)
        # obtain probability distribution
        freq = result.frequencies_compss(binary=True)
        # we call the function that computes the expectation value
        exp_value_2[s] = compute_expectation_value(
            freq, num_z_sub2 + b, shots=shots
        )
        if verbose:
            print(exp_value_2[s])

    exp_value_1 = compss_wait_on(exp_value_1)
    exp_value_2 = compss_wait_on(exp_value_2)

    reconstruction = (
        1
        / 2
        * (
            exp_value_1["I"] * (exp_value_2["0"] + exp_value_2["1"])
            + exp_value_1["X"] * (exp_value_2["+"] - exp_value_2["-"])
            + exp_value_1["Y"] * (exp_value_2["+i"] - exp_value_2["-i"])
            + exp_value_1["Z"] * (exp_value_2["0"] - exp_value_2["1"])
        )
    )
    print(
        "Expectation value after circuit cutting and reconstruction:",
        reconstruction,
    )
    return reconstruction


def analytical_solution(observables, circuit, verbose=False):
    """Calculate the analytical expected value of a whole circuit.

    :param observables: string.
    :param ciruit: Circuit.
    :param verbose: bool.
    :return: analytical expected value.
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

    print(
        "The expectation value of",
        expectation_value,
        "in the entire circuit is ",
        exp_full_circuit,
    )
    return exp_full_circuit
