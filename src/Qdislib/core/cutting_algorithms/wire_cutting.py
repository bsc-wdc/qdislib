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

import numpy as np
import qibo
from qibo import models, gates, hamiltonians  # , callbacks
import networkx as nx

import numpy as np
import networkx as nx
from qibo import models, gates
import numpy as np
import networkx as nx
from Qdislib.api import *
from qiboconnection.connection import ConnectionConfiguration
from qiboconnection.api import API
from Qdislib.api import *
from scipy.optimize import minimize
import time
import itertools

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.parameter import *

import igraph as ig
import matplotlib.pyplot as plt
import random


import networkx as nx
from qibo import gates, models

import math


'''from Qdislib.core.cutting_algorithms._pycompss_functions import (
    _first_subcircuit_basis,
    _second_subcircuit_states,
    _compute_expectation_value,
)
from Qdislib.utils.graph import (
    _DAGgraph,
    _build_dag,
    _separate_observables,
    _create_graph,
    _partition_circuit,
)
from Qdislib.core.qubit_mapping.qubit_mapping import *
import time


def wire_cutting(
    observables,
    circuit,
    gate_tuple,
    shots=30000,
    draw=False,
    verbose=False,
    sync=True,
):
    """Description
    -----------
    Implement the algorithm of circuit cutting in one function. Cuts the circuit in 2 between two specified gates, and calculates the expected value of the reconstruction.

    Parameters
    ----------
    observables: str
        Observable.
    circuit: Circuit
        Circuit object.
    gate_tuple: tuple of int
        Tuple containing two integers specifying the gates between which the circuit will be cut.
    shots: int, optional
        Number of shots for simulation. Defaults to 30000.
    draw: bool, optional
        Whether to draw the subcircuits. Defaults to False.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    reconstruction: float
        Reconstruction value.

    Example
    -------
    >>> reconstruction = wire_cutting(observables="ZZZZZ", circuit=circuit, gate_tuple=(2, 5),
    >>>                               shots=30000, draw=True, verbose=True)

    """

    # Show qibo version
    if verbose:
        print(f"qibo version: {qibo.__version__}")

    gate_tuple = gate_tuple[0]
    qubit, list_subcircuits, lst_observables = split(
        observables, circuit, gate_tuple, draw, verbose
    )
    reconstruction = simulation(
        lst_observables,
        qubit,
        list_subcircuits[0],
        list_subcircuits[1],
        shots,
        verbose,
        sync=sync,
    )
    return reconstruction


def split(observables, circuit, gate_tuple, draw=False, verbose=False):
    """Description
    -----------
    Split a circuit into two subcircuits by cutting a qubit between the specified gates.

    Parameters
    ----------
    observables: str
        String containing observables.
    circuit: Circuit
        Circuit object.
    gate_tuple: tuple of int
        Tuple containing two integers specifying the gates between which the circuit will be cut.
    draw: bool, optional
        Whether to draw the subcircuits. Defaults to False.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    qubit: list
        List containing the qubits between which the circuit was cut.
    list_subcircuits: list of Circuit
        List of subcircuits generated after splitting.
    list_observables: list
        List of dictionaries containing observables associated with each subcircuit.

    Example
    -------
    >>> qubit, list_subcricuits, list_observables = split(observables="ZZZZZ", circuit=circuit,
    >>>                                                   gate_tuple)=(2, 5), draw=True, verbose=True)

    """

    circuit = circuit.copy()
    qubits_first_gate = circuit.queue[gate_tuple[0] - 1].qubits
    qubits_second_gate = circuit.queue[gate_tuple[1] - 1].qubits
    first_gate = circuit.queue[gate_tuple[0] - 1]
    second_gate = circuit.queue[gate_tuple[1] - 1]
    if verbose:
        print(first_gate, second_gate)
    qubit = list(
        set(circuit.queue[gate_tuple[0] - 1].qubits).intersection(
            set(circuit.queue[gate_tuple[1] - 1].qubits)
        )
    )
    qubit = qubit[0]
    pos1 = qubits_first_gate.index(qubit)
    pos2 = qubits_second_gate.index(qubit)

    observable_dict = _separate_observables(circuit, observables, verbose)

    digraph = nx.Graph()
    dag = _DAGgraph()

    _build_dag(circuit, dag)
    _create_graph(dag, digraph)

    red_edges = [
        (edge[0], edge[1])
        for edge in digraph.edges.data("color")
        if edge[2] == "red"
    ]
    digraph.remove_edges_from(red_edges)

    digraph.remove_edge(gate_tuple[0], gate_tuple[1])

    subgraphs = list(nx.connected_components(digraph))
    if len(subgraphs) > 2:
        print("MORE THAN 2 SUBGRAPH")
        return None, None

    list_subcircuits, non_empty_list = _partition_circuit(
        subgraphs, dag, circuit, verbose=False
    )

    if verbose:
        print(qubits_first_gate, qubits_second_gate)
    new_qubit = [first_gate.qubits[pos1], second_gate.qubits[pos2]]
    if verbose:
        print(new_qubit)
    if verbose:
        print(first_gate.qubits, second_gate.qubits)

    obs1, obs2 = create_wire_observables(
        non_empty_list,
        list_subcircuits,
        first_gate,
        new_qubit,
        observable_dict,
        verbose,
    )

    lst_observables = [obs1, obs2]

    check_reverse(first_gate, list_subcircuits, lst_observables)

    """if draw:
        for circ in list_subcircuits:
            print(circ.draw())
            print("\n")"""
    return new_qubit, list_subcircuits, lst_observables


@task(list_subcircuits=COLLECTION_IN, lst_observables=COLLECTION_IN)
def check_reverse(first_gate, list_subcircuits, lst_observables):
    if first_gate not in list_subcircuits[0].queue:
        list_subcircuits.reverse()
        lst_observables.reverse()


@task(returns=2, non_empty_list=COLLECTION_IN, list_subcircuits=COLLECTION_IN)
def create_wire_observables(
    non_empty_list,
    list_subcircuits,
    first_gate,
    new_qubit,
    observable_dict,
    verbose=False,
):
    obs1 = {}
    obs2 = {}

    for index, x in enumerate(non_empty_list[0]):
        if first_gate in list_subcircuits[0].queue:
            if index == new_qubit[0]:
                obs1[index] = "-"
            else:
                obs1[index] = observable_dict[x]
        else:
            obs1[index] = observable_dict[x]

    for index, x in enumerate(non_empty_list[1]):
        if first_gate not in list_subcircuits[0].queue:
            if index == new_qubit[0]:
                obs2[index] = "-"
            else:
                obs2[index] = observable_dict[x]
        else:
            obs2[index] = observable_dict[x]

    if verbose:
        print(obs1)
    if verbose:
        print(obs2)
    return obs1, obs2


def simulation(
    list_observables,
    qubit,
    circuit_1,
    circuit_2=None,
    shots=30000,
    verbose=False,
    sync=True,
):
    """Description
    -----------
    Perform the execution of a circuit to calculate the expected value. It accepts one or two circuits. With one circuit, it calculates the expected value straightforwardly. With two circuits, it performs a reconstruction in order to provide the expected value.

    Parameters
    ----------
    list_observables: list of dict
        List of dictionaries containing observables.
    qubits: list
        List of qubits where the cut was performed.
    circuit_1: Circuit
        First circuit object.
    circuit_2: Circuit, optional
        Second circuit object. Defaults to None.
    shots: int, optional
        Number of shots for simulation. Defaults to 30000.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    reconstruction value: float.
        Reconstruction value.

    Example
    -------
    >>> reconstrution = simulation(list_obervables=[observables1, observables2], qubits=[0,2],
    >>>                            circuit1, circuit2, shots=30000, verbose=True)

    """

    if verbose:
        print(f"qibo version: {qibo.__version__}")

    if circuit_2 is not None:
        # SIMULATION

        basis = ["X", "Y", "Z", "I"]
        states = ["0", "1", "+", "+i"]

        observables_1 = list_observables[0]
        observables_2 = list_observables[1]
        # num_z_sub1 = observables[:(sub_circuit_1_dimension - 1)]
        # num_z_sub2 = observables[(sub_circuit_1_dimension - 1):]

        # first subcircuit:
        exp_value_1 = {}
        for b in basis:
            if verbose:
                print("Basis: ", b)
            copy_circuit = models.Circuit(circuit_1.nqubits)
            temp = circuit_1.copy()
            copy_circuit.queue = temp.queue
            copy_circuit = _first_subcircuit_basis(copy_circuit, b, qubit[0])
            result = copy_circuit.execute_compss(nshots=shots)
            # obtain probability distribution
            freq = result.frequencies_compss(binary=True)
            # we call the function that computes the e
            # new_obs = num_z_sub1[:qubit[0]] + b + num_z_sub1[qubit[0]:]
            # print(new_obs)

            # new_obs = "".join(
            #    [value for key, value in sorted(observables_1.items())]
            # )
            if verbose:
                print("OBSERVABLES 1: ", observables_1)
                print(circuit_1.draw())
                print(circuit_1.nqubits)
                print(copy_circuit.circuit.draw())
                print(copy_circuit.circuit.nqubits)
            exp_value_1[b] = _compute_expectation_value(
                freq, observables_1, shots=shots, wire_observables=True, b=b
            )

        # second subcircuit:
        exp_value_2 = {}
        if verbose:
            print("SECOND ONE")
        for s in states:
            if verbose:
                print("States: ", s)
            copy_circuit2 = models.Circuit(circuit_2.nqubits)
            temp2 = circuit_2.copy(True)
            copy_circuit2.queue = temp2.queue
            copy_circuit2 = _second_subcircuit_states(copy_circuit2, s, qubit[1])
            result = copy_circuit2.execute_compss(nshots=shots)
            # obtain probability distribution
            freq = result.frequencies_compss(binary=True)
            # new_obs2 = "".join(
            #    [value for key, value in sorted(observables_2.items())]
            # )
            # we call the function that computes the expectation value
            if verbose:
                print("OBSERVABLES 2: ", observables_2)
            exp_value_2[s] = _compute_expectation_value(
                freq, observables_2, shots=shots
            )

        # exp_value_1 = compss_wait_on(exp_value_1)
        # exp_value_2 = compss_wait_on(exp_value_2)

        reconstruction = wire_reconstruction(exp_value_1, exp_value_2)
        if sync:
            reconstruction = compss_wait_on(reconstruction)
            if verbose:
                print(
                    "Expectation value after circuit cutting and reconstruction:",
                    reconstruction,
                )
        return reconstruction
    else:
        circuit_1.add(gates.M(*range(circuit_1.nqubits)))
        result = circuit_1(nshots=shots)
        freq = dict(result.frequencies(binary=True))

        expec = 0
        for key, value in freq.items():
            ones = key.count("1")
            if ones % 2 == 0:
                expec += float(value) / shots
            else:
                expec -= float(value) / shots
        print("Expectation value for the circuit: ", expec)
        return expec


@task(exp_value_1=DICTIONARY_IN, exp_value_2=DICTIONARY_IN, returns=1)
def wire_reconstruction(exp_value_1, exp_value_2):
    reconstruction = (
        1
        / 2
        * (
            (exp_value_1["I"] + exp_value_1["Z"]) * exp_value_2["0"]
            + (exp_value_1["I"] - exp_value_1["Z"]) * exp_value_2["1"]
            + exp_value_1["X"]
            * (2 * exp_value_2["+"] - exp_value_2["0"] - exp_value_2["1"])
            + exp_value_1["Y"]
            * (2 * exp_value_2["+i"] - exp_value_2["0"] - exp_value_2["1"])
        )
    )
    return reconstruction


def execute_qc(
    connection, qubit, circuit_1, circuit_2=None, shots=30000, verbose=False
):
    """Description
    -----------
    Perform the execution of a circuit by sending it to the Quantum Computer and obtaining the job IDs. Accepts one or two circuits. With one circuit, it executes it straightforwardly. With two circuits, it executes them separately and returns the job IDs for each.

    Parameters
    ----------
    connection: object
        API configuration.
    qubit: list
        List of qubits where the cut was performed.
    circuit_1: Circuit
        First circuit object.
    circuit_2: Circuit, optional
        Second circuit object. Defaults to None.
    shots: int, optional
        Number of shots for simulation. Defaults to 30000.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    job_ids1: list
        List of job IDs for the first circuit.
    job_ids2: list
        List of job IDs for the second circuit, if provided. Otherwise, None.

    Example
    -------
    >>> job_ids1, job_ids2 = execute_qc(connection, qubits=[0,2], circuit1, circuit2, shots=30000, verbose=True)
    """

    connection.select_device_ids(device_ids=[9])
    connection.list_devices()

    if circuit_2 is not None:
        # SIMULATION

        basis = ["X", "Y", "Z", "I"]
        states = ["0", "1", "+", "+i"]

        job_ids1 = []
        for b in basis:
            if verbose:
                print("Basis: ", b)
            copy_circuit1 = models.Circuit(circuit_1.nqubits)
            copy_circuit = circuit_1.copy(True)
            copy_circuit1.queue = copy_circuit.queue
            final_circuit = _first_subcircuit_basis(copy_circuit1, b, qubit[0])
            """G = architecture_x()
            G1 = _qubit_arch(final_circuit.circuit)
            dict_arch = subgraph_matcher(G,G1)
            print(dict_arch)
            mapped_circuit = rename_qubits(final_circuit,2,dict_arch[0],'B')
            print(mapped_circuit.circuit.draw())"""
            job_ids = final_circuit.execute_qc_compss(connection, nshots=shots)
            # job_ids = connection.execute(circuit=circuit_1, nshots=1000)
            job_ids1.append(job_ids[0])

        # second subcircuit:
        job_ids2 = []
        for s in states:
            if verbose:
                print("States: ", s)
            copy_circuit2 = models.Circuit(circuit_2.nqubits)
            copy_circuit = circuit_2.copy(True)
            copy_circuit2.queue = copy_circuit.queue
            final_circuit = _second_subcircuit_states(
                copy_circuit2, s, qubit[1]
            )
            job_ids = final_circuit.execute_qc_compss(connection, nshots=shots)
            # job_ids = connection.execute(circuit=final_circuit, nshots=1000)
            job_ids2.append(job_ids[0])

        job_ids1 = compss_wait_on(job_ids1)
        job_ids2 = compss_wait_on(job_ids2)
    else:
        print("Missing subcircuit2")
        print("NOT YET IMPLEMENTED!")
        job_ids1, job_ids2 = None, None
    return job_ids1, job_ids2


def reconstruction_qc(
    connection,
    job_ids1,
    job_ids2,
    list_observables,
    shots=30000,
    verbose=False,
):
    """Description
    -----------
    Perform the reconstruction of a circuit after sending it to the Quantum Computer and retrieving the job IDs. Accepts job IDs for one or two circuits and the list of observables associated with each subcircuit. Calculates the expectation value of the reconstructed circuit.

    Parameters
    ----------
    connection: object
        API configuration.
    job_ids1: list
        List of job IDs for the first circuit.
    job_ids2: list
        List of job IDs for the second circuit, if provided. Otherwise, None.
    list_observables: list of dict
        List of dictionaries containing observables for each subcircuit.
    shots: int, optional
        Number of shots for simulation. Defaults to 30000.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    reconstruction: float
        Reconstruction value.

    Example
    -------
    >>> reconstruction = reconstruction_qc(connection, job_ids1, job_ids2,
    >>>                                    list_observables=[observables1, observables2],
    >>>                                    shots=30000, verbose=True)
    """
    connection.select_device_ids(device_ids=[9])
    connection.list_devices()

    results1 = connection.get_results(job_ids=job_ids1)
    results2 = connection.get_results(job_ids=job_ids2)

    while any(x is None for x in results1) or any(x is None for x in results2):
        print("RETRIVEING RESULTS...")
        results1 = connection.get_results(job_ids=job_ids1)
        results2 = connection.get_results(job_ids=job_ids2)
        time.sleep(1)

    basis = ["X", "Y", "Z", "I"]
    states = ["0", "1", "+", "+i"]

    observables_1 = list_observables[0]
    observables_2 = list_observables[1]

    # first subcircuit:
    exp_value_1 = {}
    for index, b in enumerate(basis):
        if verbose:
            print("Basis: ", b)
        result = results1[index]
        # obtain probability distribution
        freq = result.frequencies_compss(binary=True)
        # we call the function that computes the expectation value
        for key, value in observables_1.items():
            if value == "-":
                observables_1[key] = b
        new_obs = "".join(
            [value for key, value in sorted(observables_1.items())]
        )
        if verbose:
            print("OBSERVABLES: ", new_obs)
        exp_value_1[b] = _compute_expectation_value(freq, new_obs, shots=shots)

    # second subcircuit:
    exp_value_2 = {}
    if verbose:
        print("SECOND ONE")
    for index, s in enumerate(states):
        if verbose:
            print("States: ", s)
        result = results2[index]
        # obtain probability distribution
        freq = result.frequencies_compss(binary=True)
        new_obs2 = "".join(
            [value for key, value in sorted(observables_2.items())]
        )
        # we call the function that computes the expectation value
        if verbose:
            print("OBSERVABLES: ", new_obs2)
        exp_value_2[s] = _compute_expectation_value(
            freq, new_obs2, shots=shots
        )

    exp_value_1 = compss_wait_on(exp_value_1)
    exp_value_2 = compss_wait_on(exp_value_2)
    reconstruction = (
        1
        / 2
        * (
            (exp_value_1["I"] + exp_value_1["Z"]) * exp_value_2["0"]
            + (exp_value_1["I"] - exp_value_1["Z"]) * exp_value_2["1"]
            + exp_value_1["X"]
            * (2 * exp_value_2["+"] - exp_value_2["0"] - exp_value_2["1"])
            + exp_value_1["Y"]
            * (2 * exp_value_2["+i"] - exp_value_2["0"] - exp_value_2["1"])
        )
    )

    print(
        "Expectation value after circuit cutting and reconstruction:",
        reconstruction,
    )
    return reconstruction


def wire_cutting_QC(
    connection,
    list_observables,
    qubit,
    circuit1,
    circuit2,
    shots=10,
    verbose=False,
):
    """Description
    -----------
    Send the circuits to the Quantum Computer in order to be executed and retrieves the information to calculate the reconstruction value.

    Parameters
    ----------
    connection: string
        Connection string.
    list_observables: list of dict
        List of dictionaries containing observables.
    qubit: list
        List of qubits where the cut was performed.
    circuit1: Circuit
        First circuit object.
    circuit2: Circuit
        Second circuit object.
    shots: int, optional
        Number of shots for simulation. Defaults to 10.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    reconstruction value: float
        Reconstruction value.

    Example
    -------
    >>> reconstruction = quantum_computer(connection, list_observables=[observables1, observables2], qubit=[0,2],
    >>>                                   circuit1, circuit2, shots=10, verbose=True)
    """

    if verbose:
        print(f"qibo version: {qibo.__version__}")

    job_ids1, job_ids2 = execute_qc(
        connection, qubit, circuit1, circuit2, shots
    )
    reconstruction = reconstruction_qc(
        connection, job_ids1, job_ids2, list_observables, shots
    )
    print(
        "Expectation value after circuit cutting and reconstruction:",
        reconstruction,
    )
    return reconstruction'''


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
    >>> analytical_value = analytical_solution(observables="ZXYI", circuit=circuit, verbose=True)
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



def circuit_to_dag(circuit, num_qubits):
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
        if gate_name == "Measurement":
            continue
        
        # Get the qubits this gate acts on
        
        qubits = node_data['qubits']
        parameters = node_data['parameters']

        # Reconstruct the gate based on the gate type and qubits
        if gate_name == 'H':
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

        else:
            raise ValueError(f"Unsupported gate type: {gate_name}")

    # Optionally handle measurements, assuming all qubits are measured at the end
    for node in topo_order:
        node_data = dag.nodes[node]
        if node_data['gate'] == "Measurement":
            circuit.add(gates.M(node_data['qubit']))

    return circuit


def wire_cutting(rand_qc,cut,sync=True):
    dag = circuit_to_dag(rand_qc, num_qubits=rand_qc.nqubits)

    '''if nx.number_connected_components(dag.to_undirected()) > 1:
        S = [dag.subgraph(c).copy() for c in nx.connected_components(dag)]
        for s in S:
            if '''

    graphs = generate_wire_cutting(dag,cut , num_qubits=rand_qc.nqubits)
    if sync:
        graphs = compss_wait_on(graphs)
    final_recons = 1/(2**len(cut))*sum(graphs)
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
                            temp_list[i] = num_qubits+index-1
            
                    updated_tuple = tuple(temp_list)
                    dag.nodes[successor]['qubits'] = updated_tuple


        dag.add_node(f"O_{index}", gate='S', qubits=common_qubit, parameters=())
        
        # Add the new end node with the same properties as the target node
        dag.add_node(f"PS_{index}", gate='T', qubits=(num_qubits+index-1,), parameters=())

        dag.add_edge(source, f"O_{index}", color="blue")
        dag.add_edge(f"PS_{index}", target, color="blue")
    
        copy_dag = dag.copy()
        red_edges = []
        for ed in dag.edges:
            if dag.get_edge_data(ed[0],ed[1])["color"] == "red":
                red_edges.append(ed)
        
        copy_dag.remove_edges_from(red_edges)

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

        graph = generate_subcircuits_wire_cutting(graph, num_qubits+len(edges_to_replace),index, edges_to_replace, graph_components )

        exp_value = []
        for s in graph_components:
            s_new, highest_qubit = update_qubits(s)
            subcirc = dag_to_circuit(s_new,highest_qubit)
            expected_value = execute_subcircuits(subcirc, index)
            exp_value.append(expected_value)

        exp_value = change_sign(exp_value, index)
        reconstruction.append(exp_value)

    return reconstruction

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

    #graphs = [updated_dag]

    number = idx
    list_substitutions = []
    
    digit = number % 8  # Get the last digit
    list_substitutions.append(digit)
    number //= 8 # Move to the next digit
    
    while number != 0:
        digit = number % 8  # Get the last digit
        list_substitutions.append(digit)
        number //= 8 # Move to the next digit

    list_substitutions = list(reversed(list_substitutions))
        
        #index = idx // (8**(idx2-1))

    for idx2, index in enumerate(list_substitutions, start=0):
        idx2 = len(list_substitutions) - idx2
        
        # I 0 
        if index == 0:
            updated_dag.remove_node(f'O_{idx2}')
            updated_dag.remove_node(f'PS_{idx2}')

        # I 1
        elif index == 1:
            updated_dag.remove_node(f'O_{idx2}')
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
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'H'
            updated_dag.add_node(f"O2_{idx2}", gate='SDG', qubits=updated_dag.nodes[f'O_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'O_{idx2}', f'O2_{idx2}', color="blue")
            updated_dag.nodes[f'PS_{idx2}']['gate'] = 'S'
            updated_dag.add_node(f"PS2_{idx2}", gate='H', qubits=updated_dag.nodes[f'PS_{idx2}'].get('qubits'), parameters=())
            updated_dag.add_edge(f'PS2_{idx2}', f'PS_{idx2}', color="blue")

        # Y -i
        elif index == 5:
            updated_dag.nodes[f'O_{idx2}']['gate'] = 'H'
            updated_dag.add_node(f"O2_{idx2}", gate='SDG', qubits=updated_dag.nodes[f'O_{idx2}'].get('qubits'), parameters=())
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
        graph_components[i].add_edges_from(new_subgraph.edges(data=True))
        
    return updated_dag

@task(returns=1)
def execute_subcircuits(subcirc, index):
    qibo.set_backend("numpy")
    subcirc.add(gates.M(*range(subcirc.nqubits))) #
    result = subcirc(nshots=30000)
    freq = dict(result.frequencies(binary=True))

    expectation_value = 0
    for key, value in freq.items():
        contribution = 1
        for bit, obs in zip(key, 'Z'*subcirc.nqubits):
            if obs == "Z":
                contribution *= (-1) ** int(bit)
            elif obs == "I":
                contribution *= 1
            else:
                raise ValueError(f"Unsupported observable {obs}")
        
        # Add the contribution weighted by its frequency
        expectation_value += contribution * (value / 30000)

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
            0, gate_max, seed=10
        )  # number of single-qubit gates between CZ
        for j in range(rand_tmp):
            sel_gate = random.choice(gates_pull, seed=10)  # gate selected from the pull
            sel_qubit = random.randint(
                0, qubits - 1, seed=10
            )  # qubit selected to apply the gate
            circuit.add(sel_gate(sel_qubit))

        circuit.add(
            gates.CZ(edge_list[i][0], edge_list[i][1])
        )  # 2-qubit gate from graph

    return circuit

