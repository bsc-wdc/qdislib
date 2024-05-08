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
import networkx as nx

from Qdislib.core.cutting_algorithms._pycompss_functions import (
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

from pycompss.api.api import compss_wait_on
import time


def wire_cutting(
    observables, circuit, gate_tuple, shots=30000, draw=False, verbose=False
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
        lst_observables, qubit, list_subcircuits[0], list_subcircuits[1], shots, verbose
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

    non_empty_list = []

    list_subcircuits = _partition_circuit(
        subgraphs, dag, circuit, non_empty_list, verbose=False
    )

    if verbose:
        print(qubits_first_gate, qubits_second_gate)
    new_qubit = [first_gate.qubits[pos1], second_gate.qubits[pos2]]
    if verbose:
        print(new_qubit)
    if verbose:
        print(first_gate.qubits, second_gate.qubits)

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
    lst_observables = [obs1, obs2]

    if first_gate not in list_subcircuits[0].queue:
        list_subcircuits.reverse()
        lst_observables.reverse()
    if draw:
        for circ in list_subcircuits:
            print(circ.draw())
            print("\n")
    return new_qubit, list_subcircuits, lst_observables


def simulation(
    list_observables,
    qubit,
    circuit_1,
    circuit_2=None,
    shots=30000,
    verbose=False,
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
            copy_circuit1 = models.Circuit(circuit_1.nqubits)
            copy_circuit = circuit_1.copy(True)
            copy_circuit1.queue = copy_circuit.queue
            circuit1 = _first_subcircuit_basis(copy_circuit1, b, qubit[0])
            result = circuit1.execute_compss(nshots=shots)
            # obtain probability distribution
            freq = result.frequencies_compss(binary=True)
            # we call the function that computes the e
            # new_obs = num_z_sub1[:qubit[0]] + b + num_z_sub1[qubit[0]:]
            # print(new_obs)
            for key, value in observables_1.items():
                if value == "-":
                    observables_1[key] = b
            new_obs = "".join(
                [value for key, value in sorted(observables_1.items())]
            )
            if verbose:
                print("OBSERVABLES: ", new_obs)
            exp_value_1[b] = _compute_expectation_value(
                freq, new_obs, shots=shots
            )

        # second subcircuit:
        exp_value_2 = {}
        if verbose:
            print("SECOND ONE")
        for s in states:
            if verbose:
                print("States: ", s)
            copy_circuit2 = models.Circuit(circuit_2.nqubits)
            copy_circuit = circuit_2.copy(True)
            copy_circuit2.queue = copy_circuit.queue
            circuit2 = _second_subcircuit_states(copy_circuit2, s, qubit[1])
            result = circuit2.execute_compss(nshots=shots)
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
    return reconstruction


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
