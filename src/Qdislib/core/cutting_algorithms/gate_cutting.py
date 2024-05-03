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
from pycompss.api.parameter import COLLECTION_IN

import numpy as np
import qibo
from qibo import gates
from qibo.result import CircuitResult

from collections import Counter
from functools import reduce
from itertools import product

from Qdislib.core.cutting_algorithms._pycompss_functions import (
    _compute_expectation_value,
)
from Qdislib.utils.graph import gen_graph_circuit, _separate_observables


def _has_number_or_less(lst, number):
    """
    Check if a list contains a specific number or a smaller one.

    :param lst: list.
    :param number: int.
    :return: bool.
    """
    for num in lst:
        if num <= number:
            return True
    return False


def _has_number(lst, number):
    """
    Check if a list contains a specific number.

    :param lst: list.
    :param number: int.
    :return: bool.
    """
    for num in lst:
        if num == number:
            return True
    return False


def _gates_dict(circuit):
    """
    Convert the circuit queue gates to a dictionary
    depending on what qubits apply.

    :param circuit: Circuit.
    :return: _double_gates.
    """
    _double_gates = {}
    for index, gate in enumerate(circuit.queue):
        if len(gate.qubits) > 1:
            start, end = gate.qubits
            new_tuple = [(i, i + 1) for i in range(start, end)]
            for tupl in new_tuple:
                if tupl not in _double_gates:
                    _double_gates[tupl] = []
                    _double_gates[tupl].append(index + 1)
                else:
                    _double_gates[tupl].append(index + 1)
    return _double_gates


def _generate_circuits(combinations_list, circuit, gates_cut, draw):
    """
    Substitute the gates being cut for the equivalencies (only CZ right now) and
    generates all the combinations. After it checks subcircuits created being derivated
    from the cut.The return is a list with all the independent subcircuits.

    :param combinations_list: list.
    :param circuit: Circuit.
    :param gates_cut: int list.
    :param draw: bool
    :return: generated_circuits.
    """
    generated_circuits = []
    for index2, combination in enumerate(combinations_list):
        circuit1 = circuit.copy(True)
        target_gates = []
        for gate in gates_cut:
            target_gates.append(circuit1.queue[gate - 1])

        if any(
            not isinstance(element, type(target_gates[0]))
            for element in target_gates
        ):
            print("All the gates to cut have to be the same type")
            return None, None

        target_qubit = target_gates[0].control_qubits[0]

        for index, gate in enumerate(circuit1.queue):
            if gate in target_gates:
                if isinstance(gate, gates.CNOT):
                    idx = target_gates.index(gate)
                    control_qubit = gate.control_qubits[0]
                    target_qubit = gate.target_qubits[0]
                    if len(combination[idx]) > 2:
                        circuit1.queue[index] = combination[idx][0](
                            control_qubit, np.pi / 2
                        )
                        circuit1.queue.insert(
                            index + 1, combination[idx][1](control_qubit)
                        )
                        circuit1.queue.insert(
                            index + 2,
                            combination[idx][2](target_qubit, np.pi / 2),
                        )
                        circuit1.queue.insert(
                            index + 3, combination[idx][3](target_qubit)
                        )
                    else:
                        circuit1.queue[index] = combination[idx][0](
                            control_qubit, np.pi / 2
                        )
                        circuit1.queue.insert(
                            index + 1,
                            combination[idx][1](target_qubit, np.pi / 2),
                        )
                elif isinstance(gate, gates.CZ):
                    # print("CZ")
                    idx = target_gates.index(gate)
                    control_qubit = gate.control_qubits[0]
                    target_qubit = gate.target_qubits[0]
                    circuit1.queue[index] = combination[idx][0](control_qubit)
                    circuit1.queue.insert(
                        index + 1, combination[idx][1](target_qubit)
                    )
                    # changed_gates.append(
                    #     (circuit1.queue[index],circuit1.queue[index+1])
                    # )

        if draw:
            print("\n Circuit " + str(index2 + 1))
            print(circuit1.draw())
        generated_circuits.append(circuit1)
    return generated_circuits


def split_gates(observables, gates_cut, circuit, draw=False, verbose=False):
    """
    Description
    -----------
    Split a circuit into four subcircuits based on the provided gates' cut points.

    Parameters
    ----------
    observables: list
        List of observables.
    gates_cut: list
        List of integers indicating the cut points.
    circuit: object
        Circuit object.
    draw: bool, optional
        Whether to draw the subcircuits. Defaults to False.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    list_subcircuits: list
        List of subcircuits generated after splitting.
    list_observables: list
        List of observables associated with the subcircuits.

    Example
    -------
    >>> list_subcircuits, list_observables = split_gates(observables="ZZZZZ", gates_cut=[2, 5, 8],
    >>>                                                  circuit=circuit, draw=True, verbose=True)
    """
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    combinations_list = _generate_combinations(len(gates_cut), type_gates)

    observable_dict = _separate_observables(circuit, observables, verbose)

    generated_circuits = _generate_circuits(
        combinations_list, circuit, gates_cut, draw
    )

    list_subcircuits = []
    list_observables = []
    list_unpack = []
    for new_circuit in generated_circuits:
        new_list = gen_graph_circuit(new_circuit, observable_dict)
        list_unpack.append(new_list)

    list_unpack = compss_wait_on(list_unpack)

    for x in list_unpack:
        if len(x) > 1:
            list_subcircuits.append(x[0])
            list_observables.append(x[1])
        else:
            list_subcircuits.append(x[0])

    list_subcircuits = _concatenate_lists(list_subcircuits)
    list_observables = _concatenate_lists(list_observables)

    if draw:
        for index, sub in enumerate(list_subcircuits):
            print("\n Subcircuit " + str(index + 1))
            print(sub.draw())
    return list_subcircuits, list_observables


def _concatenate_lists(lst):
    """
    Concatenate all list inside the list converting
    it to a 1D array.

    :param lst: 2D list.
    :return: concatenated_list.
    """
    conc_list = []
    for el in lst:
        conc_list = conc_list + el
    return conc_list


@task(returns=CircuitResult)
def _gate_simulation(circuit, shots=30000):
    """
    Execute a circuit.

    :param circuit: Circuit.
    :param shots: int.
    :return: result.
    """
    result = circuit(nshots=shots)
    return result


@task(returns=list)
def _gate_frequencies(result):
    """
    Calculate frequencies from a result.

    :param result: CircuitResult.
    :return: frequencies.
    """
    freq = dict(result.frequencies(binary=True))
    return freq


def gate_reconstruction(type_gates, gates_cut, exp_values, verbose=False):
    """
    Description
    -----------
    This function calculates the reconstruction of a circuit after cutting a set of gates and executing the subcircuits. It computes the expected value based on the provided gate type, the list of cut gates, and the list of expected values.

    Parameters
    ----------
    type_gates: gates
        The type of gates that were cut from the circuit (e.g., "CZ", "CNOT").
    gates_cut: list of int
        A list containing the indices of the gates that were cut from the circuit.
    exp_values: list of float
        A list of expected values obtained from executing the subcircuits.
    verbose: bool, optional
        If True, additional information will be printed during the computation. Defaults to False.

    Returns
    -------
    reconstruction: float
        The reconstructed expected value of the circuit after gate cutting.

    Example
    -------
    >>> type_gates = gates.CZ
    >>> gates_cut = [2, 4]
    >>> exp_values = [0.5, 0.3, 0.6, 0.4, 0.7, 0.1, 0.8, 0.2]
    >>> reconstruction = gate_reconstruction(type_gates, gates_cut, exp_values, verbose=True)
    """
    num_generated = int(len(exp_values) / 2 ** len(gates_cut))
    if verbose:
        print(num_generated)
    result = [
        eval("*".join(map(str, exp_values[i : i + num_generated])))
        for i in range(0, len(exp_values), num_generated)
    ]
    if verbose:
        print(exp_values)
    if verbose:
        print(result)
    result1 = [x * 1j if i % 2 == 0 else x for i, x in enumerate(result)]
    result2 = [x * 1j if i % 2 != 0 else x for i, x in enumerate(result)]
    if type_gates == gates.CZ:
        reconstruction = (
            (1 / (1j + 1))
            * (sum(result1) + sum(result2))
            / 2 ** len(gates_cut)
        )
    elif type_gates == gates.CNOT:
        reconstruction = (
            1j
            * np.exp(-1j * np.pi / 4)
            / np.sqrt(2)
            * (sum(result1) + sum(result2))
            / 2 ** len(gates_cut)
        )
    if verbose:
        print("\n")
    if verbose:
        print("Reconstructed expected value: ", reconstruction)
    if verbose:
        print("Absolute value of reconstruction ", np.absolute(reconstruction))
    if verbose:
        print("Reconstruction value: ", reconstruction.real)
    return reconstruction.real


@task(returns=dict, dicts=COLLECTION_IN)
def _sum_dicts(dicts):
    """
    Sum all dictionaries by key resulting in only one
    dictionary.

    :param dics: dict list.
    :return: summed_dict.
    """
    summed_dict = reduce(lambda a, b: a + Counter(b), dicts, Counter())
    return dict(summed_dict)


def _generate_combinations(n, gate_type):
    """
    Generate combinations for gate cutting
    depending on the type of the gate and the
    number of gates being cut.

    :param n: int.
    :param gate_type: string
    :return: all_combinations.
    """
    objects = []
    if gate_type == gates.CZ:
        objects = [(gates.S, gates.S), (gates.SDG, gates.SDG)]
    elif gate_type == gates.CNOT:
        objects = [
            (gates.RZ, gates.RX),
            (gates.RZ, gates.Z, gates.RX, gates.X),
        ]
    all_combinations = list(product(objects, repeat=n))
    return all_combinations


def gate_cutting(
    observables,
    circuit,
    gates_cut,
    shots=30000,
    chunk=1,
    draw=False,
    verbose=False,
):
    """
    Description
    -----------
    Cut a circuit by removing a set of gates, executes each resulting subcircuit, and reconstructs the expected value.

    Parameters
    ----------
    observables: string
        Observables.
    circuit: Circuit
        Circuit object.
    gates_cut: list of int
        List of integers indicating the cut points.
    shots: int, optional
        Number of shots for simulation. Defaults to 30000.
    chunk: int, optional
        Chunk size for simulation. Defaults to 1.
    draw: bool, optional
        Whether to draw the subcircuits. Defaults to False.
    verbose: bool, optional
        Whether to print verbose output. Defaults to False.

    Returns
    -------
    reconstruction: float
        Reconstruction of the expected value.

    Example
    -------
    >>> reconstruction = gate_cutting(observables="ZZZZZ", circuit=circuit, gates_cut=[2, 5, 8],
    >>>                               shots=30000, chunk=1, draw=True, verbose=True)
    """
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    subcircuits, list_observables = split_gates(
        observables, gates_cut, circuit, draw
    )
    exp_value = []
    for index, i in enumerate(subcircuits):
        list_freq = []
        i.add(gates.M(*range(i.nqubits)))
        for _ in range(0, chunk):
            result = _gate_simulation(i, int(shots / chunk))
            freq = _gate_frequencies(result)
            # frq in a list COLLECTION
            list_freq.append(freq)
        # task per sumar dicts COLLECTIONS
        total_freq = _sum_dicts(list_freq)
        if verbose:
            print(total_freq)
        if verbose:
            print(list_observables)
        obs = list_observables[index]
        if verbose:
            print(obs)
        new_obs = "".join([value for key, value in sorted(obs.items())])
        exp_value.append(
            _compute_expectation_value(total_freq, new_obs, shots)
        )

    exp_value = compss_wait_on(exp_value)
    reconstruction = gate_reconstruction(
        type_gates, gates_cut, exp_value, verbose
    )
    return reconstruction
