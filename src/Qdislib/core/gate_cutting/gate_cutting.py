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

from collections import Counter
from functools import reduce
from itertools import product

import networkx as nx
import matplotlib.pyplot as plt

from Qdislib.api.api import *

def has_number_or_less(lst, number):
    for num in lst:
        if num <= number:
            return True
    return False

def gates_dict(circuit):
    double_gates = {}
    for index, gate in enumerate(circuit.queue):
        if len(gate.qubits) > 1:
            start, end = gate.qubits
            new_tuple = [(i, i + 1) for i in range(start, end)]
            # print(new_tuple)
            for tupl in new_tuple:
                if tupl not in double_gates:
                    double_gates[tupl] = []
                    double_gates[tupl].append(index + 1)
                else:
                    double_gates[tupl].append(index + 1)

    # print(double_gates)
    return double_gates

def generate_circuits(combinations_list,circuit,gates_cut,draw):
    generated_circuits = []
    for index2, combination in enumerate(combinations_list):
        circuit1 = circuit.copy(True)
        target_gates = []
        for gate in gates_cut:
            target_gates.append(circuit1.queue[gate - 1])

        if any(
            type(element) != type(target_gates[0]) for element in target_gates
        ):
            print("All the gates to cut have to be the same type")
            return None, None

        target_qubit = target_gates[0].control_qubits[0]

        for index, gate in enumerate(circuit1.queue):
            if gate in target_gates:
                if type(gate) == gates.CNOT:
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
                elif type(gate) == gates.CZ:
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
    # ------------------------------------
    # SPLIT IN 4 SUBCIRCUITS
    # ------------------------------------
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    combinations_list = generate_combinations(len(gates_cut), type_gates)

    observable_dict = {}
    for num_qubit in range(0, circuit.nqubits):
        observable_dict[num_qubit] = observables[num_qubit]
    if verbose:
        print(observable_dict)

    generated_circuits = generate_circuits(combinations_list,circuit,gates_cut,draw)

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

    list_subcircuits = concatenate_lists(list_subcircuits)
    list_observables = concatenate_lists(list_observables)

    if draw:
        for index, sub in enumerate(list_subcircuits):
            print("\n Subcircuit " + str(index + 1))
            print(sub.draw())
    return list_subcircuits, list_observables


def concatenate_lists(lst):
    conc_list = []
    for el in lst:
        conc_list = conc_list + el
    return conc_list


@task(returns=qibo.states.CircuitResult)
def gate_simulation(i, shots=30000):
    # ------------------------------------
    # SIMULATION
    # ------------------------------------
    result = i(nshots=shots)
    return result


@task(returns=list)
def gate_frequencies(result):
    # ------------------------------------
    # FREQUENCIES
    # ------------------------------------
    freq = dict(result.frequencies(binary=True))
    return freq


@task(returns=int)
def gate_expectation_value(freq, basis, shots):
    """This function computes the expectation value given a probability
    distribution (the output of the quantum computer) in a given basis that
    we choose.

     INPUT:
      - freq (dict): frequency distribution coming from the quantum computer.
      - basis (str): we aim to compute the expectation value of this set of
                     operators.  For example, "XYY" indicates that we
                     calculate the expectation value of X over the first qubit,
                     and the expectation value of Y over the second and
                     third qubits.
      - shots (int): Numer of times that we have runed the quantum computer,
                     needed to compute the probability in the probability
                     distribution.

      OUTPUT:
       - expectation_value (float): Final expectation value.

       This function assumes that for computing the 'X' and the 'Y' expectation
       value, the qubit state it is in the appropiate diagonal basis.

       For the moment, we only implement two types of cases for the basis:
       1) Only combinations of X, Y or/and Z. For example: 'XXYYXZ'
       2) Only a single I operator in the last position. For example 'ZYXXYI'.

    """

    expectation_value = 0
    for key, value in freq.items():
        if len(basis) != len(key):
            print("Not enough basis")
            return None
        result = "".join(char for char, bit in zip(basis, key) if bit == "1")
        not_i = len(result) - result.count("I")
        if not_i % 2 == 0:
            expectation_value += float(value) / shots
        else:
            expectation_value -= float(value) / shots

    return expectation_value


def gate_reconstruction(type_gates, gates_cut, exp_values, verbose=False):
    # --------------------------------------
    # RECONSTRUCTION
    # --------------------------------------
    num_generated = int(len(exp_values) / 2 ** len(gates_cut))
    if verbose:
        print(num_generated)
    result = [
        eval("*".join(map(str, exp_values[i: i + num_generated])))
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
def sum_dicts(dicts):
    summed_dict = reduce(lambda a, b: a + Counter(b), dicts, Counter())
    return dict(summed_dict)


def generate_combinations(n, gate_type):
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
    type_gates = type(circuit.queue[gates_cut[0] - 1])
    subcircuits, list_observables = split_gates(
        observables, gates_cut, circuit, draw
    )
    exp_value = []
    for index, i in enumerate(subcircuits):
        list_freq = []
        i.add(gates.M(*range(i.nqubits)))
        for p in range(0, chunk):
            result = gate_simulation(i, int(shots / chunk))
            freq = gate_frequencies(result)
            # frq in a list COLLECTION
            list_freq.append(freq)
        # task per sumar dicts COLLECTIONS
        total_freq = sum_dicts(list_freq)
        if verbose:
            print(total_freq)
        if verbose:
            print(list_observables)
        obs = list_observables[index]
        if verbose:
            print(obs)
        new_obs = "".join([value for key, value in sorted(obs.items())])
        exp_value.append(gate_expectation_value(total_freq, new_obs, shots))

    exp_value = compss_wait_on(exp_value)
    result = gate_reconstruction(type_gates, gates_cut, exp_value, verbose)
    return result
