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

"""
Qdislib circuit utils.

This file contains all auxiliary circuit classes and functions.
"""

import igraph
import inspect
import numpy as np
import qibo
import random
import typing
from collections import defaultdict
from qibo import models  # , callbacks
from qibo import gates
from qibo import hamiltonians
from Qdislib.utils.exceptions import QdislibException


def analytical_solution(
    observables: str, circuit: typing.Any, verbose: bool = False
) -> float:
    """Calculate the analytical expected value of a whole circuit.

    Example
    -------
    >>> analytical_value = analytical_solution(observables="ZZZZ", circuit=circuit, verbose=True)

    :param observables: String containing observables.
    :param circuit: Circuit object.
    :param verbose: Whether to print verbose output. Defaults to False.
    :return: Analytical expected value.
    """
    # TODO: circuit parameter is instantiable?

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
        new_expectation_value.expectation(state.state(numpy=True), normalize=False)
    )

    if verbose:
        print(
            "The expectation value of",
            expectation_value,
            "in the entire circuit is ",
            exp_full_circuit,
        )
    return exp_full_circuit


def random_circuit(
    qubits: int, gate_max, num_cz: typing.Optional[int], p: typing.Any
) -> models.Circuit:
    """Generate a random circuit.

    :param qubits: Number of qubits.
    :param gate_max: Maximum number of single-qubit gates between CZ.
    :param num_cz: Number of CZ qbit gates. If none, p will be taken.
    :param p: Probability of adding an edge.
    :raise QdislibException: If num_cz or p are both None.
    :return: New random circuit.
    """
    if p == None:
        graph = igraph.Graph.Erdos_Renyi(
            n=qubits, m=num_cz, directed=False, loops=False
        )
    elif num_cz == None:
        graph = igraph.Graph.Erdos_Renyi(n=qubits, p=p, directed=False, loops=False)
    else:
        raise QdislibException(
            "Error: only the number of edges or the probability of adding an edge must be specified"
        )

    # # Display the graph:
    # adj_mat=graph.get_adjacency()
    # fig, ax = plt.subplots()
    # igraph.plot(graph, target=ax, vertex_label=range(qubits))
    # graph.degree()

    edge_list = graph.get_edgelist()

    gates_pull = [gates.X, gates.H, gates.S, gates.T]  # pull of single-qubit gates
    circuit = models.Circuit(qubits)
    for edge in edge_list:
        # Number of single-qubit gates between CZ
        rand_tmp = random.randint(0, gate_max)
        for _ in range(rand_tmp):
            # Gate selected from the pull
            sel_gate = random.choice(gates_pull)
            # Qubit selected to apply the gate
            sel_qubit = random.randint(0, qubits - 1)
            circuit.add(sel_gate(sel_qubit))
        # 2-qubit gate from graph
        circuit.add(gates.CZ(edge[0], edge[1]))

    return circuit


def draw_to_circuit(
    text_draw: str, parameters: typing.Optional[typing.List[typing.Any]] = None
) -> models.Circuit:
    """Convert text circuit to circuit object.

    :param text_draw: Input text to convert.
    :param parameters: List of parameters, defaults to None
    :return: Circuit object.
    """
    split = text_draw.splitlines()
    # print(split)
    print(split)
    qubits_lst = []
    split = [element for element in split if element.strip()]
    split = [element for element in split if element != ""]
    print(split)
    for line in split:
        index = line.index("─")
        qubits_lst.append(line[index:])

    list_multiple_gates = defaultdict(list)
    # Now we will process each line to identify multi-qubit gates
    for idx, qubit_line in enumerate(qubits_lst):
        qubit_number = idx  # Line number corresponds to the qubit (q0 is index 0)
        qubit_state = list(qubit_line)

        # Boolean to track if we are inside a multi-qubit gate
        for i, symbol in enumerate(qubit_state):
            if symbol == "o":
                index = i
                for idx2, qubit in enumerate(qubits_lst[idx + 1 :]):
                    if list(qubit)[index] != "|":
                        name = list(qubit)[index]
                        if name == "Z":
                            name = "CZ"
                        elif name == "X":
                            name = "CNOT"
                        list_multiple_gates[idx].append((name, (idx, idx2 + idx + 1)))
                        qubits_lst[idx2 + idx + 1] = (
                            qubits_lst[idx2 + idx + 1][:index]
                            + "─"
                            + qubits_lst[idx2 + idx + 1][index + 1 :]
                        )
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
            # for idx2, char in enumerate(qubit_state):
            if char != "─" and char != "|":
                if char != "o":
                    if qubit_state[step + 1] == "─" and qubit_state[step - 1] == "─":
                        tmp = char
                        # print("Add gate: ", tmp, " qubit ", (idx,))
                        # circuit.add(getattr(gates, tmp)(idx))
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

                    elif qubit_state[step - 1] == "─" and qubit_state[step + 1] != "─":
                        tmp = ""
                        print(qubit_state)
                        print(qubit_state[step + 1])
                        print(range(step, num_steps))
                        for i in range(step, num_steps):
                            print(qubit_state[i])
                            if qubit_state[i + 1] == "─":
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
                                param_count = (
                                    len(signature.parameters) - 1
                                )  # exclude 'self'

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

                elif char == "o":
                    saved_qubit.append(idx)

        for idx in saved_qubit:
            # if list_multiple_gates[idx]:
            print(
                "Add gate: ",
                list_multiple_gates[idx][0][0],
                " qubit ",
                list_multiple_gates[idx][0][1],
            )
            circuit.add(
                getattr(gates, list_multiple_gates[idx][0][0])(
                    *list_multiple_gates[idx][0][1]
                )
            )
            list_multiple_gates[idx].remove(list_multiple_gates[idx][0])

    print(circuit.draw())
    return circuit
