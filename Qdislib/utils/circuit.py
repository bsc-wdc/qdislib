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
    circuit: typing.Any, observables: str, verbose: bool = False
) -> float:
    """
    Compute the analytical expectation value of a given quantum circuit.

    This function evaluates the expectation value of a quantum circuit with respect
    to a set of observables, using a symbolic Hamiltonian approach. It assumes
    the circuit returns the full quantum statevector.

    :param observables:
        A string representing the sequence of observables to apply (e.g., "XZI").
        Each character corresponds to a Pauli operator ("X", "Y", "Z", "I") applied
        on successive qubits.
    :param circuit:
        A callable object (such as a Qibo circuit) that, when executed, returns the final quantum state.
    :param verbose:
        If True, prints intermediate symbolic and numerical expectation values.
    :return:
        The computed expectation value as a floating-point number.

    :example:

    Suppose we have a simple quantum circuit and wish to calculate the expectation value of the Pauli operators "XZI".

    Example usage:

    .. code-block:: python

        from qibo import models, gates
        from mymodule import analytical_solution

        # Create a quantum circuit
        circuit = models.Circuit(3)
        circuit.add(gates.H(0))         # Apply Hadamard on qubit 0
        circuit.add(gates.CNOT(0, 1))    # Apply CNOT between qubit 0 and 1

        # Define the observables
        observables = "XZI"  # Apply Pauli-X to qubit 0, Z to qubit 1, Identity to qubit 2

        # Compute the analytical expectation value
        result = analytical_solution(observables, circuit, verbose=True)

        print(f"Expectation value: {result}")

    In this example:
    - The `observables` string `"XZI"` indicates that we are applying the Pauli X operator on qubit 0, Z operator on qubit 1, and the identity operator on qubit 2.
    - The `circuit` is a simple 3-qubit quantum circuit, with a Hadamard gate on the first qubit followed by a CNOT gate.

    The function will compute the expectation value of the given observables for the state produced by the circuit and return the result as a floating-point number.
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


def _random_circuit(
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

    This function takes a string representation of a quantum circuit and converts it
    into a Qibo `models.Circuit` object. The input string should represent the circuit in a
    "text drawing" format, where qubits are represented by lines and gates are represented
    by various symbols. The function processes these text representations and constructs
    the equivalent quantum circuit object.

    :param text_draw:
        A string representing the quantum circuit in a text drawing format, where
        qubits are represented by lines and gates are represented by symbols.
    :param parameters:
        An optional list of parameters to be passed to the gates, defaults to None.
    :return:
        A Qibo `models.Circuit` object corresponding to the input text drawing.

    :example:

    .. code-block:: python

        from qibo.models import Circuit
        from mymodule import draw_to_circuit

        # Define the text representation of the quantum circuit
        text_rep = '''
        q0: ─H─S─X───SDG───H───
        q1: ─X─S───o─X─T─o─S───
        q2: ───────Z─T───|───o─
        q3: ─────────────|─X─|─
        q4: ──RX─────RX──Z───Z─
        '''

        # Convert the text drawing into a quantum circuit
        circuit = draw_to_circuit(text_rep)

        # Print the resulting circuit
        print(circuit.draw())
    """
    split = text_draw.splitlines()
    qubits_lst = []
    split = [element for element in split if element.strip()]
    split = [element for element in split if element != ""]
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
                        for i in range(step, num_steps):
                            if qubit_state[i + 1] == "─":
                                tmp = tmp + qubit_state[i]

                                gate_name = tmp
                                qubits = idx
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
                                    param = parameters[idx][parameter_tracker][1]
                                    # Pass qubits and parameters if the gate requires both
                                    circuit.add(gate_class(qubits, param))
                                    parameter_tracker += 1
                                    break
                                else:
                                    circuit.add(gate_class(qubits))
                                    break

                            else:
                                tmp = tmp + qubit_state[i]

                elif char == "o":
                    saved_qubit.append(idx)

        for idx in saved_qubit:
            # if list_multiple_gates[idx]:
            '''print(
                "Add gate: ",
                list_multiple_gates[idx][0][0],
                " qubit ",
                list_multiple_gates[idx][0][1],
            )'''
            circuit.add(
                getattr(gates, list_multiple_gates[idx][0][0])(
                    *list_multiple_gates[idx][0][1]
                )
            )
            list_multiple_gates[idx].remove(list_multiple_gates[idx][0])
    return circuit


def parse_qsim(fname: str, depth:int = 22):
        """Produces a Qibo :class:`Circuit` based on a .qsim description of a circuit.

        This function parses a .qsim file to create a quantum circuit based on the gate instructions
        provided in the file. The .qsim file should contain gate operations for the quantum circuit.
        The file should start with the number of qubits, followed by a series of lines specifying
        gate operations in the format: `<moment> <gate> <qubit> [additional parameters]`. The function
        also allows setting a `depth` parameter to limit the number of moments to process.

        :param fname: Path to the .qsim file containing the quantum circuit description.
        :param depth: Maximum number of moments (steps) to process from the .qsim file. Default is 22.

        :return: A :class:`models.Circuit` object from Qibo containing the parsed gate operations.

        :raises ValueError: If the first line of the .qsim file cannot be converted to an integer
                            representing the number of qubits.

        :example:

        Suppose you have a `.qsim` file with the following content:

        .. code-block:: text

            3
            0 rz 0 0.5
            1 x_1_2 1
            2 rz 2 1.0
            3 hz_1_2 0
            4 fs 1 2

        You can parse this file and create a circuit like this:

        .. code-block:: python

            from mymodule import parse_qsim

            # Parse the .qsim file to create a quantum circuit
            circuit = parse_qsim("path_to_qsim_file.qsim", depth=10)

            # Print the resulting circuit
            print(circuit.draw())
        """

        with open(fname, 'r') as f:
            data = f.read()
        print(data)
        lines = data.strip().splitlines()

        try:
            qcount = int(lines.pop(0))
        except ValueError:
            raise ValueError('First line should be qubit count')

        c = models.Circuit(qcount)

        for l in lines:
            l = l.strip()
            if l == '': continue
            gdesc = l.split(' ')

            q = int(gdesc[2])
            moment = int(gdesc[0])
            if moment >= depth:
                break

            random_param = np.random.uniform(0, 2*np.pi)
            if gdesc[1] == 'rz':
                phase = float(gdesc[3]) / np.pi
                c.add(gates.RZ(q, phase))
            elif gdesc[1] == 'hz_1_2':
                c.add(gates.H(q))
            elif gdesc[1] == 'x_1_2':
                c.add(gates.RX(q,theta=random_param))
            elif gdesc[1] == 'y_1_2':
                c.add(gates.RY(q, theta=random_param))
            elif gdesc[1] == 'fs':
                q1 = int(gdesc[3])
                #theta = float(gdesc[4]) / np.pi
                #phi = float(gdesc[5]) / np.pi
                c.add(gates.CZ(q, q1))


        #c.gates = gates
        return c
