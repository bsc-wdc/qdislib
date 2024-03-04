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

import numpy as np
import qibo
from qibo import models, gates, hamiltonians  # , callbacks
from qibo.symbols import Z, X, Y, I


from Qdislib.circuit_classes import NewCircuit


# Now we define the two subcircuits of the above example
@task(returns=NewCircuit)
def first_subcircuit(circuit, basis, numgates, sub_circuit_1_dimension):
    """
    Input:
    - basis (str): ['X', 'Y' or 'Z']. Basis in which we will measure the cut qubit
    - shots (int): Number of times that we obtain a state for doing stadistics

    Output:
    - expectation_value (float): expectation value of the three qubits when we measure ZxZxO_i, being O_i the 'basis' selected
    """

    # We introduce a quantum gate that is not inherently supported by the qibo framework.
    S_dagger = np.array([[1, 0], [0, -1j]])

    sub_circuit_1 = models.Circuit(sub_circuit_1_dimension)

    for i in range(0, numgates):
        sub_circuit_1.add(circuit.queue[i])

    if basis == "X":
        sub_circuit_1.add(gates.H(sub_circuit_1_dimension - 1))
    elif basis == "Y":
        sub_circuit_1.add(gates.H(sub_circuit_1_dimension - 1))
        sub_circuit_1.add(gates.Unitary(S_dagger, (sub_circuit_1_dimension - 1)))
    elif basis in ["Z", "I"]:
        pass
    else:
        raise Exception("Error - Unsupported Basis")

    # measurement gates
    sub_circuit_1.add(gates.M(*range(sub_circuit_1_dimension)))

    new_circuit = NewCircuit(sub_circuit_1)

    return new_circuit


@task(returns=NewCircuit)
def second_subcircuit(circuit, initial, numgates, sub_circuit_1_dimension):
    """
    Input:
    - initial (str): initial state that the first qubit will be initialized in. ['0', '1', '+', '-', '+i', '-i']
    - shots (int): Number of times that we obtain a state for doing stadistics

    Output:
    - expectation_value (float): expectation value of the three qubits when we measure ZxZ when we initialize the first qubit in 'initial'

    """

    sub_circuit_2 = models.Circuit(circuit.nqubits - sub_circuit_1_dimension + 1)

    # We prepare the different initial states
    if initial == "1":
        sub_circuit_2.add(gates.X(0))
    elif initial == "+":
        sub_circuit_2.add(gates.H(0))

    elif initial == "-":
        sub_circuit_2.add(gates.X(0))
        sub_circuit_2.add(gates.H(0))

    elif initial == "+i":
        sub_circuit_2.add(gates.H(0))
        sub_circuit_2.add(gates.S(0))

    elif initial == "-i":
        sub_circuit_2.add(gates.X(0))
        sub_circuit_2.add(gates.H(0))
        sub_circuit_2.add(gates.S(0))

    for j in range(numgates, len(circuit.queue)):
        circuit_copy = circuit.copy(True)
        temp = circuit_copy.queue[j]

        if len(temp.qubits) > 1:
            control = temp.qubits[0] - sub_circuit_1_dimension + 1
            temp._set_control_qubits((control,))
            target = temp.qubits[1] - sub_circuit_1_dimension + 1
            temp._set_target_qubits((target,))
        else:
            target = temp.qubits[0] - sub_circuit_1_dimension + 1
            temp._set_target_qubits((target,))
        sub_circuit_2.add(temp)

    sub_circuit_2.add(gates.M(*range(circuit.nqubits - sub_circuit_1_dimension + 1)))

    new_circuit_2 = NewCircuit(sub_circuit_2)

    return new_circuit_2


@task(returns=int)
def compute_expectation_value(freq, basis, shots):
    """This function computes the expectation value given a probability distribution (the output of the quantum computer)
    in a given basis that we choose.
     INPUT:
      - freq (dict): frequency distribution coming from the quantum computer.
      - basis (str):   we aim to compute the expectation value of this set of operators.  For example, "XYY" indicates that we
        calculate the expectation value of X over the first qubit, and the expectation value of Y over the second and third qubits.
      - shots (int): Numer of times that we have runed the quantum computer, needed to compute the probability in the probability distribution.

      OUTPUT:
       - expectation_value (float): Final expectation value.

       This function assumes that for computing the 'X' and the 'Y' expectation value, the qubit state it is in the appropiate diagonal basis.

       For the moment, we only implement two types of cases for the basis:
       1) Only combinations of X, Y or/and Z. For example: 'XXYYXZ'
       2) Only a single I operator in the last position. For example 'ZYXXYI'.

    """

    expectation_value = 0
    for key, value in freq.items():
        if basis.count("I") == 0:
            ones = key.count("1")
            if ones % 2 == 0:
                expectation_value += float(value) / shots
            else:
                expectation_value -= float(value) / shots

        if basis[-1] == "I":
            ones = key[:2].count(("1"))
            if ones % 2 == 0:
                expectation_value += float(value) / shots
            else:
                expectation_value -= float(value) / shots

    return expectation_value


@task(returns=NewCircuit)
def first_subcircuit_basis(circuit1, basis, qubit):
    """
    Input:
    - basis (str): ['X', 'Y' or 'Z']. Basis in which we will measure the cut qubit
    - shots (int): Number of times that we obtain a state for doing stadistics

    Output:
    - expectation_value (float): expectation value of the three qubits when we measure ZxZxO_i, being O_i the 'basis' selected
    """

    # We introduce a quantum gate that is not inherently supported by the qibo framework.
    S_dagger = np.array([[1, 0], [0, -1j]])

    dimension = circuit1.nqubits

    if basis == "X":
        circuit1.add(gates.H(qubit))
    elif basis == "Y":
        circuit1.add(gates.H(qubit))
        circuit1.add(gates.Unitary(S_dagger, (qubit)))
    elif basis in ["Z", "I"]:
        pass
    else:
        raise Exception("Error - Unsupported Basis")

    # measurement gates
    circuit1.add(gates.M(*range(dimension)))

    new_circuit = NewCircuit(circuit1)

    return new_circuit


@task(returns=NewCircuit)
def second_subcircuit_states(circuit2, initial, qubit):
    """
    Input:
    - initial (str): initial state that the first qubit will be initialized in. ['0', '1', '+', '-', '+i', '-i']
    - shots (int): Number of times that we obtain a state for doing stadistics

    Output:
    - expectation_value (float): expectation value of the three qubits when we measure ZxZ when we initialize the first qubit in 'initial'

    """

    dimension = circuit2.nqubits

    # We prepare the different initial states
    if initial == "1":
        circuit2.queue.insert(0, gates.X(qubit))
        # circuit2.add(gates.X(0))
    elif initial == "+":
        # circuit2.add(gates.H(0))
        circuit2.queue.insert(0, gates.H(qubit))

    elif initial == "-":
        circuit2.queue.insert(0, gates.H(qubit))
        circuit2.queue.insert(0, gates.X(qubit))
        # circuit2.add(gates.X(0))
        # circuit2.add(gates.H(0))

    elif initial == "+i":
        circuit2.queue.insert(0, gates.S(qubit))
        circuit2.queue.insert(0, gates.H(qubit))
        # circuit2.add(gates.H(0))
        # circuit2.add(gates.S(0))

    elif initial == "-i":
        circuit2.queue.insert(0, gates.S(qubit))
        circuit2.queue.insert(0, gates.H(qubit))
        circuit2.queue.insert(0, gates.X(qubit))
        # circuit2.add(gates.X(0))
        # circuit2.add(gates.H(0))
        # circuit2.add(gates.S(0))

    circuit2.add(gates.M(*range(dimension)))

    new_circuit_2 = NewCircuit(circuit2)

    return new_circuit_2
