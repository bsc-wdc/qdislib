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
import qibo
from qibo import gates

from Qdislib.classes.circuit_classes import _NewCircuit


@task(returns=int)
def _compute_expectation_value(
    freq, basis, shots, wire_observables=False, b=None, gpu=False, gpu_counter=0
):
    """
    Compute the expectation value given a probability
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
    if gpu:
        qibo.set_device(f"/GPU:{gpu_counter}")

    if wire_observables:
        for key, value in basis.items():
            if value == "-":
                basis[key] = b

    basis = "".join([value for key, value in sorted(basis.items())])
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


@task(returns=_NewCircuit)
def _first_subcircuit_basis(circuit1, basis, qubit):
    """
    First subcircuit basis.

    INPUT:
     - basis (str): ['X', 'Y' or 'Z']. Basis in which we will measure the
                    cut qubit.
     - shots (int): Number of times that we obtain a state for doing
                    stadistics.

    OUTPUT:
     - expectation_value (float): expectation value of the three qubits when
                                  we measure ZxZxO_i, being O_i the 'basis'
                                  selected.

    """

    # We introduce a quantum gate that is not inherently supported by
    # the qibo framework.
    # s_dagger = np.array([[1, 0], [0, -1j]])

    dimension = circuit1.nqubits

    if basis == "X":
        circuit1.add(gates.H(qubit))
    elif basis == "Y":
        circuit1.add(gates.H(qubit))
        circuit1.add(gates.SDG(qubit))
    elif basis in ["Z", "I"]:
        pass
    else:
        raise ValueError("Error - Unsupported Basis")

    # measurement gates
    circuit1.add(gates.M(*range(dimension)))

    new_circuit = _NewCircuit(circuit1)

    return new_circuit


@task(returns=_NewCircuit)
def _second_subcircuit_states(circuit2, initial, qubit):
    """
    Second subcircuit states.

    INPUT:
     - initial (str): initial state that the first qubit will be initialized
                      in ['0', '1', '+', '-', '+i', '-i'].
     - shots (int): Number of times that we obtain a state for doing
                    stadistics.

    OUTPUT:
     - expectation_value (float): expectation value of the three qubits when
                                  we measure ZxZ when we initialize the first
                                  qubit in 'initial'.

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

    new_circuit_2 = _NewCircuit(circuit2)

    return new_circuit_2
