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
import networkx as nx
from qibo import models, gates
import numpy as np
import networkx as nx

from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.api import compss_barrier
from pycompss.api.parameter import *

from Qdislib.core.optimal_cut.optimal_cut import optimal_cut, execute_optimal_cut
from Qdislib.core.cutting_algorithms._pycompss_functions import _compute_expectation_value

class QAOASolver(object):
    def __init__(self, input_graph: nx.Graph, num_shots=8000, seed=10):
        self.G = input_graph
        self.E = list(input_graph.edges())
        self.V = list(input_graph.nodes())
        np.random.seed(seed)
        self.gamma = np.random.uniform()
        np.random.seed(seed*2)
        self.beta = np.random.uniform()
        self.num_shots = num_shots
        self.qaoa_solution = 0

def create_qaoa_circuit(G, p: int=1, seed=10):
    """
    Creates QAOA circuit self.qaoa_circuit for graph self.G of p layers
    Args:
        p: int, the number of layers of QAOA circuit
    Return:
        None
    """
    qaoa_solver = QAOASolver(G, seed)

    qaoa_circuit = models.Circuit(len(qaoa_solver.V))

    for i in range(0, qaoa_circuit.nqubits):
        qaoa_circuit.add(gates.H(i))

    for _ in range(p):
        for e in qaoa_solver.E:
            qaoa_circuit.add(gates.H(e[1]))
            qaoa_circuit.add(gates.CZ(e[0], e[1]))
            qaoa_circuit.add(gates.H(e[1]))
            qaoa_circuit.add(gates.RZ(e[1], qaoa_solver.gamma))
            qaoa_circuit.add(gates.H(e[1]))
            qaoa_circuit.add(gates.CZ(e[0], e[1]))
            qaoa_circuit.add(gates.H(e[1]))

        for i in range(0,qaoa_circuit.nqubits):
            qaoa_circuit.add(gates.RX(i, 2 * qaoa_solver.beta))
    return qaoa_circuit

def execute_qaoa_expected_value(qaoa_circuit, cut, observables, num_shots=8000,chunk=1,verbose=False,sync=True,gpu=False,gpu_counter=0):
    reconstruction = execute_optimal_cut(observables=observables,circuit=qaoa_circuit,cut=cut,shots=num_shots,chunk=chunk,verbose=verbose,sync=sync,gpu=gpu,gpu_counter=gpu_counter)
    return reconstruction

def execute_qaoa(qaoa_circuit, num_shots=8000):
    qaoa_circuit.add(gates.M(*range(qaoa_circuit.nqubits)))
    result = qaoa_circuit(nshots=num_shots)
    qaoa_counts = dict(result.frequencies(binary=True))
    return qaoa_counts

def find_maxcut(qaoa_solver, qaoa_counts):
    for s in qaoa_counts:
        qaoa_solver.qaoa_solution = max(qaoa_solver.qaoa_solution, _compute_cut(qaoa_solver, s))
    return qaoa_solver.qaoa_solution

def _compute_cut(qaoa_solver, s: str):
    """
    Computes the cut value of cut/partition string s (little-endian) of graph self.G
    Args:
        s: str, the string of cut/partition of nodes in little-endian order
    Return:
        cut_value: int
    """

    cut_value = 0
    s = s[::-1]

    for e in qaoa_solver.E:
        if s[e[0]] != s[e[1]]:
            cut_value += 1
    
    return cut_value

def optimize_parameters(self):
    pass

def solve_maxcut(G, p=1, num_shots=8000, seed=10):
    qaoa_solver = QAOASolver(G, num_shots, seed)
    qaoa_circuit = create_qaoa_circuit(qaoa_solver, p)
    qaoa_counts = execute_qaoa(qaoa_circuit, num_shots)
    qaoa_solution = find_maxcut(qaoa_solver, qaoa_counts)
    return qaoa_solution

def solve_maxcut_expected_value(G, observables, p=1, num_shots=8000, max_qubits=None,gate_cut=True,wire_cut=True,draw=False, seed=10, chunk=1,verbose=False,sync=True,gpu=False,gpu_counter=0):
    if not wire_cut and not gate_cut:
        qaoa_circuit = create_qaoa_circuit(G, p, seed)
        qaoa_circuit.add(gates.M(*range(qaoa_circuit.nqubits)))
        freq = execute_qaoa(qaoa_circuit,num_shots)
        return _compute_expectation_value(freq, observables, num_shots)
    else:
        qaoa_circuit = create_qaoa_circuit(G, p, seed)
        cut = optimal_cut(circuit=qaoa_circuit,max_qubits=max_qubits,gate_cut=gate_cut,wire_cut=wire_cut,draw=draw)
        qaoa_reconstruction = execute_qaoa_expected_value(qaoa_circuit=qaoa_circuit, cut=cut, observables=observables, num_shots=num_shots, chunk=chunk,verbose=verbose,sync=sync,gpu=gpu,gpu_counter=gpu_counter)
        return qaoa_reconstruction