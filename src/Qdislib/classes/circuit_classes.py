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

"""
Circuit cutting classes.

This file contains all auxiliary classes that wrap the qibo Circuit objects.
"""

from pycompss.api.task import task


class _NewCircuitResult:
    """CircuitResult (qibo.states.CircuitResult) class wrapper.

    Contains an internal CircuitResult object instead of extending the
    CircuitResult class.
    """

    __slots__ = ["result"]

    def __init__(self, result=None):
        """Circuit result wrapper constructor.

        :param result: Qibo result object.
        """
        self.result = result

    @task(returns=dict)
    def frequencies_compss(self, binary=True, registers=False):
        """Calculate the frequencies function task wrapper.

        :param binary: Binary.
        :param registers: Registers.
        :return: Frequencies.
        """
        return dict(self.result.frequencies(binary, registers))


class _NewCircuit:
    """Circuit class wrapper.

    Contains an internal Circuit object instead of extending the Circuit class.
    """

    __slots__ = ["circuit"]

    def __init__(self, circuit=None):
        """Circuit wrapper constructor.

        :param circuit: Qibo circuit object.
        """
        self.circuit = circuit

    @task(returns=_NewCircuitResult)
    def execute_compss(self, initial_state=None, nshots=None):
        """Execute function task wrapper.

        :param initial_state: Initial state.
        :param nshots: Number of shots.
        :return: Circuit results.
        """
        result = self.circuit.execute(initial_state, nshots)
        new_result = _NewCircuitResult(result)
        return new_result

    @task(returns=int)
    def execute_qc_compss(self, connection, initial_state=None, nshots=None):
        """Execute function task wrapper.

        :param initial_state: Initial state.
        :param nshots: Number of shots.
        :return: Circuit results.
        """
        job_ids = connection.execute(self.circuit, initial_state, nshots)
        # result = connection.get_results(job_ids=job_ids)
        # new_result = _NewCircuitResult(result)
        return job_ids
