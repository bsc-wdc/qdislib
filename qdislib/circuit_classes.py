"""
Circuit cutting classes.

This file contains all auxiliary classes that wrap the qibo Circuit objects.
"""

from pycompss.api.task import task


class NewCircuitResult:
    """CircuitResult (qibo.states.CircuitResult) class wrapper.

    Contains an internal CircuitResult object instead of extending the CircuitResult class.
    """

    __slots__ = ["result"]

    def __init__(self, result=None):
        """Circuit result wrapper constructor.

        :param result: Qibo result object.
        """
        self.result = result

    @task(returns=dict)
    def frequencies_COMPSs(self, binary=True, registers=False):
        """Calculate the frequencies function task wrapper.

        :param binary: Binary.
        :param registers: Registers.
        :return: Frequencies.
        """
        return dict(self.result.frequencies(binary, registers))


class NewCircuit:
    """Circuit class wrapper.

    Contains an internal Circuit object instead of extending the Circuit class.
    """

    __slots__ = ["circuit"]

    def __init__(self, circuit=None):
        """Circuit wrapper constructor.

        :param circuit: Qibo circuit object.
        """
        self.circuit = circuit

    @task(returns=NewCircuitResult)
    def execute_COMPSs(self, initial_state=None, nshots=None):
        """Execute function task wrapper.

        :param initial_state: Initial state.
        :param nshots: Number of shots.
        :return: Circuit results.
        """
        result = self.circuit.execute(initial_state, nshots)
        new_result = NewCircuitResult(result)
        return new_result
    
    @task(returns=NewCircuitResult)
    def execute_qc_COMPSs(self, connection, initial_state=None, nshots=None):
        """Execute function task wrapper.

        :param initial_state: Initial state.
        :param nshots: Number of shots.
        :return: Circuit results.
        """
        job_ids = connection.execute(self, initial_state, nshots)
        result = connection.get_results(job_ids=job_ids)
        new_result = NewCircuitResult(result)
        return new_result