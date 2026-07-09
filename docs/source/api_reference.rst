API Reference
=============

Qdislib.api: Quantum Distributed Library API
--------------------------------------------


This module provides the primary interface for interacting with the Quantum Distributed Library.


Qdislib.core.cutting_algorithms: Circuit Cutting
------------------------------------------------

:class:`cutting_algorithms.gate_cutting <Qdislib.core.cutting_algorithms.gate_cutting>`
    Distributed gate cutting (measure-and-prepare); cuts any two-qubit gate.

:class:`cutting_algorithms.wire_cutting <Qdislib.core.cutting_algorithms.wire_cutting>`
    Distributed wire cutting algorithm.

:class:`cutting_algorithms.hadamard_cut <Qdislib.core.cutting_algorithms.hadamard_cut>`
    Ancilla / double-Hadamard-test gate cut (arXiv:2403.01018): any gate, optimal
    overhead, no mid-circuit measurement (``gate_cutting_hadamard`` + sampling / three-step).

:class:`cutting_algorithms.mixed_cutting <Qdislib.core.cutting_algorithms.mixed_cutting>`
    Mixed wire + gate cutting on one circuit (``gate_wire_cutting``).


Qdislib.core.graph_algorithms: Graph Algorithms
-----------------------------------------------

:class:`graph_algorithms.random_graph_generator <Qdislib.core.graph_algorithms.random_graph_generator>`
    Random graph generation utility.


Qdislib.core.find_cut: Find Cut Circuit Cutting
-------------------------------------------------

:class:`find_cut.find_cut <Qdislib.core.find_cut.find_cut>`
    Quantum circuit find-cut algorithms (options: ``max_qubits``, ``max_cuts``, ``seed``, ``forbidden_gates``, ``forbidden_qubits``).


Qdislib.core.estimate: One-call estimation
------------------------------------------

:class:`estimate <Qdislib.core.estimate>`
    ``estimate`` (cut, run and reconstruct in one call) and ``cut_cost`` (preview a cut's subcircuit count / sampling overhead before running).


Qdislib.core.backend: Execution backends
-----------------------------------------

:class:`backend.base_backend <Qdislib.core.backend.base_backend>`
    ``BaseBackend`` — the interface every execution target (simulator, GPU, QPU) implements.

:class:`backend.registry <Qdislib.core.backend.registry>`
    ``build_execution_backend`` / ``register_backend`` / ``available_backends``.

:class:`backend.unified <Qdislib.core.backend.unified>`
    ``QdBackEnd`` — the load-balancing front-end over a pool of backends.


Qdislib.utils: Utility functions
--------------------------------

:meth:`utils.backends <Qdislib.utils.backends>`
    Backend-agnostic helpers, including ``expectation_value`` and ``CutSubcircuits``.

:meth:`utils.circuit_cache <Qdislib.utils.circuit_cache>`
    Semantic subcircuit cache: ``default_circuit_cache``, ``semantic_key``.

:meth:`utils.circuit <Qdislib.utils.circuit>`
    Circuit object and methods.

:meth:`utils.exceptions <Qdislib.utils.exceptions>`
    Qdislib exceptions.

:meth:`utils.graph_qibo <Qdislib.utils.graph_qibo>`
    Qibo graph object and methods.

:meth:`utils.graph_qiskit <Qdislib.utils.graph_qiskit>`
    Qiskit graph object and methods.

:meth:`utils.graph_cudaq <Qdislib.utils.graph_cudaq>`
    CUDA-Q graph object and methods.
