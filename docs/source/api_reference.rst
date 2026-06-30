API Reference
=============

Qdislib.api: Quantum Distributed Library API
--------------------------------------------


This module provides the primary interface for interacting with the Quantum Distributed Library.


Qdislib.core.cutting_algorithms: Circuit Cutting
------------------------------------------------

:class:`cutting_algorithms.gate_cutting <Qdislib.core.cutting_algorithms.gate_cutting>`
- Distributed Gate cutting algorithm.

:class:`cutting_algorithms.wire_cutting <Qdislib.core.cutting_algorithms.wire_cutting>`
- Distributed Wire cutting algorithm.


Qdislib.core.graph_algorithms: Graph Algorithms
-----------------------------------------------

:class:`graph_algorithms.random_graph_generator <Qdislib.core.graph_algorithms.random_graph_generator>`
- Random graph generation utility.


Qdislib.core.find_cut: Find Cut Circuit Cutting
-------------------------------------------------

:class:`find_cut.find_cut <Qdislib.core.find_cut.find_cut>`
- Quantum Circuit find cut algorithms.


Qdislib.utils: Utility functions
--------------------------------

:meth:`utils.backends <Qdislib.utils.backends>`
- Backend-agnostic helpers, including ``expectation_value`` and ``CutSubcircuits``.

:meth:`utils.circuit <Qdislib.utils.circuit>`
- Circuit object and methods.

:meth:`utils.exceptions <Qdislib.utils.exceptions>`
- Qdislib exceptions.

:meth:`utils.graph_qibo <Qdislib.utils.graph_qibo>`
- Qibo graph object and methods.

:meth:`utils.graph_qiskit <Qdislib.utils.graph_qiskit>`
- Qiskit graph object and methods.

:meth:`utils.graph_cudaq <Qdislib.utils.graph_cudaq>`
- CUDA-Q graph object and methods.
