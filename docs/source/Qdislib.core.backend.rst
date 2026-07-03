Qdislib.core.backend package
============================

The pluggable execution layer. Every target -- a CPU simulator (Qibo / Qiskit /
CUDA-Q), a GPU, or a real QPU -- is a :class:`~Qdislib.core.backend.base_backend.BaseBackend`;
the cutting algorithms never branch on the backend, they hand each subcircuit to
the object that :func:`~Qdislib.core.backend.registry.build_execution_backend`
selects from the ``software`` / ``gpu`` / ``qpu`` arguments.

Submodules
----------

Qdislib.core.backend.base\_backend module
-----------------------------------------

.. automodule:: Qdislib.core.backend.base_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.registry module
-------------------------------------

.. automodule:: Qdislib.core.backend.registry
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.unified module
-----------------------------------

.. automodule:: Qdislib.core.backend.unified
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.circuit\_converter module
----------------------------------------------

.. automodule:: Qdislib.core.backend.circuit_converter
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: Qdislib.core.backend
   :members:
   :undoc-members:
   :show-inheritance:
