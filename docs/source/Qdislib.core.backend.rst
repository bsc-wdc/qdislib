Qdislib.core.backend package
============================

The pluggable execution layer. Every target -- a simulator, a GPU, or a real QPU --
is a :class:`~Qdislib.core.backend.base_backend.BaseBackend`; the cutting algorithms
never branch on the backend, they hand each subcircuit to the object that
:func:`~Qdislib.core.backend.registry.build_execution_backend` builds from the
``backend=`` argument. A new backend needs one method: see :doc:`guide/backends`.

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

Qdislib.core.backend.aer\_backend module
----------------------------------------

.. automodule:: Qdislib.core.backend.aer_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.qibo\_backend module
-----------------------------------------

.. automodule:: Qdislib.core.backend.qibo_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.cudaq\_backend module
------------------------------------------

.. automodule:: Qdislib.core.backend.cudaq_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.pennylane\_backend module
----------------------------------------------

.. automodule:: Qdislib.core.backend.pennylane_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.qilisdk\_backend module
--------------------------------------------

.. automodule:: Qdislib.core.backend.qilisdk_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.iqm\_backend module
----------------------------------------

.. automodule:: Qdislib.core.backend.iqm_backend
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.backend.qpu\_backends module
-----------------------------------------

.. automodule:: Qdislib.core.backend.qpu_backends
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: Qdislib.core.backend
   :members:
   :undoc-members:
   :show-inheritance:
