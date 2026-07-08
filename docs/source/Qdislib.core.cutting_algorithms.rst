Qdislib.core.cutting\_algorithms package
========================================

Submodules
----------

Qdislib.core.cutting\_algorithms.gate\_cutting module
-----------------------------------------------------

.. automodule:: Qdislib.core.cutting_algorithms.gate_cutting
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.cutting\_algorithms.wire\_cutting module
-----------------------------------------------------

.. automodule:: Qdislib.core.cutting_algorithms.wire_cutting
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.cutting\_algorithms.hadamard\_cut module
-----------------------------------------------------

Ancilla / double-Hadamard-test gate cutting (Harrow & Lowe, `arXiv:2403.01018
<https://arxiv.org/abs/2403.01018>`_): cuts **any** two-qubit gate with one local
ancilla per fragment and no mid-circuit measurement, at the optimal overhead.
Provides :func:`~Qdislib.core.cutting_algorithms.hadamard_cut.gate_cutting_hadamard`
(all-in-one), the ``_subcircuits`` three-step form, and the Monte-Carlo
``_sampling`` variant.

.. automodule:: Qdislib.core.cutting_algorithms.hadamard_cut
   :members:
   :undoc-members:
   :show-inheritance:

Qdislib.core.cutting\_algorithms.mixed\_cutting module
------------------------------------------------------

Mixed wire + gate cutting (circuit knitting): cut some gates *and* some wires on the
same circuit in one integrated, fully parallel workflow via
:func:`~Qdislib.core.cutting_algorithms.mixed_cutting.gate_wire_cutting`.

.. automodule:: Qdislib.core.cutting_algorithms.mixed_cutting
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: Qdislib.core.cutting_algorithms
   :members:
   :undoc-members:
   :show-inheritance:
