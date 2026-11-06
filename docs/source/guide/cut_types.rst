Cut types
=========

The notebook :doc:`../notebooks/10_cut_types` runs everything on this page. The examples
below use this circuit:

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(6)
    circuit.ry(0.6, range(6))
    for q in range(5):
        circuit.cz(q, q + 1)
    circuit.rx(0.3, range(6))
    circuit.ry(0.2, range(6))
    circuit.cz(1, 2)

Measure-and-prepare gate cut
----------------------------

:func:`~Qdislib.gate_cutting` cuts a two-qubit gate with a measure-and-prepare
quasiprobability decomposition. A gate other than CZ or CX is cut as the ``rzz`` rotations
it decomposes into: one for ``cp``, ``crz``, ``cry``, ``rxx``, ``ryy``, ``rzz``, and up to
three for a general gate. Each is cut at its own angle (Mitarai and Fujii), with 6 terms
and ``gamma = 1 + 2|sin(phi)|`` for ``rzz(phi)``, so a weak gate is cheap: 3 for a
CZ-like angle, 1.2 for ``rzz(0.1)``.

A side of the cut has a mid-circuit measurement. A simulator defers it onto an extra qubit
read at the end, so it samples the final state once instead of running every shot.

.. code-block:: python

    value = qd.gate_cutting(circuit, ["CZ_9"], "Z" * 6, shots=8192, seed=1)
    print(value, qd.analytical_solution(circuit, "Z" * 6))

Hadamard gate cut
-----------------

:func:`~Qdislib.gate_cutting_hadamard` (Harrow and Lowe,
`arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) replaces the mid-circuit
measurement with one local ancilla per fragment running a Hadamard test. It cuts any
two-qubit gate as a single cut. For one rotation its ``gamma`` is that of the
measure-and-prepare cut, which needs no ancilla and so is the better choice. For a gate with
several interactions the rotations' ``gamma``\ s multiply (``iswap`` 9, ``swap`` 27) while
the Hadamard cut has 7. Qdislib warns when it is clearly cheaper, and says how to ask for it.

.. note::

   Cutting a gate other than CZ, CX, ``cp`` or ``rzz`` with the Hadamard test needs Qibo installed, also for a Qiskit circuit (``pip install "qdislib[qibo]"``). A dense ``unitary`` gate cannot be cut with it in this release; use ``gate_cutting``.

With ``joint=True``, the default, the gates that cross the **same boundary** are cut
together: one shared ancilla per fragment runs a single Hadamard test over all of them
(Schmitt, Piveteau and Sutter, Quantum 9, 1634, 2025). For ``k`` CZ or CX gates that is
``gamma = 2**(k+1) - 1`` instead of ``3**k`` (7 instead of 9 at ``k = 2``, 15 instead of 27
at ``k = 3``), ``2**k (2**(k+1) - 1)`` settings instead of ``6**k``, and one ancilla per
side instead of ``k``. The gates may sit at different times and be any two-qubit gate; no
classical communication is needed.

.. code-block:: python

    joint = qd.gate_cutting_hadamard(circuit, ["CZ_8", "CZ_24"], "Z" * 6, shots=None)
    separate = qd.gate_cutting_hadamard(circuit, ["CZ_8", "CZ_24"], "Z" * 6, shots=None, joint=False)
    print(joint, separate, qd.cutting_overhead(circuit, {"hadamard": ["CZ_8", "CZ_24"]}))

The shared ancilla has to stay coherent from the boundary's first cut gate to its last, and
needs a controlled gate with each of them. On hardware with short T2 or sparse connectivity,
cutting gates far apart in time separately can be better.

``gamma`` bounds the shots per sample: the bound is reached when every sampled configuration
gets a single shot, where the joint cut needs ``(3**k / (2**(k+1)-1))**2`` times fewer
samples (3.2 at ``k = 3``). With exact or many-shot fragment values the actual variance of
either method is well below its bound and depends on the circuit.

Wire cut
--------

:func:`~Qdislib.wire_cutting` cuts a qubit's line between two gates, with a measure-and-prepare
decomposition of 8 terms and ``gamma`` 4.

With ``optimized=True`` the wires that leave the same fragment are measured together:
``2**k + 1`` measurement circuits for ``k`` wires instead of ``3**k`` (5 instead of 9 at
``k = 2``, 9 instead of 27 at ``k = 3``). Each circuit reads a commuting group of Paulis, so
every Pauli is read once. The value is the same. At equal total shots the error is usually lower
(always in tests with 3 wires), but with 2 wires per-wire measurement is sometimes better. Some
of the joint circuits add CZ gates between the fragment's cut qubits; ``joint=False`` measures
each wire on its own, for hardware with costly two-qubit gates. The preparation side (``4**k``)
does not change, and in general a wire cut's sampling overhead cannot drop below ``4**n``
without classical communication between fragments.

.. code-block:: python

    wire = [("CZ_9", "CZ_10")]
    value = qd.wire_cutting(circuit, wire, "Z" * 6, shots=8192, seed=1, optimized=True)
    per_wire = qd.wire_cutting(circuit, wire, "Z" * 6, shots=8192, seed=1, optimized=True, joint=False)
    print(value, per_wire)

**Golden cut points.** When the fragment that sends a wire has a real state (every gate a real
matrix such as ``ry``, ``cz``, ``cx``, ``h``, no gate cut, and every wire it receives golden
too) and sends no other wire, the wire's two ``Y`` terms are zero and are left out: 6 terms
instead of 8, ``gamma`` 3 instead of 4, and no ``Y`` measurement or ``|+i>`` preparation (Chen et
al., "golden cutting points"). It is exact and automatic, and ``cutting_overhead`` reports it. A
``Y`` in the observable, a complex gate (``rx``, ``rz``, ``s``, ...) or a second wire from the same
fragment keeps all 8.

Any mix of cuts
---------------

:func:`~Qdislib.circuit_cutting` takes measure-and-prepare gate cuts, Hadamard gate cuts and wire
cuts on the same circuit, in any number. It is the engine the other three functions run on, so a
cut costs and behaves the same through any of them. ``cuts`` is a dict with any of ``"gate"``,
``"hadamard"`` and ``"wire"``, which is what ``find_cut(..., mixed=True)`` returns, or a list of
gate names and wire edges, with ``hadamard=True`` to cut its gates with the Hadamard test.

.. code-block:: python

    cuts = {"gate": ["CZ_7"], "hadamard": ["CZ_8", "CZ_24"], "wire": [("CZ_10", "CZ_11")]}
    print(qd.cutting_overhead(circuit, cuts))
    print(qd.cutting_overhead(circuit, cuts, optimized=True))

    value = qd.circuit_cutting(circuit, cuts, "Z" * 6, shots=None)
    dist = qd.circuit_cutting_sampling(circuit, cuts, shots=None, qubits=[0, 5])
    print(value, qd.analytical_solution(circuit, "Z" * 6), dist)

``gate_wire_cutting`` from 1.x still works, with a deprecation warning.

How the cost grows
------------------

``gamma`` is the product of the cuts': 3 per measure-and-prepare CZ, 3 per Hadamard CZ or
``2**(k+1) - 1`` for a joint group of ``k``, 4 per wire. The number of circuits depends on
``optimized``:

- ``optimized=False`` runs every fragment for every term of the combined decomposition:
  ``6**g * 8**w`` configurations, times the fragments.
- ``optimized=True`` runs each fragment once per setting of *the cuts that touch it* and reuses
  it in every term. A product only builds up inside a fragment; across fragments the counts add.
  Per cut and side: 5 circuits for a measure-and-prepare gate, since two of its 6 terms share
  each side's circuit; 3 per wire sent (``2**k + 1`` for ``k`` sent jointly); 4 per wire
  received.

Circuits for ``k`` cuts between two fragments:

===  ===========  ==========  ===========  ============  ==========  =============
k    gate, std    gate, opt   Hadamard,    wire, std     wire, opt   wire, opt +
                              joint                                  joint
===  ===========  ==========  ===========  ============  ==========  =============
1    12           10          12           16            7           7
2    72           50          56           128           25          21
3    432          250         240          1,024         91          73
4    2,592        1,250       992          8,192         337         273
===  ===========  ==========  ===========  ============  ==========  =============

Three fragments, ``g`` gate cuts between the first two and ``w`` wire cuts between the last
two. The middle fragment holds both, so its settings multiply:

=====  =========  ===========  =============  ==================  ============
g, w   std        opt          opt + joint    Hadamard joint +    gamma
                                              wire joint          (mp / joint)
=====  =========  ===========  =============  ==================  ============
1, 1   144        24           24             28                  12 / 12
2, 2   6,912      266          166            184                 144 / 112
3, 3   331,776    3,564        1,314          1,264               1,728 / 960
=====  =========  ===========  =============  ==================  ============

Every cut has to separate
-------------------------

A cut is valid only if its two sides end up in different fragments, judged with all the other
cuts applied. Two gates across one boundary, or a gate cut that only comes loose together with a
wire cut, are valid. A cut that another, uncut path bypasses separates nothing, since summing its
terms in one fragment just rebuilds what was cut, and raises ``InvalidCutError`` naming it, even
when the other cuts do split the circuit. A cut set that splits nothing at all also raises an
error. An empty cut list runs the circuit as it is.
