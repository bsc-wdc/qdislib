Finding cuts
============

The notebook :doc:`../notebooks/03_finding_cuts` runs everything on this page.

The default search
------------------

:func:`~Qdislib.find_cut` scores several graph-partitioning strategies (Kernighan-Lin,
Girvan-Newman, spectral, METIS) and keeps the best. It returns a list of gate names, or of wire
edges when it searches only wires, which you pass to the cutting functions.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(8)
    circuit.ry(0.6, range(8))
    for layer in range(2):
        for q in range(0, 7, 2):
            circuit.cz(q, q + 1)
        for q in range(1, 7, 2):
            circuit.cz(q, q + 1)
    circuit.ry(0.4, range(8))

    print(qd.find_cut(circuit, max_qubits=4, wire_cut=False))        # gate cuts only
    print(qd.find_cut(circuit, max_qubits=4, gate_cut=False)[:2])    # wire cuts only
    print(qd.find_cut(circuit, max_qubits=4, max_cuts=2))            # at most 2 cuts
    print(qd.find_cut(circuit, max_qubits=4, seed=0))                # reproducible
    print(qd.find_cut(circuit, max_qubits=4, forbidden_gates=["CZ_12"]))   # never these gates
    print(qd.find_cut(circuit, max_qubits=4, forbidden_qubits=[0, 1]))     # never on these qubits

The limits are checked as each cut will be applied: every fragment of a wire cut must fit in
``max_qubits`` qubits, and no cut touches a ``forbidden_qubits`` wire. ``max_cuts=0`` allows no
cut at all, and ``max_qubits`` must be at least 1. ``[]`` means no cut meets the limits, or that
the circuit's qubits already split into independent groups and need none. The other way round,
the cutting functions raise ``ValueError`` for a cut that does not split the circuit.

``top_n=`` returns several candidates with their scores, ``weights=`` retunes the loss, and
``qd.set_log_level("DEBUG")`` prints the candidates and the strategy that won.

Searching by cost: ``mixed=True``
---------------------------------

With ``mixed=True``, ``find_cut`` minimises ``gamma`` over every kind of cut, priced as
:func:`~Qdislib.circuit_cutting` runs it:

- each gate cut takes the cheaper of measure-and-prepare and the Hadamard test; both follow the
  gate's angle, and on a tie measure-and-prepare wins, as it needs no ancilla;
- the Hadamard cuts between two fragments are priced jointly (``2**(k+1) - 1`` for ``k`` CZ or CX,
  one shared ancilla per side);
- a wire cut (``gamma`` 4, 3 at a golden point) splits a qubit's timeline where that is cheaper.

Every fragment fits ``max_qubits`` *with its Hadamard ancillas* (``None`` asks for halves), every
cut separates its two sides, and the search returns the best cut it found within ``time_limit``
seconds, 10 by default. Up to 40 qubits it also seeds the search with the exact optimal gate
partition, an integer program solved with SciPy's HiGHS (``exact=``). It returns the dict
``{"gate", "hadamard", "wire"}``:

.. code-block:: python

    cut = qd.find_cut(circuit, max_qubits=5, mixed=True)
    print(cut, qd.cutting_overhead(circuit, cut))
    print(qd.circuit_cutting(circuit, cut, "Z" * 8, shots=None))

``hadamard=False`` keeps to measure-and-prepare gate cuts, ``wire_cut=False`` or ``gate_cut=False``
switch a kind off, and ``joint=False`` prices Hadamard cuts one by one, as
``circuit_cutting(..., joint=False)`` runs them.

The changelog reports what this gave on QAOA circuits, from ``gamma`` 729 to about 32 at 10 qubits,
and from ``10**13.4`` to ``10**4.4`` at 100 qubits in two halves, with the search taking about 1-8 s at
100 qubits. On circuits where measure-and-prepare is already optimal it returns the same cut.

For a given QPU: ``coupling_map``
---------------------------------

A gate between qubits that are not neighbours on the device needs SWAPs. With ``coupling_map``,
``find_cut`` cuts only such gates, preferring the longest-routing ones, and returns the cheapest
gate cut of them that **separates** the circuit, or ``[]`` if those gates cannot split it. This
follows the SparseCut idea of Guo et al.,
`arXiv:2511.05492 <https://arxiv.org/abs/2511.05492>`_. ``layout`` is the logical-to-physical map
if your circuit is not placed identity-first.

.. code-block:: python

    line = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]

    remote = QuantumCircuit(6)
    remote.ry(0.5, range(6))
    for a, b in [(0, 1), (1, 2), (3, 4), (4, 5)]:
        remote.cz(a, b)
    remote.cz(0, 5)
    remote.cz(1, 4)

    print(qd.find_cut(remote, coupling_map=line))
    print(qd.find_cut(remote, coupling_map=line, layout={0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5}))

Seeing a cut
------------

:func:`~Qdislib.draw_cut` draws the circuit with each gate's index, every fragment in its own
colour and every cut marked: a measure-and-prepare gate cut in green, a Hadamard cut in purple, a wire
cut in red, as scissors on the wire. It takes any cut ``circuit_cutting`` takes. A cut that does not
separate is drawn with a cross and named in the title rather than raising. The title gives the
number of fragments and ``gamma``. Pass a filename to save it, or an ``ax=`` to draw into.

.. code-block:: python

    qd.draw_cut(circuit, qd.find_cut(circuit, max_qubits=4, mixed=True), "cut.png")

Before you run: ``cutting_overhead``
------------------------------------

:func:`~Qdislib.cutting_overhead` reports what any cut costs without running it, and
``qd.shots_for`` the shots for a target error (see :doc:`shots`).

.. code-block:: python

    print(qd.cutting_overhead(circuit, qd.find_cut(circuit, max_qubits=4, wire_cut=False)))

``qd.estimate(circuit, observable, max_qubits=..., shots=...)`` finds a cut and runs it in one call.
