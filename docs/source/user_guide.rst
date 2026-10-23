User guide
==========

This guide covers the basics of using Qdislib, as well as details on the main algorithms included in the library for circuit partitioning, simulation, and distributed execution.

Qdislib supports *gate cutting* and *wire cutting* techniques to decompose large quantum circuits and execute them on distributed and heterogeneous environments. All computations are parallelized transparently using `PyCOMPSs <https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar/>`_.

Example
-------

Below is a basic example demonstrating how to construct a quantum circuit and apply gate cutting using Qdislib:

.. code-block:: python

    from qibo import models, gates
    import Qdislib as qd  # qd is Qdislib

    # Build a 10-qubit circuit
    def entire_circuit():
        c = models.Circuit(10)
        c.add(gates.H(0))
        c.add(gates.RX(1, 0.3))
        c.add(gates.RY(2, 0.4))
        c.add(gates.RZ(3, 0.5))
        c.add(gates.CZ(0, 1))
        c.add(gates.CZ(1, 2))
        c.add(gates.CZ(2, 3))
        c.add(gates.CZ(3, 4))
        c.add(gates.CZ(4, 5))
        c.add(gates.CZ(5, 6))
        return c

    circuit = entire_circuit()

    # Find a gate to cut
    cut = qd.find_cut(circuit)

    # Apply gate cutting
    result = qd.gate_cutting(circuit, cut)
    print("Reconstructed expectation value:", result)

It is worth noting that, although the code above looks completely sequential,
all Qdislib algorithms and operations are internally parallelized using PyCOMPSs.

How to run Qdislib
------------------

Qdislib can be installed and used as a regular Python library, and runs
**serially** out of the box (``python my_Qdislib_application.py``). PyCOMPSs is
optional: when it is installed, Qdislib transparently parallelizes its tasks,
and you launch applications with the `runcompss` or `enqueue_compss` tools:

.. code-block:: bash

    runcompss my_Qdislib_application.py

Alternatively, inside Jupyter notebooks, you can use the `ipycompss` interface:

.. code-block:: python

    import ipycompss
    ipycompss.start(project="project.xml", resources="resources.xml", monitor=True, graph=True)

Refer to the :doc:`quickstart guide <quickstart>` and :doc:`examples <examples>` for details on PyCOMPSs setup and configuration.

Calling convention
------------------

Every cut type -- gate, wire, Hadamard, and any mix of them (``circuit_cutting``) --
has the same four functions with the same arguments:

.. code-block:: python

    X(circuit, cuts, observable=None, *, backend=None, shots=1024, seed=None, sync=True, ...)
    X_sampling(circuit, cuts, *, backend=None, shots=1024, seed=None, sync=True, ...)
    X_subcircuits(circuit, cuts, observable=None, *, software=None, ...)
    X_sampling_subcircuits(circuit, cuts, *, software=None, ...)

where ``X`` is ``gate_cutting``, ``wire_cutting``, ``gate_cutting_hadamard`` or
``circuit_cutting``.

- ``X`` returns the expectation value of ``observable`` (a Pauli string or a
  weighted Pauli sum; ``None`` is ``Z`` on every qubit); ``X_sampling`` returns the
  output **distribution** ``{bitstring: probability}``.
- ``backend`` is **where** the subcircuits run: an SDK name (``"qibo"``,
  ``"qiskit"``, ``"cudaq"``, ``"pennylane"``, ``"qilisdk"`` -- its local simulator), a registered
  backend name (``"IBM_BACKEND"``, ``"IQM_BACKEND"``, ``"BSC_QPU_BACKEND"``, ...) or a
  backend object from ``qd.build_backend(...)``. ``None`` uses the circuit's own SDK.
- ``shots=None`` gives the **exact** result (statevector) on every cut type.
- ``seed`` makes simulator runs reproducible (each subcircuit gets its own stream);
  nothing is random with ``shots=None``.
- ``sync=False`` returns a future (PyCOMPSs), to overlap several cuts.
- Every option has the same name and default everywhere: ``optimized=False``
  (gate, wire, ``circuit_cutting``), ``joint=True`` (wire, Hadamard,
  ``circuit_cutting``), ``num_samples=None`` (every cut type), ``hadamard=False``
  (``circuit_cutting``).
- Only ``circuit``, ``cuts`` and ``observable`` go by position; the rest by keyword.

The three-step forms take ``software`` -- the **language** the subcircuits are
written in (``None``: the circuit's own) -- and return a list of
:class:`~Qdislib.Subcircuit` objects, the same for every cut type, plus
``reconstruct``. See :ref:`three-step <three-step-label>`.

``qd.cutting_overhead(circuit, cuts, ...)`` reports the cost of any cut without
running it: ``{"subcircuits", "configurations", "fragments", "gamma"}``.

The 1.x names (``observables=``, ``software=`` on the all-in-one calls, ``qpu=``,
``qpu_dict=``, ``gpu=``, ``method=``, ``gates_cut=``, ``gate_wire_cutting``,
``cut_cost`` and the standalone ``*_reconstruction`` functions) still work in this
release, with a ``DeprecationWarning`` naming the replacement.

Circuits, SDKs and naming the gates to cut
------------------------------------------

A circuit can come from **Qibo**, **Qiskit**, **CUDA-Q**, **PennyLane** or
**QiliSDK** (``pip install qdislib[qilisdk]``); subcircuits can be written in any of
them (``software=``) and run on any of them (``backend=``), whatever the input was.

Inside, every circuit becomes the same DAG, whose gates are named and defined as in
**OpenQASM 3** (``stdgates.inc``: ``h``, ``x``, ``rz``, ``p``, ``u``, ``cx``, ``cz``,
``cp``, ``cu``, ``swap``, ``ccx``, ... plus ``rxx``, ``ryy``, ``rzz``, ``iswap``).
Each SDK's conventions are converted on the way in -- qibo's ``CU3`` (a controlled
SU(2) gate) is ``cu`` with the matching phase -- so the same circuit gives the same
DAG in every SDK.

A cut names a gate by its **name and position** in the circuit, counting from 1:
``"CZ_8"`` is the 8th instruction, a CZ. Any SDK's spelling of the gate works --
``"CX_7"``, ``"CNOT_7"`` and ``"cx_7"`` are the same gate -- and naming a gate that
is not there says what gate 7 actually is. :func:`~Qdislib.find_cut` returns these
names, and :func:`~Qdislib.draw_cut` prints each gate's index on the diagram. A wire
cut is the edge between two gates on the same qubit, ``("RY_3", "CZ_8")``.

Another SDK can be added in one place (``Qdislib/utils/sdks.py``): a converter to
and from the DAG, and a backend, registered with
``Qdislib.utils.sdks.register_sdk``.

Algorithms
----------

Qdislib includes several circuit cutting strategies:

- **Gate Cutting** (``gate_cutting``): Cuts selected two-qubit gates — a measure-and-prepare quasiprobability decomposition — to divide a circuit into smaller, independently executable subcircuits, then post-processes them to reconstruct the original expectation value. Cuts **any** two-qubit gate: ``CZ``/``CX`` natively, and other gates by automatically decomposing them into ``CZ`` + single-qubit gates. Cost is ``6`` terms per cut.

- **Hadamard-test Gate Cutting** (``gate_cutting_hadamard``): An alternative gate cut (Harrow & Lowe, `arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) that replaces the mid-circuit measurement with **one local ancilla** per fragment running a Hadamard test. It cuts **any** two-qubit gate as a *single* cut at the **optimal** overhead (``gamma(CZ) = 3``, tending to ``1`` for weak interactions), which makes it a good fit for superconducting hardware where mid-circuit measurement is costly. Like the other cuts, it has an all-in-one call and a three-step ``_subcircuits`` workflow for expectation values, and ``_sampling`` variants for the output distribution. ``num_samples=N`` draws the configurations at random instead of enumerating them all (the native Harrow–Lowe estimator), which bounds the cost for many cuts.

  With ``joint=True``, the gates that cross the **same boundary** are cut together:
  one shared ancilla per fragment runs a single Hadamard test over all of them
  (Schmitt, Piveteau & Sutter, Quantum 9, 1634, 2025). For ``k`` CZ/CX gates that is
  ``gamma = 2**(k+1) - 1`` instead of ``3**k`` (7 instead of 9 at ``k = 2``, 15
  instead of 27 at ``k = 3``), ``2**k (2**(k+1) - 1)`` settings instead of ``6**k``,
  and one ancilla per side instead of ``k``. The gates may sit at different times and
  be any two-qubit gate; no mid-circuit measurement or classical communication is
  needed. It is the default (``joint=False`` cuts each gate on its own);
  ``qd.cutting_overhead(circuit, {"hadamard": cut})`` returns ``gamma``.

  .. code-block:: python

      value = qd.gate_cutting_hadamard(circuit, ["CZ_5", "CZ_9"], "ZZZZ")
      separate = qd.gate_cutting_hadamard(circuit, ["CZ_5", "CZ_9"], "ZZZZ", joint=False)

  The shared ancilla must stay coherent from the boundary's first cut gate to its
  last, and needs a controlled gate with each of them, so on hardware with short T2
  or sparse connectivity, cutting gates far apart in time separately can be better.
  ``gamma`` bounds the shots per sample: the bound is reached when every sampled
  configuration gets a single shot, where the joint cut needs ``(3**k / (2**(k+1)-1))**2``
  times fewer samples (3.2x at ``k = 3``). With exact or many-shot fragment values the
  actual variance of either method is well below its bound and depends on the circuit.

- **Wire Cutting** (``wire_cutting``): Cuts "wires" (qubit connections) between gates to disconnect logical qubits, via a measure-and-prepare decomposition (``8`` terms per cut, or ``3**n + 4**n`` with ``optimized=True``).

  With ``optimized=True``, the cut wires that leave the same fragment are measured
  together: ``2**k + 1`` measurement circuits for ``k`` such wires instead of ``3**k``
  (5 instead of 9 at ``k = 2``, 9 instead of 27 at ``k = 3``). Each circuit reads a
  whole group of commuting Paulis of those wires, so every Pauli is read exactly
  once. Where the wires arrive does not matter, and it applies to
  ``wire_cutting``, ``wire_cutting_subcircuits`` and ``wire_cutting_sampling``. The
  value is the same; at equal total shots the error is usually lower (always lower
  in tests with 3 wires), but with 2 wires per-wire measurement is sometimes
  better. Some of the joint circuits add ``CZ`` gates between the fragment's cut
  qubits; pass ``joint=False`` to measure each wire on its own, with no entangling
  gates, e.g. on hardware with costly two-qubit gates. The preparation side
  (``4**k``) is unchanged, and a wire cut's sampling overhead cannot drop below
  ``4**n`` without classical communication between fragments.

  .. code-block:: python

      value = qd.wire_cutting(circuit, cut, "ZZZZ", optimized=True)
      per_wire = qd.wire_cutting(circuit, cut, "ZZZZ", optimized=True, joint=False)

- **Any combination of cuts** (``circuit_cutting``): measure-and-prepare gate cuts,
  Hadamard gate cuts and wire cuts on the **same circuit**, in any number. Their
  quasiprobability decompositions tensor, and every fragment runs in parallel. It has
  the same four forms as the other cuts: ``circuit_cutting`` (expectation value),
  ``circuit_cutting_subcircuits`` (three-step), ``circuit_cutting_sampling`` (output
  distribution) and ``circuit_cutting_sampling_subcircuits``.

  .. code-block:: python

      cuts = {"gate": ["CX_7"], "hadamard": ["CZ_3"], "wire": [("U_5", "CZ_9")]}
      value = qd.circuit_cutting(circuit, cuts, "ZZZZZ", optimized=True)
      dist = qd.circuit_cutting_sampling(circuit, cuts, optimized=True)
      qd.cutting_overhead(circuit, cuts, optimized=True)
      # e.g. {'subcircuits': 81, 'configurations': 2304, 'fragments': 4, 'gamma': 144.0}

  ``cuts`` is a dict (any subset of ``"gate"``, ``"hadamard"``, ``"wire"``, as
  ``find_cut(..., mixed=True)`` returns) or a list of gate names and wire edges, with
  ``hadamard=True`` to cut its gates with the Hadamard test. ``joint`` (default
  ``True``) covers both the joint wire measurement and the joint Hadamard cut.
  ``shots=None`` gives the exact result, mid-circuit measurements included.
  ``gate_wire_cutting`` (1.x) still works, with a deprecation warning.

  **How the cost grows.** The sampling overhead ``gamma`` is the product of the cuts'
  (3 per measure-and-prepare gate, 3 per Hadamard gate or ``2**(k+1) - 1`` for a joint
  group of ``k``, 4 per wire), and the shots grow as ``gamma**2``. The number of
  circuits depends on ``optimized``:

  - ``optimized=False`` runs every fragment for every term of the combined
    decomposition: ``6**g * 8**w * ...`` configurations, times the fragments.
  - ``optimized=True`` runs each fragment once per setting of *the cuts that touch
    it* and reuses it in every term. A product only builds up inside a fragment; across
    fragments the counts add. Per cut and side: 5 circuits for a measure-and-prepare
    gate (two of its 6 terms share each side's circuit; a general gate counts as the
    CZs it is decomposed into), 3 per wire sent (``2**k + 1`` for ``k`` sent
    jointly), 4 per wire received.

  **Every cut has to separate.** In every algorithm, a cut is valid only if its two
  sides end up in different fragments, judged with all the other cuts applied: two
  gates across the same boundary, or a gate cut that only comes loose together with
  a wire cut, are valid. A cut that another, uncut path bypasses separates nothing
  (summing its terms in one fragment just rebuilds what was cut) and raises
  ``InvalidCutError`` naming it, even when the other cuts do split the circuit.

  ``k`` cuts between two fragments:

  ===  ===========  ==========  ===========  ============  ==========  =============
  k    gate, std    gate, opt   Hadamard,    wire, std     wire, opt   wire, opt +
                                joint                                  joint
  ===  ===========  ==========  ===========  ============  ==========  =============
  1    12           10          12           16            7           7
  2    72           50          56           128           25          21
  3    432          250         240          1,024         91          73
  4    2,592        1,250       992          8,192         337         273
  ===  ===========  ==========  ===========  ============  ==========  =============

  Three fragments, ``g`` gate cuts between the first two and ``w`` wire cuts between
  the last two (the middle fragment holds both, so its settings multiply):

  =====  =========  ===========  =============  ==================  ==========
  g, w   std        opt          opt + joint    Hadamard joint +    gamma
                                                wire joint          (mp / joint)
  =====  =========  ===========  =============  ==================  ==========
  1, 1   144        24           24             28                  12 / 12
  2, 2   6,912      266          166            184                 144 / 112
  3, 3   331,776    3,564        1,314          1,264               1,728 / 960
  =====  =========  ===========  =============  ==================  ==========

Every subcircuit is an independent PyCOMPSs task, so all methods fan out and run in parallel; for large fan-outs, ``sync=False`` overlaps independent cuts and ``reduce_chunk=n`` folds the results through a balanced reduction tree.


**Find Cut Function**

To simplify the cutting process, Qdislib provides the `find_cut` utility to automatically suggest suitable locations for gate or wire cutting:

.. code-block:: python

    cut = qd.find_cut(circuit)                       # Gate and wire cutting suggestion
    wire_cut = qd.find_cut(circuit, gate_cut=False)  # Wire cutting suggestion
    gate_cut = qd.find_cut(circuit, wire_cut=False)  # Gate cutting suggestion

The returned list contains candidate gate names or gate pairs that are suitable cut points.

``find_cut`` also accepts a few constraints and knobs:

.. code-block:: python

    cut = qd.find_cut(circuit, max_qubits=5, max_cuts=2)  # bound subcircuit size / #cuts
    cut = qd.find_cut(circuit, seed=0)                    # reproducible (same seed -> same cut)
    cut = qd.find_cut(circuit, forbidden_gates=["CZ_3"])  # never cut these gates
    cut = qd.find_cut(circuit, forbidden_qubits=[0, 1])   # never cut on these qubits

**Hardware-aware cut selection (SparseCut).** When you pass the target QPU's
``coupling_map``, ``find_cut`` switches to a hardware-aware mode: it cuts only gates
whose two qubits are not adjacent on the device, preferring the **longest-routing**
ones — the gates that would otherwise need the most SWAPs (and so cause the most
crosstalk / error) — and returns the cheapest such gate cut that **separates** the
circuit (``[]`` if those gates cannot split it). This follows the *SparseCut* idea of
Guo et al., `arXiv:2511.05492 <https://arxiv.org/abs/2511.05492>`_.

.. code-block:: python

    line = [(0, 1), (1, 2), (2, 3), (3, 4)]           # the QPU's connectivity
    cut = qd.find_cut(circuit, coupling_map=line, max_cuts=2)
    cut = qd.find_cut(circuit, coupling_map=line, layout={0: 2, 1: 0, 2: 4})

The result is a gate cut you feed to ``gate_cutting`` as usual; pass a ``layout``
(logical → physical qubit map) if your circuit is not placed identity-first.

**Cost-aware cut discovery.** With ``mixed=True``, ``find_cut`` searches every kind of
cut at once and prices it the way :func:`~Qdislib.circuit_cutting` runs it, minimising
the sampling overhead ``gamma`` (shots grow as ``gamma**2``):

- each gate cut takes the cheaper of measure-and-prepare (``gamma = 3`` per CZ the gate
  decomposes into) and the Hadamard test, whose ``gamma`` follows the gate's angle —
  ``1 + 2|sin θ|``-like for a rotation, so weak gates are nearly free;
- the Hadamard cuts between two fragments are priced jointly (``2**(k+1) - 1`` for k
  CZ/CX, one shared ancilla per side);
- a wire cut (``gamma = 4``) splits a qubit's timeline where that is cheaper.

Every fragment fits ``max_qubits`` *with its Hadamard ancillas* (``None`` asks for
halves), every cut separates its two sides, and the search returns the best cut it
found within ``time_limit`` seconds (10 by default). For up to 40 qubits it also seeds
the search with the exact optimal gate partition (an integer program, ``exact=``).

.. code-block:: python

    cut = qd.find_cut(circuit, max_qubits=12, mixed=True)  # {"gate", "hadamard", "wire"}
    qd.cutting_overhead(circuit, cut)                       # gamma, circuits
    val = qd.circuit_cutting(circuit, cut, "ZZZZ")

``hadamard=False`` keeps to measure-and-prepare gate cuts, ``wire_cut=False`` /
``gate_cut=False`` switch a kind off, and ``joint=False`` prices Hadamard cuts one by
one (as ``circuit_cutting(..., joint=False)`` runs them).

**Seeing a cut.** ``qd.draw_cut(circuit, cuts, filename=None)`` draws the circuit as a
diagram with each gate's index (the name to cut), every fragment the cut leaves in its
own colour, and every cut marked: a measure-and-prepare gate cut in red, a Hadamard cut
in purple, a wire cut as a ✂ on the wire. It takes any cut ``circuit_cutting`` takes,
and draws a cut that does not separate with ✗ rather than raising, naming it in the
title with the fragment count and ``gamma``:

.. code-block:: python

    qd.draw_cut(circuit, qd.find_cut(circuit, max_qubits=4, mixed=True), "cut.png")

.. tip::

   Each cut multiplies the number of subcircuit evaluations — a **gate cut by 6**,
   a **wire cut by 8** — so the cost grows *exponentially* in the number of cuts
   (``optimized=True`` turns much of that product into a sum, see above). Preview it
   before running with ``qd.cutting_overhead``:

   .. code-block:: python

       qd.cutting_overhead(circuit, {"gate": cut})
       # {'subcircuits': 72, 'configurations': 36, 'fragments': 2, 'gamma': 9.0}

.. _three-step-label:

**Recommended workflow: generate → evaluate → reconstruct**

The recommended way to use Qdislib is the **three-step workflow**: generate the
subcircuits, evaluate each one in your own loop, then recombine the results.
Owning the evaluation loop is the point — that is where you dispatch to a
simulator, a GPU, or a real QPU, batch, or cache, while Qdislib handles the
cutting and the reconstruction. This enables:

- Custom execution of subcircuits (e.g., submitting to specific quantum devices).
- Saving/loading subcircuits for delayed simulation.
- Manual post-processing outside Qdislib.

.. code-block:: python

    # 1. Cut the circuit into subcircuits, written in the SDK you choose
    #    ("qibo", "qiskit", "cudaq", "pennylane", "qilisdk"; None = the circuit's own).
    subcircuits = qd.gate_cutting_subcircuits(circuit, gate_cut, "ZZZZ", software="qiskit")
    # (the same for wire_cutting_, gate_cutting_hadamard_ and circuit_cutting_subcircuits)

    # 2. Run each subcircuit — this loop is yours to control.
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

    # 3. Recombine into the original circuit's expectation value.
    result = subcircuits.reconstruct(values)

Every item is a :class:`~Qdislib.Subcircuit`, whatever the cut type:

- ``s.circuit`` -- the native circuit in ``software`` (a Qiskit ``QuantumCircuit``
  above), to inspect, draw, transpile or submit yourself;
- ``s.software``, ``s.num_qubits``;
- ``s.observable`` -- the Z/I string read at its end, or ``None`` when every qubit is
  read in Z (the ``_sampling_subcircuits`` lists).

Run the items of a ``*_subcircuits`` list with ``qd.expectation_value(s, shots=...,
backend=..., seed=...)`` and those of a ``*_sampling_subcircuits`` list with
``qd.sample_subcircuit(s, shots=..., backend=..., seed=...)``; ``backend`` takes the
same values as in the all-in-one calls, and ``None`` runs the subcircuit on its own
SDK's simulator. Choosing it per subcircuit is how you route by size:

.. code-block:: python

    ibm = qd.build_backend("IBM_BACKEND", qd.BackendConfig(...))
    values = [
        qd.expectation_value(s, shots=4000, backend=ibm if s.num_qubits <= 127 else None)
        for s in subcircuits
    ]

``subcircuits.reconstruct(results)`` is the one way to recombine (the 1.x
``*_reconstruction`` functions are deprecated).

If you prefer a single call, ``qd.gate_cutting(...)`` / ``qd.wire_cutting(...)``
run the cut, execution and reconstruction together. This modular structure gives
users fine-grained control over hybrid quantum-classical workflows.


Observables
-----------

The ``observable`` argument accepts a per-qubit Pauli string over
``X``, ``Y``, ``Z`` and ``I`` (identity), e.g. ``"ZIZ"``. ``X`` and ``Y`` terms
are measured by rotating into the Z basis automatically, on any backend.

You can also pass a **weighted sum of Pauli operators** (a Hamiltonian),
evaluated as ``<sum_i c_i P_i> = sum_i c_i <P_i>`` — either as a list of
``(coefficient, pauli_string)`` pairs or as a ``{pauli_string: coefficient}``
mapping:

.. code-block:: python

    # <0.5 ZZ> + <1.0 XI> - <0.3 IZ>
    hamiltonian = [(0.5, "ZZ"), (1.0, "XI"), (-0.3, "IZ")]
    result = qd.gate_cutting(circuit, cut, hamiltonian)

    # equivalently, as a mapping
    result = qd.gate_cutting(circuit, cut, {"ZZ": 0.5, "XI": 1.0, "IZ": -0.3})

The same forms work with every cut type (and its ``_subcircuits`` variant) and
``qd.analytical_solution``.


Drawing configurations instead of running them all (``num_samples``)
--------------------------------------------------------------------

A cut is a weighted sum over configurations -- ``6**n`` for ``n`` gate cuts, ``8**n``
for ``n`` wire cuts, products of the Hadamard settings. By default every
configuration runs with ``shots`` each. With ``num_samples=N``, Qdislib instead uses
``N`` configurations in proportion to ``|weight| / gamma`` and runs only the fragment
circuits they use, each with ``shots`` per draw that uses it -- the quasiprobability
Monte-Carlo estimate. It is unbiased.

The draws are **allocated**, not drawn at random (Dai & Koczor: a share in proportion
to the weight is never worse than random draws). When all configurations fit in
``N``, each gets its share, which is the same as running them all; otherwise the cuts
that fit are allocated and the rest drawn.

The error has two parts. Which configurations are drawn contributes at most
``gamma / sqrt(N)``, **whatever the shots**; the shot noise is the rest. When the
configurations do not all fit in ``N``, the first part usually dominates, so more
shots per draw do not help: for a given budget, use **more draws with fewer shots
each** (``shots=1`` is fine). When they all fit, the error falls with the shots, as
when running every configuration.

.. code-block:: python

    value = qd.gate_cutting(circuit, cut, "ZZZZ", shots=1, num_samples=100_000)
    value = qd.circuit_cutting(circuit, cuts, "ZZZZ", shots=1, num_samples=100_000, seed=3)
    dist = qd.wire_cutting_sampling(circuit, wire_cut, shots=1, num_samples=200_000)

    # three-step: each subcircuit's ``draws`` says how many draws use it
    subs = qd.circuit_cutting_subcircuits(circuit, cuts, "ZZZZ", num_samples=50_000)
    values = [qd.expectation_value(s, shots=s.draws) for s in subs]
    value = subs.reconstruct(values)

When to use it:

- **too many configurations to list** -- e.g. a dozen gate cuts (``6**12``) or a
  mix of gate and wire cuts: the cost is set by ``N``, not by the number of cuts;
- **unequal weights** (Hadamard cuts of rotations): the shots go where the weight
  is, which can save orders of magnitude over equal shots per configuration.

With few enough configurations to run them all, ``num_samples`` at least as large as
their number gives the same accuracy as running them all with the same total shots
(each gets its share). Sampling does not change ``gamma``: the draws and shots for a
given precision still grow as ``gamma**2``.

Sampling the output distribution
--------------------------------

``gate_cutting`` / ``wire_cutting`` / ``gate_cutting_hadamard`` return a single
**expectation value**. To reconstruct the full **output distribution**
``{bitstring: probability}`` instead — e.g. to read out an image or inspect a state —
use the ``_sampling`` variant of the same cut, which has the same all-in-one and
three-step forms:

.. code-block:: python

    dist = qd.gate_cutting_sampling(circuit, cut, shots=200_000)
    # or  qd.wire_cutting_sampling(circuit, wire_cut, shots=200_000)
    # or  qd.gate_cutting_hadamard_sampling(circuit, cut, shots=200_000)
    # or  qd.circuit_cutting_sampling(circuit, cuts, shots=200_000, optimized=True)

    # three-step: sample each subcircuit yourself, then reconstruct
    subs = qd.gate_cutting_hadamard_sampling_subcircuits(circuit, cut)
    dist = subs.reconstruct([qd.sample_subcircuit(s, shots=8192) for s in subs])

For the Hadamard cut, every fragment is read in the Z basis anyway, so the
distribution costs the same circuits as one expectation value. With ``shots=None``
every cut type gives the exact distribution.

Reconstructing a distribution carries a quasiprobability overhead, so it needs
**more shots** than an expectation value. Pass ``seed=`` for a reproducible
distribution:

.. code-block:: python

    dist = qd.gate_cutting_sampling(circuit, cut, shots=200_000, seed=0)


Caching repeated subcircuits
----------------------------

Cutting produces many *functionally equivalent* subcircuits (repeated fragments,
or the same ansatz across variational iterations). The **semantic cache** keys
each subcircuit by a ZX-reduced graph hash and simulates an equivalent fragment
only once. Install it with ``pip install "qdislib[cache]"`` and pass a store to
any cutting call:

.. code-block:: python

    cache = qd.default_circuit_cache("memory")   # or "lmdb" (on disk) / "redis"
    value = qd.gate_cutting(circuit, cut, cache=cache)
    print(cache.stats())                         # hits, misses, hit-rate, size

Bound the in-memory store for long runs: ``max_entries`` evicts the
least-recently-used entry, and ``ttl`` expires entries after a number of seconds
(useful after an SDK upgrade):

.. code-block:: python

    cache = qd.default_circuit_cache("memory", max_entries=5000, ttl=3600)

Under PyCOMPSs, use the ``"redis"`` store so the cache is shared across workers.


Choosing a backend
-------------------

The ``backend`` argument selects where the subcircuits run -- an SDK's simulator
(``"qibo"``, ``"qiskit"``, ``"cudaq"``, ``"pennylane"``, ``"qilisdk"``), a registered backend by name,
or a configured backend object -- and is accepted by every cutting call and by
``expectation_value`` / ``sample_subcircuit``. Under the hood a pluggable **backend
layer** turns it into a concrete backend:

.. code-block:: python

    qd.available_backends()                 # registered backends
    qd.build_execution_backend("qiskit")    # the single place that chooses one
    qd.QdBackEnd(qpu_dict={...})            # load-balance across a pool of devices

You can register your own backend (an in-house simulator or a cloud QPU) by
subclassing :class:`~Qdislib.core.backend.base_backend.BaseBackend` and calling
``qd.register_backend(...)``; the cutting code then selects it via
``backend="MY_BACKEND"`` with no other change. See notebook *08 · Backends under
the hood*.


Resource Allocation
-------------------

Qdislib is designed to run efficiently across CPUs, GPUs, and quantum devices depending on the user’s configuration and available resources. By default, simulation is performed on CPU backends.

.. _gpu-support-label:

Using GPUs
----------

Qdislib runs subcircuits on Qiskit Aer's GPU simulator (``device="GPU"``,
cuStateVec-enabled) through a GPU-configured Aer backend; send only the large
subcircuits there from the three-step loop:

.. code-block:: python

    gpu = qd.build_backend("AER_BACKEND", qd.BackendConfig(extra={"use_gpu": True}))
    value = qd.gate_cutting(circuit, cut, backend=gpu)

    subs = qd.gate_cutting_subcircuits(circuit, cut, software="qiskit")
    values = [
        qd.expectation_value(s, backend=gpu if s.num_qubits >= 18 else None)
        for s in subs
    ]

GPU usage is transparent: backend selection and memory transfers are handled
internally. The GPU path supports the same cases as the CPU path, including an
identity observable on a cut qubit.

Running on real QPUs
--------------------

Pass the device as ``backend`` -- by its registered name, or as a configured
backend object:

.. code-block:: python

    # BSC (MareNostrum Ona) via the shared-filesystem bridge:
    value = qd.gate_cutting(circuit, cut, backend="BSC_QPU_BACKEND")

    # IBM Quantum via qiskit-ibm-runtime; it connects once, here (token, channel,
    # instance in extra= or the IBM_QUANTUM_* environment variables):
    ibm = qd.build_backend(
        "IBM_BACKEND",
        qd.BackendConfig(num_qubits=156, extra={"backend_name": "ibm_marrakesh"}),
    )
    value = qd.gate_cutting(circuit, cut, backend=ibm)

Every subcircuit then runs on that device. To send only the subcircuits that fit
(and simulate the rest), choose the backend per subcircuit in the three-step loop
(see above). A real device returns only aggregated counts, so an identity
observable on a cut qubit that keeps evolving *after* the cut is not supported on a
QPU (it is on a simulator).
