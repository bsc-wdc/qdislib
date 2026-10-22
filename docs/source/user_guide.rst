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
  needed. ``qd.hadamard_cut_overhead(circuit, cut, joint=True)`` returns the overhead
  per boundary.

  .. code-block:: python

      value = qd.gate_cutting_hadamard(circuit, ["CZ_5", "CZ_9"], "ZZZZ", joint=True)

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

      value = qd.wire_cutting(circuit, cut, observables="ZZZZ", optimized=True)
      per_wire = qd.wire_cutting(circuit, cut, observables="ZZZZ", optimized=True, joint=False)

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
      qd.circuit_cutting_overhead(circuit, cuts, optimized=True)
      # e.g. {'subcircuits': 81, 'configurations': 2304, 'fragments': 4, 'gamma': 144.0}

  ``cuts`` is a dict (any subset of ``"gate"``, ``"hadamard"``, ``"wire"``, as
  ``find_cut(..., mixed=True)`` returns) or a list of gate names and wire edges, with
  ``hadamard=True`` to cut its gates with the Hadamard test. ``joint`` turns on the
  joint wire measurement and the joint Hadamard cut (``None``, the default, keeps each
  cut's own default). ``shots=None`` gives the exact result (statevector, qibo or
  Qiskit), mid-circuit measurements included. ``gate_wire_cutting`` is still there and
  forwards to it when given ``hadamard``, ``optimized`` or ``joint``.

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
    cut = qd.find_cut(circuit, implementation="ibm-ckt")  # IBM's cut finder instead

**Hardware-aware cut selection (SparseCut).** When you pass the target QPU's
``coupling_map``, ``find_cut`` switches to a hardware-aware mode that ranks the
two-qubit gates by the physical distance between their qubits and cuts the
**longest-routing** ones — removing exactly the gates that would otherwise need the
most SWAPs (and so cause the most crosstalk / error) on a connectivity-limited
device. This follows the *SparseCut* idea of Guo et al.,
`arXiv:2511.05492 <https://arxiv.org/abs/2511.05492>`_.

.. code-block:: python

    line = [(0, 1), (1, 2), (2, 3), (3, 4)]           # the QPU's connectivity
    cut = qd.find_cut(circuit, coupling_map=line, max_cuts=2)
    cut = qd.find_cut(circuit, coupling_map=line, layout={0: 2, 1: 0, 2: 4})

The result is a gate cut you feed to ``gate_cutting`` as usual; pass a ``layout``
(logical → physical qubit map) if your circuit is not placed identity-first.

**Mixed cut discovery.** With ``mixed=True``, ``find_cut`` partitions the circuit once
and classifies each severed connection as a gate cut or a wire cut, returning a
combined ``{"gate": [...], "wire": [...]}`` you feed straight to
:func:`~Qdislib.circuit_cutting`. The partition is **cost-aware** — it weighs a gate
cut (``6`` terms) against a wire cut (``8``) and prefers whichever is cheaper per
severed connection — so a circuit with both a gate bottleneck and a wire bottleneck is
cut with the cheapest mix rather than forced into one type.

.. code-block:: python

    cut = qd.find_cut(circuit, mixed=True)            # {"gate": [...], "wire": [...]}
    val = qd.circuit_cutting(circuit, cut, "ZZZZ")

.. note::

   With the *optimized* wire cut the per-cut cost is closer to ``4`` than ``8``, which
   can tip the balance toward wire cuts; the cost model used for the default mixed
   discovery assumes the standard (``8``-term) wire cut, i.e. ``optimized=False``.

.. tip::

   Each cut multiplies the number of subcircuit evaluations — a **gate cut by 6**,
   a **wire cut by 8** — so the cost grows *exponentially* in the number of cuts.
   Preview it before running with ``qd.cut_cost``:

   .. code-block:: python

       qd.cut_cost(cut)   # {'num_cuts': 2, 'kind': 'gate', 'base': 6, 'terms': 36}

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

    # 1. Cut the circuit into independently executable subcircuits.
    #    software is one of "qibo", "qiskit" or "cudaq".
    subcircuits = qd.gate_cutting_subcircuits(circuit, gate_cut, software="qibo")
    # (or qd.wire_cutting_subcircuits(circuit, wire_cut, software="qibo"))

    # 2. Evaluate each subcircuit — this loop is yours to control.
    values = [qd.expectation_value(s, shots=8192) for s in subcircuits]

    # 3. Recombine into the original circuit's expectation value.
    result = subcircuits.reconstruct(values)

The returned ``subcircuits`` object behaves like a list of subcomponents (so it
can be iterated, indexed and executed independently using Qibo, Qiskit, CUDA-Q,
PennyLane, or submitted to QPUs) and additionally carries the cut metadata needed by
:meth:`~Qdislib.utils.backends.CutSubcircuits.reconstruct`.

If you prefer a single call, ``qd.gate_cutting(...)`` / ``qd.wire_cutting(...)``
run the cut, execution and reconstruction together. This modular structure gives
users fine-grained control over hybrid quantum-classical workflows.


Observables
-----------

The ``observables`` argument accepts a per-qubit Pauli string over
``X``, ``Y``, ``Z`` and ``I`` (identity), e.g. ``"ZIZ"``. ``X`` and ``Y`` terms
are measured by rotating into the Z basis automatically, on any backend.

You can also pass a **weighted sum of Pauli operators** (a Hamiltonian),
evaluated as ``<sum_i c_i P_i> = sum_i c_i <P_i>`` — either as a list of
``(coefficient, pauli_string)`` pairs or as a ``{pauli_string: coefficient}``
mapping:

.. code-block:: python

    # <0.5 ZZ> + <1.0 XI> - <0.3 IZ>
    hamiltonian = [(0.5, "ZZ"), (1.0, "XI"), (-0.3, "IZ")]
    result = qd.gate_cutting(circuit, cut, observables=hamiltonian)

    # equivalently, as a mapping
    result = qd.gate_cutting(circuit, cut, observables={"ZZ": 0.5, "XI": 1.0, "IZ": -0.3})

The same forms work with ``qd.wire_cutting``, ``qd.gate_cutting_hadamard`` (and
its ``_subcircuits`` variant) and ``qd.analytical_solution``.


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
    # or  qd.gate_cutting_hadamard_sampling(circuit, cut, shots=200_000, joint=True)

    # three-step: sample each subcircuit yourself, then reconstruct
    subs = qd.gate_cutting_hadamard_sampling_subcircuits(circuit, cut)
    dist = subs.reconstruct([qd.sample_subcircuit(s, shots=8192) for s in subs])

For the Hadamard cut, every fragment is read in the Z basis anyway, so the
distribution costs the same circuits as one expectation value, and with
``shots=None`` (``software="qibo"``) it is exact.

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

The ``software`` argument selects the simulator the subcircuits run on —
``"qibo"``, ``"qiskit"`` or ``"cudaq"`` — and is accepted by every cutting call.
Under the hood a pluggable **backend layer** turns your ``software`` / ``gpu`` /
``qpu`` arguments into a concrete backend:

.. code-block:: python

    qd.available_backends()                 # registered backends
    qd.build_execution_backend("qiskit")    # the single place that chooses one
    qd.QdBackEnd(qpu_dict={...})            # load-balance across a pool of devices

You can register your own backend (an in-house simulator or a cloud QPU) by
subclassing :class:`~Qdislib.core.backend.base_backend.BaseBackend` and calling
``qd.register_backend(...)``; the cutting code then selects it via
``software="MY_BACKEND"`` with no other change. See notebook *08 · Backends under
the hood*.


Resource Allocation
-------------------

Qdislib is designed to run efficiently across CPUs, GPUs, and quantum devices depending on the user’s configuration and available resources. By default, simulation is performed on CPU backends.

.. _gpu-support-label:

Using GPUs
----------

Qdislib offloads eligible subcircuit simulations to Qiskit Aer's GPU simulator
(``device="GPU"``, cuStateVec-enabled), for both ``gate_cutting`` and
``wire_cutting``. Pass ``gpu=True`` and a ``gpu_min_qubits`` threshold — only
subcircuits at least that large go to the GPU, the rest stay on CPU:

.. code-block:: python

    reconstruction = qd.gate_cutting(circuit, cut, gpu=True, gpu_min_qubits=18)

GPU usage is transparent: backend selection and memory transfers are handled
internally. The GPU path supports the same cases as the CPU path, including an
identity observable on a cut qubit.

Running on real QPUs
--------------------

Set ``qpu=True`` and a ``qpu_dict`` mapping a device name to its qubit count to
route each subcircuit to real hardware instead of a simulator:

.. code-block:: python

    # BSC (MareNostrum Ona) via the shared-filesystem bridge:
    value = qd.gate_cutting(circuit, cut, qpu=True, qpu_dict={"MN_Ona": 5})

    # IBM Quantum via qiskit-ibm-runtime:
    value = qd.gate_cutting(circuit, cut, qpu=True, qpu_dict={"IBM_Quantum": 127})

Qdislib routes each subcircuit to a large-enough configured QPU. A real device
returns only aggregated counts, so an identity observable on a cut qubit that
keeps evolving *after* the cut is not supported on QPU (it is on CPU/GPU). The
three-step workflow lets you submit each subcircuit to hardware yourself.
