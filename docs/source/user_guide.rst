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

- **Hadamard-test Gate Cutting** (``gate_cutting_hadamard``): An alternative gate cut (Harrow & Lowe, `arXiv:2403.01018 <https://arxiv.org/abs/2403.01018>`_) that replaces the mid-circuit measurement with **one local ancilla** per fragment running a Hadamard test. It cuts **any** two-qubit gate as a *single* cut at the **optimal** overhead (``gamma(CZ) = 3``, tending to ``1`` for weak interactions), which makes it a good fit for superconducting hardware where mid-circuit measurement is costly. Available as an all-in-one call, a three-step ``_subcircuits`` workflow, and a Monte-Carlo ``_sampling`` variant.

- **Wire Cutting** (``wire_cutting``): Cuts "wires" (qubit connections) between gates to disconnect logical qubits, via a measure-and-prepare decomposition (``8`` terms per cut, or ``3**n + 4**n`` with ``optimized=True``).

- **Mixed Gate + Wire Cutting** (``gate_wire_cutting``): Cuts some gates *and* some wires on the **same circuit** in one integrated workflow. The two quasiprobability decompositions tensor, and every fragment across the combined configuration space runs in parallel.

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
:func:`~Qdislib.gate_wire_cutting`. The partition is **cost-aware** — it weighs a gate
cut (``6`` terms) against a wire cut (``8``) and prefers whichever is cheaper per
severed connection — so a circuit with both a gate bottleneck and a wire bottleneck is
cut with the cheapest mix rather than forced into one type.

.. code-block:: python

    cut = qd.find_cut(circuit, mixed=True)            # {"gate": [...], "wire": [...]}
    val = qd.gate_wire_cutting(
        circuit, gates_cut=cut["gate"], wire_cut=cut["wire"], observable="ZZZZ"
    )

.. note::

   With the *optimized* wire cut the per-cut cost is closer to ``4`` than ``8``, which
   can tip the balance toward wire cuts; the cost model used for the default mixed
   discovery assumes the standard (``8``-term) wire cut that ``gate_wire_cutting`` runs.

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

The same forms work with ``qd.wire_cutting`` and ``qd.analytical_solution``.


Sampling the output distribution
--------------------------------

``gate_cutting`` / ``wire_cutting`` return a single **expectation value**. To
reconstruct the full **output distribution** ``{bitstring: probability}`` instead
— e.g. to read out an image or inspect a state — use the sampling variants, which
have the same all-in-one and three-step forms:

.. code-block:: python

    dist = qd.gate_cutting_sampling(circuit, cut, shots=200_000)
    # or  qd.wire_cutting_sampling(circuit, wire_cut, shots=200_000)

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
