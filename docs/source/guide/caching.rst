Caching
=======

The notebook :doc:`../notebooks/07_caching` runs everything on this page.

A cut runs many fragment circuits, and much of that work comes back: the terms of a Hamiltonian share every fragment
their basis changes do not touch, a variational loop changes only some fragments per step, and a run is often repeated.
The **cache** keeps what each fragment circuit returned, its joint counts or its exact correlators, so a circuit already
run is read back instead of run again, and one entry serves every observable later read from that circuit.

Turning it on
-------------

Use ``cache=`` on every cutting function for one call, or :func:`~Qdislib.use_cache` for a session:

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd

    circuit = QuantumCircuit(5)
    circuit.ry(0.5, range(5))
    for q in range(4):
        circuit.cz(q, q + 1)
    cut = ["CZ_7"]
    hamiltonian = [(0.5, "ZZZZZ"), (1.0, "XIIII")]

    value = qd.circuit_cutting(circuit, cut, hamiltonian, cache=True)

    qd.use_cache(True)                             # every call from now on
    qd.use_cache(None)                             # off again

    with qd.use_cache("memory"):
        value = qd.gate_cutting(circuit, cut)      # only inside the block
        off = qd.gate_cutting(circuit, cut, cache=False)   # not for this call

``cache=`` and ``use_cache`` take:

- ``True``: in memory; under PyCOMPSs, LMDB in the working directory or ``$QDISLIB_CACHE_DIR``, so every worker shares it;
- ``"memory"``;
- ``"lmdb:<path>"``: persistent on disk (``pip install "qdislib[cache]"``);
- ``"redis"`` or ``"redis://host:port/db"``;
- a store object from :func:`~Qdislib.default_circuit_cache`. ``store.stats()`` gives hits, misses and size.

Bound the in-memory store for long runs: ``max_entries`` evicts the least recently used entry, and ``ttl`` expires entries after
some seconds, which is useful after an SDK upgrade.

.. code-block:: python

    store = qd.default_circuit_cache("memory", max_entries=5000, ttl=3600)
    value = qd.gate_cutting(circuit, cut, cache=store)
    print(store.stats())

What counts as the same circuit
-------------------------------

``cache_key=`` chooses:

- ``"exact"``, the default: the same fragment circuit (gates, qubits and the terms it runs) on the same backend with the same
  shots. It is computed on the master without building the circuit, so a hit never becomes a task, and it never matches a
  different circuit.
- ``"semantic"``: also equivalent circuits, by a ZX-reduced graph hash computed in the task (``pip install "qdislib[cache]"`` for
  PyZX). ZX reduction is not a canonical form, so an equivalent circuit can still be missed, and is then run. The hash also looks only a few steps around each node, so two different circuits can, rarely, share a key, and the second then gets the first one's result, with no error. Qdislib warns about this once. Use ``"exact"`` where a wrong value would matter.

.. warning::

   ``"semantic"`` can, rarely, give two different circuits the same key, and the second then
   gets the first one's result, with no error. Its hash only looks a few steps around each ZX
   node, so layered circuits that repeat a block are the likeliest. Qdislib warns once when it
   is used. Use ``"exact"`` where a wrong value would matter.

A Redis store takes its password from ``$QDISLIB_REDIS_PASSWORD``, set where the tasks run, and
refuses one in a ``redis://`` URL: PyCOMPSs writes task arguments to files.

A shot result is reused by later runs: the key holds the backend and the shots, not the seed, so a repeated run returns the same
samples and costs no device time. ``cache_refresh=True`` runs every circuit again and replaces what is stored, for fresh samples,
for instance to repeat a run for error bars. An exact result is the same on every simulator and is shared between them. Fragments
wider than 18 qubits are not stored exactly, as their correlator vector would be too large.

The three-step workflow
-----------------------

``qd.expectation_value(s, cache=...)`` caches per subcircuit, as ``subs.reconstruct`` needs nothing extra:

.. code-block:: python

    subs = qd.gate_cutting_subcircuits(circuit, cut, "ZZZZZ")
    store = qd.default_circuit_cache("memory")
    values = [qd.expectation_value(s, shots=2000, seed=1, cache=store) for s in subs]
    print(subs.reconstruct(values), store.stats())

Under PyCOMPSs, on the paths that take ``cache_mode`` (``expectation_value`` on a native
circuit or a list's subcircuit, and the deprecated cut options), a task that finds a
circuit reserved by another waits for its result. ``$QDISLIB_PENDING_TIMEOUT`` on the
workers is how long, in seconds (default 600, ten minutes): past it the run raises and
the reservation is removed, so one left by a killed run costs at most that wait. Raise
it for fragments that take longer.
