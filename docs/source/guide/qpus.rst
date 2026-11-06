Running on real QPUs
====================

.. note::

   The code on this page is not run by the documentation: it needs credentials and a device. The notebook
   :doc:`../notebooks/13_packing` runs the packing part on a fake backend with no credentials.

Pass the device as ``backend``, by its registered name or as a configured backend object. Every subcircuit then runs on that device. To send only the subcircuits that fit and simulate the rest, choose the backend per subcircuit in the :doc:`three-step loop <three_step>`.

Device noise is amplified by the cutting overhead, and more so with every extra cut, so expect hardware results to drift further from the exact value as cuts are added.

IBM Quantum
-----------

The backend connects once, here, through ``qiskit-ibm-runtime``. Token, channel and instance go in ``extra=`` or in the ``IBM_QUANTUM_TOKEN``, ``IBM_QUANTUM_CHANNEL`` and ``IBM_QUANTUM_INSTANCE`` environment variables. Without ``backend_name`` it takes the least busy device.

.. not run: needs IBM credentials

.. code-block:: python

    ibm = qd.build_backend(
        "IBM_BACKEND",
        qd.BackendConfig(num_qubits=156, extra={"backend_name": "ibm_marrakesh"}),
    )
    value = qd.gate_cutting(circuit, cut, backend=ibm)

Each subcircuit is one plain job. ``mode="batch"`` or ``mode="session"`` in ``extra``, or ``IBM_QUANTUM_MODE``, opens a Batch or Session where your plan allows one. ``extra={"options": {...}}`` is passed to ``SamplerV2``, for twirling or dynamical decoupling. The joint final and cut-point counts of a gate cut's mid-circuit measurement are read per shot.

Under PyCOMPSs the backend reaches a task as its device's name and settings, not as a connection or a token, because task arguments are written to files. Each worker connects again, reading ``IBM_QUANTUM_TOKEN`` from its environment. The workers need a route to IBM: nodes without internet access, such as MareNostrum 5's, cannot submit these jobs.

IQM
---

``pip install "qdislib[iqm]"``. ``token``, ``tokens_path``, ``url`` and ``quantum_computer`` each fall back to their ``IQM_*`` environment variable. The changelog records a 4-qubit circuit cut into 12 subcircuits on IQM Resonance's Garnet, which reconstructed ``<ZZZZ> = 0.852`` against an exact 1.0, the gap being hardware noise.

.. not run: needs IQM credentials

.. code-block:: python

    iqm = qd.build_backend("IQM_BACKEND", qd.BackendConfig(quantum_computer="garnet"))
    subs = qd.gate_cutting_subcircuits(circuit, ["CZ_3"], "ZZZZ")
    values = [qd.expectation_value(s, shots=2000, backend=iqm) for s in subs]
    value = subs.reconstruct(values)

BSC QPUs
--------

At BSC the classical cut job and the quantum-control job run on separate clusters that cannot share a scheduler allocation, so they communicate through a shared directory, ``$QDISLIB_SUBCIRCUITS_DIR``, which you have to set: it has no default. The cut job writes each subcircuit there and waits, bounded by ``$QDISLIB_QPU_TIMEOUT``, for its result; a listener on the quantum cluster, ``examples/qpu_listener.py``, runs it on a QPU and writes the result back.

.. code-block:: bash

    # On the quantum cluster (pip install "qdislib[qpu]" there)
    export QDISLIB_SUBCIRCUITS_DIR=/gpfs/projects/.../subcircuits
    export RUNCARD_RED=/path/to/red_runcard.yml
    python examples/qpu_listener.py --qpu red --nshots 1000

.. not run: needs the BSC QPU

.. code-block:: python

    # On the classical cluster, with the same QDISLIB_SUBCIRCUITS_DIR
    value = qd.gate_cutting(circuit, cut, backend="BSC_QPU_BACKEND")

A device returns only aggregated counts, so a cut point whose qubit keeps evolving is deferred to the end, on a fresh qubit, which is exact and costs one qubit more.

Many subcircuits at once: packing
---------------------------------

A device usually has many more qubits than one subcircuit uses. ``qd.pack(circuits, backend)`` compiles each circuit onto its own group of well calibrated, connected qubits, with groups at least ``buffer`` idle qubits apart, chosen from the backend's calibrated errors, and puts several of them in one device-wide circuit. ``unpack`` splits the results back, keyed as if each circuit had run alone:

.. not run: needs IBM credentials

.. code-block:: python

    from qiskit_ibm_runtime import SamplerV2

    circuits = []
    for c in subcircuits.circuits("qiskit", defer_measurements=True):
        c = c.copy()
        c.measure_all()
        circuits.append(c)

    packed = qd.pack(circuits, backend)   # a Qiskit BackendV2
    results = [SamplerV2(mode=backend).run([c], shots=8192).result()[0]
               for c in packed.circuits]  # one plain job per packed circuit
    value = subcircuits.reconstruct_counts(packed.unpack(results), bit_order="q0_last")

``packed.draw()`` shows where they went: the device's qubits and couplers, each circuit's group in its own colour, one job at a time with ``job=``. ``packed.circuits`` are already compiled for the device: run them as they are. ``unpack`` also takes counts dicts. ``packed.slots`` says where each circuit went and ``packed.two_qubit_gates`` what it cost once compiled. How the groups are chosen may still change between versions, so the same call can place circuits differently; ``regions=`` fixes them. So far packing has been checked on one device (ibm_kobe) with shallow 3-5-qubit circuits: packed and alone gave the same fidelity there, but deeper circuits may show more crosstalk, so compare a few circuits alone first, with ``qd.pack([c], backend, regions=[group])``.

The IBM backend does all of this itself. With ``pack=True`` each group of ``circuits_per_task`` subcircuits runs side by side, in plain jobs. A group holds one fragment's circuits.

.. not run: needs IBM credentials

.. code-block:: python

    ibm = qd.build_backend("IBM_BACKEND", qd.BackendConfig(
        extra={"backend_name": "ibm_kobe", "pack": True}))
    dist = qd.circuit_cutting_sampling(qc, cuts, backend=ibm, shots=4096, circuits_per_task=15)
