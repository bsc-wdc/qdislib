Running with PyCOMPSs and on HPC
================================

The notebooks :doc:`../notebooks/09_distributed_and_qpu` and :doc:`../notebooks/Qdislib_try_it` were run under real COMPSs.

Qdislib runs serially out of the box. Every fragment circuit is a PyCOMPSs task, so when `PyCOMPSs <https://compss-doc.readthedocs.io/>`__ is installed, the same code runs them in parallel. It needs COMPSs 3.3 or later. Without it the tasks run one after another in the same process.

How to run
----------

.. code-block:: bash

    python my_app.py                                # serial
    runcompss --pythonpath=$(pwd) my_app.py         # on this machine, in parallel
    runcompss -g --pythonpath=$(pwd) my_app.py      # and write the task graph
    enqueue_compss --num_nodes=2 ... my_app.py      # a batch job on a cluster

In a notebook, start the runtime **before importing Qdislib**: importing it first loads PyCOMPSs outside the runtime, and ``start()`` then fails.

.. not run: needs PyCOMPSs

.. code-block:: python

    import pycompss.interactive as ipycompss
    ipycompss.start(graph=True, cpu_affinity="disabled")

    import Qdislib as qd
    ...
    ipycompss.stop(sync=True)

``cpu_affinity="disabled"``, or ``--cpu_affinity=disabled`` on ``runcompss``, skips COMPSs' CPU pinning, which fails on Apple Silicon, also in a container, leaving the run waiting.

What runs as a task
-------------------

Each fragment circuit is built and run in its own task, on the worker; the master ships each fragment once, plus a small spec per circuit. An expectation value then has one more task, which contracts the results. A distribution reduces each fragment's counts in its own task and contracts them in a last one, as do an expectation value from counts (``error=``, ``total_shots=``, ``result_from_counts``). In the three-step form each ``expectation_value`` and ``sample_subcircuit`` is a task, and ``reconstruct`` takes the futures.

Overlapping calls: ``sync=False``
---------------------------------

Every all-in-one call waits for its result by default. With ``sync=False`` it returns a future at once, so several calls overlap and all their tasks are in one pool. Use it whenever a program runs more than one cut.

.. code-block:: python

    from qiskit import QuantumCircuit
    import Qdislib as qd
    from Qdislib.utils.runtime import compss_wait_on

    circuit = QuantumCircuit(4)
    circuit.ry(0.5, range(4))
    for q in range(3):
        circuit.cz(q, q + 1)
    cut = ["CZ_6"]

    futures = [qd.gate_cutting(circuit, cut, shots=4096, seed=s, sync=False) for s in range(3)]
    print(compss_wait_on(futures))

``Qdislib.utils.runtime`` holds ``task`` and ``compss_wait_on`` with serial fallbacks, so this code also runs without PyCOMPSs.

Small circuits: ``circuits_per_task=``
--------------------------------------

Every cutting function can run its fragment circuits in groups instead of one task each: fewer, longer tasks when the circuits are small and many, where a task's overhead outweighs the simulation. The changelog gives about 70 ms for a 9-qubit subcircuit in Aer. It has the same seeds and the same results, but is not yet available with ``cache=`` or ``total_shots=``. A group returns one future, and stays inside one fragment. PyCOMPSs handles every element of a task's collection on the master, so grouping also lowers the number of serialised parameters: the changelog's case is 240 subcircuits per tile, 6 parameters instead of 240.

.. code-block:: python

    print(qd.circuit_cutting(circuit, cut, shots=4096, seed=1, circuits_per_task=4))

The cache under PyCOMPSs
------------------------

``cache=True`` is an LMDB store in the working directory, or ``$QDISLIB_CACHE_DIR``, which every worker shares. Point it at a shared filesystem on a cluster. See :doc:`caching`.

GPU tasks
---------

COMPSs assigns a GPU only to tasks that ask for one. A backend **explicitly** on a GPU, with ``use_gpu`` for Aer, ``target="nvidia"`` for CUDA-Q, or a custom one overriding ``uses_gpu``, runs its circuits in tasks that ask for one GPU and one CPU core. Nothing else is needed.

Run these on GPU nodes: on a cluster without GPUs, those tasks never start. A backend that only *may* use a GPU does not count, for that reason.

Backends in tasks
-----------------

A backend object goes to a task by pickle, as any argument. Your own backend's module must be importable on the workers, so put it on ``--pythonpath``. The IBM backend goes as its device's name and settings, and each worker connects again reading ``IBM_QUANTUM_TOKEN`` from its environment. See :doc:`qpus`.

MareNostrum 5
-------------

PyCOMPSs and Qdislib are installed on MareNostrum 5. Load the environment one module at a time and not through a pipe, since a piped ``module load`` runs in a subshell and is lost. On GPP:

.. code-block:: bash

    export COMPSS_PYTHON_VERSION=3.12.1
    module load hdf5 python/3.12.1
    module load COMPSs/3.4.post2603

    REPO=/path/to/your/qdislib
    export PYTHONPATH=$REPO:/apps/GPP/COMPSs/LIBRARIES/QDISLIB/lib:$PYTHONPATH
    enqueue_compss --num_nodes=1 --exec_time=10 --qos=<qos> --project_name=<account> \
        --lang=python --pythonpath="$PYTHONPATH" my_app.py

Compute nodes have no internet access: they cannot install packages or reach IBM. Pass an ``--env_script`` for what a worker needs in its environment, such as ``OMP_NUM_THREADS=1``.

Checking an installation
~~~~~~~~~~~~~~~~~~~~~~~~

``scripts/mn5/check.sh`` submits Qdislib's real COMPSs checks, for what CI cannot run. Copy the repository to MN5 and run it from its root on the matching login node:

.. code-block:: bash

    scripts/mn5/check.sh engine <account> <qos>    # on a GPP login: the engine's options
    scripts/mn5/check.sh gpu <account> <qos>       # on an ACC login: GPU tasks

It writes to ``$QDISLIB_CHECK_DIR``, by default ``./mn5-check``, and prints the job. Each log ends with ``RESULT: PASS`` or ``RESULT: FAIL ...``: grep for it in ``compss-<job>.out``. ``COMPSS_MODULE`` chooses the COMPSs module. Look before trusting a version: an install can lack ``enqueue_compss`` or its configuration.

In a container
--------------

``binder/Dockerfile`` is the image behind the notebooks' Binder links. It runs as a normal user, and serves JupyterLab:

.. code-block:: bash

    docker build -f binder/Dockerfile -t qdislib-try .
    docker run --rm -p 8888:8888 qdislib-try \
        jupyter lab --ip=0.0.0.0 --port=8888 --no-browser --IdentityProvider.token=qdislib

``Dockerfile.compss`` builds an image on the COMPSs tutorial image, for testing without a cluster:

.. code-block:: bash

    docker build -f Dockerfile.compss -t qdislib-compss .
    docker run --rm qdislib-compss       # scripts/run_compss_tests.sh

COMPSs starts its local workers over ssh, so a shell in the container has to start ``sshd`` first. ``runcompss`` exits 0 even when an inner test fails, so the script checks the output for ``RESULT: PASS`` and the unittest ``OK`` line.

Checks for the distributed path are in ``examples/``: see :doc:`../examples`.
