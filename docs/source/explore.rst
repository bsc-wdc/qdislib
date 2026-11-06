Explore Qdislib
===============

The whole library, on one page, with nothing to install. The first section computes everything in your browser. The others show stored results of real Qdislib runs, and say so. To run the real code, open :doc:`the notebooks <examples>` in Colab or Binder, or :doc:`try the poster's four methods <try>`.

.. raw:: html

   <link rel="stylesheet" href="_static/try.css">
   <link rel="stylesheet" href="_static/explore.css">
   <script src="_static/explore_engine.js"></script>
   <script src="_static/explore_data.js"></script>

Circuit playground
------------------

Build a circuit of up to 6 qubits and cut it. **This runs live, in your browser:** a statevector simulation of each fragment, with the actual decompositions of a gate cut, a Hadamard cut and a wire cut. ``docs/check_explore_engine.py`` checks these numbers against real Qdislib on random circuits. Hadamard cuts here use one ancilla per gate, Qdislib's ``joint=False``.

.. raw:: html

   <div id="qdx-play" class="qd-try"><p>This section needs JavaScript.</p></div>

Finding cuts
------------

Pick a circuit and a device size, and see what :func:`~Qdislib.find_cut` chooses. *Stored from real runs.*

.. raw:: html

   <div id="qdx-find" class="qd-try"><p>This section needs JavaScript.</p></div>

Planning a run
--------------

What a number of cuts costs for each kind of cut, before anything runs. *Stored from* :func:`~Qdislib.cutting_overhead` *and* :func:`~Qdislib.shots_for`.

.. raw:: html

   <div id="qdx-plan" class="qd-try"><p>This section needs JavaScript.</p></div>

Distributions
-------------

The ``_sampling`` functions return the output distribution: all of it, the qubits you name, or over more than 24 qubits a :class:`~Qdislib.FactoredDistribution`. *Stored from real runs.*

.. raw:: html

   <div id="qdx-dist" class="qd-try"><p>This section needs JavaScript.</p></div>

Error bars
----------

``error=True`` gives a standard error from one run. Is it right? Compare it with the spread of many runs. *Stored from real runs.*

.. raw:: html

   <div id="qdx-err" class="qd-try"><p>This section needs JavaScript.</p></div>

Noise
-----

An Aer noise model on the cut and the uncut circuit. *Stored from real runs.*

.. raw:: html

   <div id="qdx-noise" class="qd-try"><p>This section needs JavaScript.</p></div>

Any SDK
-------

The same circuit and cut in five SDKs. *Each snippet was run; its results are stored.*

.. raw:: html

   <div id="qdx-sdk" class="qd-try"><p>This section needs JavaScript.</p></div>

Packing on a QPU
----------------

:func:`qd.pack <Qdislib.pack>` runs several subcircuits side by side, on groups of a device's qubits. *Stored from a real call on a fake device.*

.. raw:: html

   <div id="qdx-pack" class="qd-try"><p>This section needs JavaScript.</p></div>
   <script src="_static/explore.js"></script>

.. only:: latex

   The interactive page is at https://qdislib.readthedocs.io/en/latest/explore.html.
