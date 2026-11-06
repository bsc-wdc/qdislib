Hardware and noise
==================

How subcircuits share a chip, and how cuts behave under noise.

.. raw:: html

   <link rel="stylesheet" href="_static/try.css">
   <link rel="stylesheet" href="_static/explore.css">
   <script src="_static/explore_engine.js"></script>
   <script src="_static/explore_data.js"></script>

Packing on a QPU
----------------

:func:`qd.pack <Qdislib.pack>` runs several subcircuits side by side, on groups of a device's qubits. *Stored from real calls on real device layouts: nothing ran on them. Each layout says where its calibration comes from.*

.. raw:: html

   <div id="qdx-pack" class="qd-try"><p>This section needs JavaScript.</p></div>

Noise
-----

How each kind of cut behaves as noise grows, next to the uncut circuit. *Stored from real Aer runs.*

.. raw:: html

   <div id="qdx-noise" class="qd-try"><p>This section needs JavaScript.</p></div>

.. raw:: html

   <script src="_static/explore.js"></script>

.. only:: latex

   The interactive page is at https://qdislib.readthedocs.io/en/latest/explore_hardware.html.
