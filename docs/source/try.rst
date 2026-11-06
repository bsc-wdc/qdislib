.. Landing page of the printed QR code: https://qdislib.readthedocs.io/en/latest/try.html
   This path is printed: do not move this file.

Try Qdislib
===========

Cut a circuit into pieces, run them, and get the uncut result back. It all runs in your browser.

.. raw:: html

   <link rel="stylesheet" href="_static/try.css">
   <link rel="stylesheet" href="_static/explore.css">
   <script src="_static/explore_engine.js"></script>
   <script src="_static/explore_sdk_data.js"></script>
   <div id="qdx-play" class="qd-try"><p>This playground needs JavaScript. The notebook below does the same with the real library.</p></div>

.. only:: latex

   The live playground is at https://qdislib.readthedocs.io/en/latest/try.html.

Any SDK
-------

The same circuit and cut in five SDKs. *Each snippet was run; its results are stored.*

.. raw:: html

   <div id="qdx-sdk" class="qd-try"><p>This section needs JavaScript.</p></div>
   <script src="_static/explore.js"></script>

Run the real code
-----------------

A notebook with the real library, in a free, temporary environment, with nothing to install on your device.

.. raw:: html

   <div class="qd-links">
     <a href="https://mybinder.org/v2/gh/bsc-wdc/qdislib/main?urlpath=lab%2Ftree%2Fnotebooks%2FQdislib_try_it.ipynb"><b>Binder, with real PyCOMPSs</b><span>No account. The subcircuits run as parallel tasks, and you see the task graph. Starting can take 1&ndash;3 minutes, or longer if the service is busy.</span></a>
     <a href="https://colab.research.google.com/github/bsc-wdc/qdislib/blob/main/notebooks/Qdislib_try_it.ipynb"><b>Colab</b><span>Needs a Google account. Its first cell installs Qdislib, and the subcircuits run one after another, without PyCOMPSs.</span></a>
     <a href="notebooks/Qdislib_try_it.html"><b>Just read it</b><span>The notebook with the results of a real run under COMPSs.</span></a>
   </div>

.. only:: latex

   :doc:`notebooks/Qdislib_try_it` is the notebook, with its results.

Install it yourself
-------------------

.. code-block:: bash

   pip install "qdislib[qiskit]"
   uv pip install "qdislib[qiskit]"     # the same, with uv, typically much faster

Qdislib runs serially as it is. For parallel runs add `PyCOMPSs <https://compss-doc.readthedocs.io/en/stable/Sections/01_Installation_Configuration.html>`__: see the :doc:`Quickstart <quickstart>`.

More
----

:doc:`Explore the whole library <explore>`, read :doc:`the notebooks <examples>`, or see the source on `GitHub <https://github.com/bsc-wdc/qdislib>`__.
