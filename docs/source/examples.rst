Examples
========

A hands-on **learning path** — the notebooks build on each other, so read them
in order. Each is self-contained (its own imports and example circuit) and the
core ones run serially with just a Qibo install; later notebooks introduce
alternative backends (Qiskit / CUDA-Q), sampling, caching, and distribution.

.. toctree::
   :maxdepth: 2
   :caption: Learning path:

   notebooks/01_quickstart
   notebooks/02_cutting_workflow
   notebooks/03_finding_cuts
   notebooks/04_backends
   notebooks/05_observables
   notebooks/06_sampling
   notebooks/07_caching
   notebooks/08_backend_internals
   notebooks/09_distributed_and_qpu

End-to-end **applications** — larger, self-contained studies that put the whole
learning path to work (image reconstruction, a VQE chemistry run, and scaling to
bigger circuits). These are standalone, not part of the ordered path.

.. toctree::
   :maxdepth: 2
   :caption: Applications:

   notebooks/applications/01_image_reconstruction
   notebooks/applications/02_vqe_h2
   notebooks/applications/03_scaling
