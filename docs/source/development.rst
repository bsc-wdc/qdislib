Developer's guide
=================

Contributing
------------
See the `Contributing Guide <https://github.com/bsc-wdc/qdislib/blob/main/CONTRIBUTING.md>`_.


Tests
-----

.. code-block:: bash

    pytest tests                         # serial, with a strict fake PyCOMPSs
    bash scripts/run_compss_tests.sh     # under runcompss, inside the COMPSs container
    scripts/mn5/check.sh engine <account> <qos>   # on MareNostrum 5, see the HPC guide

See :doc:`guide/hpc` for the last two.

Documentation
-------------

The pages are in ``docs/source``, and ``docs/requirements.txt`` is what Read the Docs installs. nbsphinx also needs ``pandoc``, which Read the Docs has.

.. code-block:: bash

    pip install -r docs/requirements.txt
    sphinx-build -b html docs/source /tmp/qdislib-docs

``docs/check_code_blocks.py`` runs every Python block of the pages. Mark one it cannot run with the comment ``.. not run: <why>`` just before it.

The notebooks are the source of what they show, stored with their outputs; Read the Docs does not run them. ``notebooks/run_notebooks.py`` runs them in place and cleans what they store. Run it in an environment without CUDA-Q, and without PyZX except for ``07_caching``: importing CUDA-Q prints a warning with a local path, and importing PyZX in a notebook turns its figures into SVG, which the PDF build cannot use. ``09_distributed_and_qpu`` has its outputs from a real COMPSs run, by ``notebooks/run_under_compss.py``, whose docstring has the command.

The Colab and Binder links in the notebooks name a branch of the GitHub repository (``main``), which is a mirror updated at each release, so a new notebook's links work from then on. ``binder/`` configures Binder.

Drafting new releases
---------------------

1. Ensure that the tests pass and that the
   `latest version of the documentation <https://qdislib.readthedocs.io/en/latest/>`_
   builds.

2. Decide whether to issue a minor or a major release following this
   `guide <https://semver.org/>`_.

3. Create and switch to a new branch named ``release-X.Y``.

4. Update the release number in ``Qdislib/__init__.py`` (the ``__version__`` attribute); the
   build reads it dynamically through ``pyproject.toml``.

5. Update the required PyCOMPSs version in ``QUICKSTART.md`` if necessary.

6. Update ``CHANGELOG.md``, and the version in the notebooks' install cell if the major version changes.

7. Run ``scripts/mn5/check.sh`` on MareNostrum 5.

8. Merge the branch into the main branch, and tag it ``vX.Y.Z``. The GitHub repository is synchronised from there at the release.

9. Build the package and upload it to PyPI:

   .. code-block:: bash

      pip3 install --upgrade build twine
      ./build.sh
      python3 -m twine upload dist/qdislib-X.Y.Z*
