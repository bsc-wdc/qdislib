
## Quickstart guide

There are two ways in which you can get started with Qdislib. You can perform
a manual installation, or you can download our ready-to-use docker image.

### Manual installation

#### Dependencies

Qdislib currently requires:

* PyCOMPSs >= 3.3
* rustworkx >= 0.15.0
* numpy >= 1.17,<3
* scipy >= 1.5
* sympy >= 1.3
* dill >= 0.3
* python-dateutil >= 2.8.0
* stevedore >= 3.0.0
* typing-extensions
* symengine >= 0.11,<0.14
* igraph >= 0.11.8
* qiskit >= 1.1.0
* qiskit_aer >= 0.15.1
* qibo >= 0.1.12
* networkx >= 3.1
* matplotlib >= 3.4.3
* pymetis >= 2023.1.1

While in order to use GPUs, cupy and/or pytorch are also required.

#### Installation steps

1. Check which PyCOMPSs version to install.
    * Latest Qdislib release requires **PyCOMPSs 3.3** or greater (check [Qdislib releases](https://github.com/bsc-wdc/Qdislib/releases) for information about other releases).

2. Install PyCOMPSs following these [instructions](https://compss-doc.readthedocs.io/en/stable/Sections/01_Installation.html).

3. Install the latest Qdislib version with ``pip3 install Qdislib``.
   * **IMPORTANT:** Qdislib requires the ``pycompss`` Python module. However, this command will **NOT** install the module automatically. The module should be available after manually installing PyCOMPSs following the instructions in step 2. For more information on this, see [this issue](https://github.com/bsc-wdc/dislib/issues/190).

4. You can check that everything works fine by running one of our examples:

    * Download the latest source code [here](https://github.com/bsc-wdc/Qdislib/releases/latest).

    * Extract the contents of the tar package.

    ```bash
    tar xzvf Qdislib-X.Y.Z.tar.gz
    ```

    * Run an example application.

    ```bash
    runcompss Qdislib-X.Y.Z/examples/app.py
    ```
