
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

### Using docker

#### 1. Install Docker and docker-py

**Warning:** requires docker version >= 17.12.0-ce

1. Follow these instructions

    * [Docker for Mac](https://store.docker.com/editions/community/docker-ce-desktop-mac). Or, if you prefer to use [Homebrew](https://brew.sh/).
    * [Docker for Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce-1).
    * [Docker for Arch Linux](https://wiki.archlinux.org/index.php/Docker#Installation).

    Be aware that the docker package has been renamed from `docker` to `docker-ce` for some distributions.
    Make sure you install the new package.

2. Add user to docker group to run Qdislib as a non-root user: [example](https://docs.docker.com/install/linux/linux-postinstall/)

3. Check that docker is correctly installed.

    ```bash
    docker --version
    docker ps # this should be empty as no docker processes are yet running.
    ```

4. Install [docker-py](https://docker-py.readthedocs.io/en/stable/)

    ```bash
    pip3 install docker
    ```

#### 2. Install Qdislib

```bash
pip3 install Qdislib
```

This should add the Qdislib executable to your path.

#### 3. Start PyCOMPSs in your development directory

Initialize PyCOMPSs where your source code will be (you can re-init anytime).
This will allow docker to access your local code and run it inside the container.

**Note** that the first time PyCOMPSs needs to download the docker image from the registry, and it may take a while.

```bash
# Without a path it operates on the current working directory.
pycompss init

# You can also provide a path
pycompss init /home/user/replace/path/
```

#### 4. Running applications

**Note**: running the docker PyCOMPSs does not work with applications with GUI or with visual plots such as `examples/app.py`).

First clone Qdislib repo and checkout release branch vX.Y.Z (docker version and Qdislib code should preferably be the same to avoid inconsistencies):

```bash
git clone https://github.com/bsc-wdc/Qdislib.git
```

Init the PyCOMPSs environment in the root of the repo.
The source files path are resolved from the init directory which sometimes can be confusing.
As a rule of thumb, initialize the library in a current directory and check the paths are correct running the file with `python3 path_to/file.py` (in this case `python3 examples/app.py`).

```bash
cd Qdislib
pycompss init
pycompss exec examples/app.py
```

The log files of the execution can be found at `$HOME/.COMPSs`.

You can also init PyCOMPSs inside the examples folder.
This will mount the examples directory inside the container so you can execute it without adding the path:

```bash
cd Qdislib/examples
pycompss init
pycompss exec app.py
```

#### 5. Running Jupyter notebooks

Notebooks can be run using the `pycompss jupyter` command. Run the
following snippet from the root of the project:

```bash
pycompss init
pycompss jupyter ./notebooks
```

An alternative and more flexible way of starting jupyter is using the
`pycompss run` command in the following way:

```bash
pycompss run jupyter-notebook ./notebooks --ip=0.0.0.0  --allow-root
```

Access your notebook by CTRL-clicking or copy pasting into the browser the link shown on the CLI (e.g. `http://127.0.0.1:8888/?token=TOKEN_VALUE`).

If the notebook process is not properly closed, you might get the following warning when trying to start jupyter notebooks again:

`The port 8888 is already in use, trying another port.`

To fix it, just restart the PyCOMPSs container with `pycompss init`.


#### 6. Adding more nodes

**Note**: adding more nodes is still in beta phase. Please report issues, suggestions, or feature requests on [Github](https://github.com/bsc-wdc/Qdislib).

To add more computing nodes, you can either let docker create more workers for you or manually create and config a custom node.

For docker just issue the desired number of workers to be added. For example, to add 2 docker workers:

```bash
pycompss components add worker 2
```

You can check that both new computing nodes are up with:

```bash
pycompss components list
```

If you want to add a custom node it needs to be reachable through ssh without user.
Moreover, PyCOMPSs will try to copy the `working_dir` there, so it needs write permissions for the scp.

For example, to add the local machine as a worker node:

```bash
pycompss components add worker '127.0.0.1:6'
```

* '127.0.0.1': is the IP used for ssh (can also be a hostname like 'localhost' as long as it can be resolved).
* '6': desired number of available computing units for the new node.

**Please be aware** that `pycompss components` will not list your custom nodes because they are not docker processes and thus it can't be verified if they are up and running.
