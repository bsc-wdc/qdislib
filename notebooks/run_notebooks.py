#!/usr/bin/env python3
"""Execute the notebooks in place and clean what they store.

    python notebooks/run_notebooks.py                     # every notebook
    python notebooks/run_notebooks.py 01_quickstart.ipynb applications/02_scaling.ipynb

Run it from an environment with `pip install ".[qibo,qiskit]" nbclient ipykernel
ipywidgets`, without CUDA-Q (importing it prints a FutureWarning with a local path), and
without PyZX except for 07_caching.ipynb, which needs `.[cache]`: importing PyZX makes a
notebook's figures SVG, which the PDF build of the documentation cannot use.

Qdislib_try_it.ipynb and 09_distributed_and_qpu.ipynb are left out unless named: their
stored outputs come from a real COMPSs run, in the Binder image, as Binder runs it:

    docker build -f binder/Dockerfile -t qdislib-try .
    docker run --rm --user 1000 -v "$PWD/notebooks":/home/ubuntu/notebooks qdislib-try \
        python3 notebooks/run_notebooks.py Qdislib_try_it.ipynb 09_distributed_and_qpu.ipynb

Besides running each notebook this removes the live-widget outputs (cells tagged
`live-widget` show nothing on a static page) and the COMPSs runtime's banner, log path and
messages (cells tagged `compss-start` and `compss-stop`), keeps the notebook metadata to the
kernel, and fails if an output holds a local path.
"""

import argparse
import pathlib
import re
import sys
import time

import nbformat
from nbclient import NotebookClient
from nbformat.v4 import new_output

ROOT = pathlib.Path(__file__).resolve().parent
LOCAL_PATH = re.compile(r"/(Users|home|private|tmp|gpfs|root|opt)/[\w.-]+|site-packages")
RUN_UNDER_COMPSS = {"09_distributed_and_qpu.ipynb", "Qdislib_try_it.ipynb"}


def clean(notebook):
    """Strip metadata and live-widget outputs; return the local paths still in outputs."""
    notebook.metadata = {
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
        "language_info": {"name": "python"},
    }
    found = []
    under_compss = any(
        "compss-start" in cell.metadata.get("tags", []) for cell in notebook.cells
    )
    for cell in notebook.cells:
        metadata = {k: v for k, v in cell.metadata.items() if k in ("tags", "nbsphinx")}
        cell.metadata = metadata
        if cell.cell_type != "code":
            continue
        tags = metadata.get("tags", [])
        if "live-widget" in tags:
            cell.outputs = []
        if under_compss:
            for output in list(cell.outputs):
                if output.get("output_type") != "stream":
                    continue
                output.text = output.text.replace("INFO: Storage API successfully imported.\n", "")
                if output.name == "stderr" or "compss-stop" in tags or not output.text:
                    cell.outputs.remove(output)
            if "compss-start" in tags and "PyCOMPSs is not installed" not in str(cell.outputs):
                cell.outputs = [new_output("stream", name="stdout", text="PyCOMPSs runtime started.\n")]
        for output in cell.outputs:
            text = output.get("text", "")
            data = output.get("data", {})
            text += "".join(str(data.get(kind, "")) for kind in ("text/plain", "text/html"))
            text = re.sub(r"data:image/png;base64,[A-Za-z0-9+/=]+", "", text)
            text += "\n".join(output.get("traceback", []))
            found += LOCAL_PATH.findall(text) and [LOCAL_PATH.search(text).group(0)]
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("notebooks", nargs="*", help="paths relative to notebooks/")
    parser.add_argument("--timeout", type=int, default=600, help="per-cell seconds")
    args = parser.parse_args()
    paths = (
        [ROOT / name for name in args.notebooks]
        or [
            path
            for path in sorted(ROOT.glob("*.ipynb")) + sorted(ROOT.glob("applications/*.ipynb"))
            if path.name not in RUN_UNDER_COMPSS
        ]
    )
    failed = False
    for path in paths:
        notebook = nbformat.read(path, as_version=4)
        start = time.time()
        NotebookClient(
            notebook, timeout=args.timeout, kernel_name="python3",
            resources={"metadata": {"path": str(path.parent)}},
        ).execute()
        found = clean(notebook)
        nbformat.write(notebook, path)
        print(f"{time.time() - start:7.1f} s  {path.relative_to(ROOT)}", flush=True)
        if found:
            print(f"          local paths in outputs: {sorted(set(found))}", flush=True)
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
