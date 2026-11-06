"""Run 09_distributed_and_qpu.ipynb under real COMPSs, storing what it shows.

In the container built from Dockerfile.compss, with the repository mounted:

    docker run --rm --platform linux/amd64 -v "$PWD":/code --entrypoint bash qdislib-compss -c \
        '/usr/sbin/sshd; sleep 1; cd /tmp; \
         PYTHONPATH=/code:/opt/COMPSs/Bindings/python/3/ ipython --no-banner \
         /code/notebooks/run_under_compss.py /code/notebooks/09_distributed_and_qpu.ipynb'

It runs inside IPython, not through a Jupyter kernel, because a kernel can fail to start
in an emulated container. It leaves out of the outputs the runtime's banner and log
path, and its messages on stderr. Pass a second path to write the result elsewhere.
"""

import base64
import os
import sys
import time

import nbformat
from IPython import get_ipython
from IPython.utils.capture import capture_output
from nbformat.v4 import new_output


def outputs_of(cell, shell):
    """Run one cell; return its outputs, and whether it failed."""
    with capture_output() as captured:
        result = shell.run_cell(cell.source, store_history=False)
    if "compss-start" in cell.metadata.get("tags", []):
        stdout = "PyCOMPSs runtime started.\n"
    elif "STOPPING PyCOMPSs" in captured.stdout:
        stdout = ""
    else:
        stdout = captured.stdout.replace("INFO: Storage API successfully imported.\n", "")
    outputs = [new_output("stream", name="stdout", text=stdout)] if stdout else []
    for display in captured.outputs:
        outputs.append(
            new_output(
                "display_data",
                data={
                    kind: value if isinstance(value, str) else base64.b64encode(value).decode()
                    for kind, value in display.data.items()
                },
            )
        )
    failed = result.error_in_exec or result.error_before_exec
    if failed:
        print(captured.stdout[-1500:], captured.stderr[-1500:], file=sys.__stdout__, flush=True)
    return outputs, failed


def main():
    source, target = sys.argv[1], sys.argv[-1]
    notebook = nbformat.read(source, as_version=4)
    shell = get_ipython()
    shell.run_line_magic("matplotlib", "inline")
    count = 0
    for cell in notebook.cells:
        if cell.cell_type != "code":
            continue
        count += 1
        start = time.time()
        cell.outputs, failed = outputs_of(cell, shell)
        cell.execution_count = count
        print(f"cell {count}: {time.time() - start:.1f} s", file=sys.__stdout__, flush=True)
        if failed:
            print("FAILED:", repr(failed), file=sys.__stdout__, flush=True)
            os._exit(1)
    nbformat.write(notebook, target)
    print("DONE", file=sys.__stdout__, flush=True)
    os._exit(0)  # the runtime's threads would keep the process alive


main()
