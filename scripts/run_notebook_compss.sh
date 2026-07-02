#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#
# Execute a Qdislib learning-path notebook end-to-end under PyCOMPSs inside a
# COMPSs container (see Dockerfile.compss). The notebook drives the runtime
# itself (``ipycompss.start()`` / ``ipycompss.stop()``); this wrapper only makes
# sure sshd is up -- COMPSs' local workers connect over ssh to localhost, same as
# scripts/run_compss_tests.sh -- and runs the notebook with nbconvert.
#
#   bash scripts/run_notebook_compss.sh [notebook.ipynb]

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}" || exit 1
NB="${1:-notebooks/07_image_reconstruction.ipynb}"

if command -v sshd >/dev/null 2>&1 && ! pgrep -x sshd >/dev/null 2>&1; then
    /usr/sbin/sshd 2>/dev/null || service ssh start 2>/dev/null || true
    sleep 1
fi

export PYTHONPATH="${REPO_ROOT}:/opt/COMPSs/Bindings/python/3/:${PYTHONPATH}"
export MPLBACKEND=Agg   # headless: plt.show() becomes a no-op

# nbconvert --execute needs a live Jupyter kernel, which is flaky under QEMU
# emulation ("Kernel didn't respond"). Instead export the notebook to a plain
# script and run it under IPython -- ipycompss.start() needs the IPython context
# (get_ipython), and IPython provides it for a script too. This executes the exact
# same notebook code end-to-end (encode -> entangle -> cut -> reconstruct).
case "${NB}" in
    *.py) script="${NB}" ;;
    *)
        script="/tmp/$(basename "${NB}" .ipynb).py"
        jupyter nbconvert --to script "${NB}" --stdout >"${script}"
        ;;
esac

echo "==> Running ${script} under PyCOMPSs (ipython)"
ipython "${script}"
status=$?
[ "${status}" -eq 0 ] && echo "==> OK" || echo "==> FAILED (exit ${status})"
exit "${status}"
