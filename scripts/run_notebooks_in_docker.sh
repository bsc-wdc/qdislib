#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#

# Execute the learning-path notebooks *with their outputs* inside the
# qdislib-compss container (PyCOMPSs + CUDA-Q + Jupyter; see Dockerfile.compss),
# writing the executed .ipynb back in place. Unlike run_notebook_compss.sh (which
# converts to a plain script), this keeps the Jupyter cell outputs -- so the
# notebooks can be committed/rendered with results embedded.
#
# The PyCOMPSs notebooks (07-09) drive the runtime themselves
# (ipycompss.start()/stop()) and need COMPSs' local ssh workers; the light ones
# (01-05) run without it. nbconvert's kernel can be flaky under QEMU emulation
# (amd64 image on an arm64 host) -- bump NB_TIMEOUT and re-run stragglers.
#
#   bash scripts/run_notebooks_in_docker.sh                     # all notebooks
#   bash scripts/run_notebooks_in_docker.sh notebooks/01_quickstart.ipynb
#   NB_TIMEOUT=7200 bash scripts/run_notebooks_in_docker.sh notebooks/08_vqe_h2.ipynb

set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${QDISLIB_COMPSS_IMAGE:-qdislib-compss:latest}"
TIMEOUT="${NB_TIMEOUT:-3600}"          # per-cell execution timeout (s)
STARTUP_TIMEOUT="${NB_STARTUP_TIMEOUT:-300}"   # kernel-start timeout (s); QEMU is slow

NBS=("$@")
if [ ${#NBS[@]} -eq 0 ]; then
    NBS=(
        notebooks/01_quickstart.ipynb
        notebooks/02_cutting_workflow.ipynb
        notebooks/03_backends.ipynb
        notebooks/04_observables.ipynb
        notebooks/05_caching.ipynb
        notebooks/07_image_reconstruction.ipynb
        notebooks/08_vqe_h2.ipynb
        notebooks/09_scaling_up.ipynb
    )
fi

echo "Image: ${IMAGE} | timeout: ${TIMEOUT}s | notebooks: ${NBS[*]}"

# The heredoc runs inside the container: set up the paths, start COMPSs' ssh
# workers (harmless for the non-PyCOMPSs notebooks), then execute each notebook
# in place with nbconvert, keeping outputs.
docker run --rm -v "${REPO}:/work" -w /work "${IMAGE}" bash -lc '
    set -u
    export PYTHONPATH="/work:/opt/COMPSs/Bindings/python/3/:${PYTHONPATH:-}"
    export MPLBACKEND=Agg
    if command -v sshd >/dev/null 2>&1 && ! pgrep -x sshd >/dev/null 2>&1; then
        (/usr/sbin/sshd || service ssh start) >/dev/null 2>&1 || true
        sleep 1
    fi
    rc=0
    for nb in "$@"; do
        echo "=== executing ${nb} ==="
        if jupyter nbconvert --to notebook --execute --inplace "${nb}" \
                --ExecutePreprocessor.timeout='"${TIMEOUT}"' \
                --ExecutePreprocessor.startup_timeout='"${STARTUP_TIMEOUT}"' \
                --ExecutePreprocessor.kernel_name=python3; then
            echo "OK   ${nb}"
        else
            echo "FAIL ${nb}"
            rc=1
        fi
    done
    exit ${rc}
' _ "${NBS[@]}"
