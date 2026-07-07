#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#
# Run the Hadamard-test gate cut under `runcompss -g` and export the generated task
# dependency graph. Intended to run inside the COMPSs container (Dockerfile.compss),
# with the repo mounted at /Qdislib:
#
#     docker run --rm -v "$(pwd)":/Qdislib qdislib-compss \
#         bash /Qdislib/scripts/run_hadamard_graph.sh
#
# The rendered graph + the raw .dot land in /Qdislib/examples/ (i.e. on the host).

set -u
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}" || exit 1
OUT_DIR="${REPO_ROOT}/examples"

# COMPSs NIO workers connect over ssh to localhost; start sshd if it is not up.
if command -v sshd >/dev/null 2>&1 && ! pgrep -x sshd >/dev/null 2>&1; then
    /usr/sbin/sshd 2>/dev/null || service ssh start 2>/dev/null || true
    sleep 1
fi

APP_LOG="$(mktemp)"
echo "==> Hadamard gate cut under runcompss -g"
runcompss -g \
    --pythonpath="${REPO_ROOT}" \
    --python_interpreter=python3 \
    "${REPO_ROOT}/examples/hadamard_cut_compss.py" 2>&1 | tee "${APP_LOG}"

grep -q "RESULT: PASS" "${APP_LOG}" || echo "WARNING: app did not report PASS"

# Locate the freshest complete_graph.dot COMPSs just wrote.
DOT="$(ls -t "${HOME}"/.COMPSs/*/monitor/complete_graph.dot 2>/dev/null | head -1)"
if [ -z "${DOT}" ]; then
    echo "ERROR: no complete_graph.dot found under ~/.COMPSs/*/monitor/"
    exit 3
fi
echo "==> Found graph: ${DOT}"
cp "${DOT}" "${OUT_DIR}/hadamard_compss_graph.dot"

# Render to PNG (graphviz) and/or PDF (compss_gengraph), whichever is available.
if command -v dot >/dev/null 2>&1; then
    dot -Tpng "${DOT}" -o "${OUT_DIR}/hadamard_compss_graph.png" \
        && echo "==> wrote examples/hadamard_compss_graph.png"
fi
if command -v compss_gengraph >/dev/null 2>&1; then
    compss_gengraph "${DOT}" >/dev/null 2>&1 || true
    PDF="${DOT%.dot}.pdf"
    [ -f "${PDF}" ] && cp "${PDF}" "${OUT_DIR}/hadamard_compss_graph.pdf" \
        && echo "==> wrote examples/hadamard_compss_graph.pdf"
fi
echo "==> exported graph artifacts to examples/"
