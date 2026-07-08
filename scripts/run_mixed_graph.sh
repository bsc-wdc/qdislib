#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#
# Run mixed wire+gate cutting under `runcompss -g` and export the task graph.
# Intended for the COMPSs container with the repo mounted at /Qdislib:
#     docker run --rm -v "$(pwd)":/Qdislib qdislib-compss \
#         bash /Qdislib/scripts/run_mixed_graph.sh

set -u
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}" || exit 1
OUT_DIR="${REPO_ROOT}/examples"

if command -v sshd >/dev/null 2>&1 && ! pgrep -x sshd >/dev/null 2>&1; then
    /usr/sbin/sshd 2>/dev/null || service ssh start 2>/dev/null || true
    sleep 1
fi

APP_LOG="$(mktemp)"
echo "==> Mixed wire+gate cutting under runcompss -g"
runcompss -g \
    --pythonpath="${REPO_ROOT}" \
    --python_interpreter=python3 \
    "${REPO_ROOT}/examples/mixed_cutting_compss.py" 2>&1 | tee "${APP_LOG}"

grep -q "RESULT: PASS" "${APP_LOG}" || echo "WARNING: app did not report PASS"

DOT="$(ls -t "${HOME}"/.COMPSs/*/monitor/complete_graph.dot 2>/dev/null | head -1)"
if [ -z "${DOT}" ]; then
    echo "ERROR: no complete_graph.dot found under ~/.COMPSs/*/monitor/"
    exit 3
fi
echo "==> Found graph: ${DOT}"
cp "${DOT}" "${OUT_DIR}/mixed_compss_graph.dot"
if command -v compss_gengraph >/dev/null 2>&1; then
    compss_gengraph png "${DOT}" >/dev/null 2>&1 \
        && cp "${DOT%.dot}.png" "${OUT_DIR}/mixed_compss_graph.png" \
        && echo "==> wrote examples/mixed_compss_graph.png"
elif command -v dot >/dev/null 2>&1; then
    dot -Tpng "${DOT}" -o "${OUT_DIR}/mixed_compss_graph.png"
fi
echo "==> exported graph artifacts to examples/"
