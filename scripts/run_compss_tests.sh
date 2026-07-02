#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#
# Exercise Qdislib's PyCOMPSs (@task / compss_wait_on) parallel path by running
# it under `runcompss`. The local unit tests use the serial fallback, so this is
# the only thing that validates the distributed execution path (circuit
# serialisation across tasks, reconstruction synchronisation, etc.).
#
# Intended to run inside a COMPSs container (see Dockerfile.compss):
#     bash scripts/run_compss_tests.sh

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO_ROOT}" || exit 1

if ! command -v runcompss >/dev/null 2>&1; then
    echo "ERROR: 'runcompss' not found. Run this inside a COMPSs environment."
    exit 2
fi

# COMPSs starts its (local) NIO workers over SSH to localhost, so sshd must be
# running -- otherwise the worker fails with "ssh: connect to host localhost
# port 22: Connection refused". Start it if it is not already up.
if command -v sshd >/dev/null 2>&1 && ! pgrep -x sshd >/dev/null 2>&1; then
    /usr/sbin/sshd 2>/dev/null || service ssh start 2>/dev/null || true
    sleep 1
fi

status=0
APP_LOG="$(mktemp)"
SAMPLE_LOG="$(mktemp)"
TEST_LOG="$(mktemp)"

echo "==> [1/3] Example application (expectation) under runcompss"
# -g generates the task dependency graph under $HOME/.COMPSs/<app>/monitor/
runcompss -g \
    --pythonpath="${REPO_ROOT}" \
    --python_interpreter=python3 \
    "${REPO_ROOT}/examples/compss_app.py" 2>&1 | tee "${APP_LOG}"
grep -q "RESULT: PASS" "${APP_LOG}" || {
    echo "Example app did not report PASS"
    status=1
}

echo "==> [2/3] Sampling reconstruction under runcompss (qiskit, cudaq)"
# Exercises the gate_cutting_sampling fan-out (_sample_gate_subcircuit @task) and
# checks the reconstruction against the exact distribution -- the distributed
# counterpart of tests/test_gate_cutting_sampling.py.
QDISLIB_SOFTWARE="qiskit,cudaq" QDISLIB_SHOTS=8192 runcompss -g \
    --pythonpath="${REPO_ROOT}" \
    --python_interpreter=python3 \
    "${REPO_ROOT}/examples/verify_sampling_pycompss.py" 2>&1 | tee "${SAMPLE_LOG}"
grep -q "RESULT: PASS" "${SAMPLE_LOG}" || {
    echo "Sampling verification did not report PASS"
    status=1
}

# The full suite under runcompss is the slow part; set COMPSS_SKIP_SUITE=1 to run
# only the fast example + sampling checks (e.g. when iterating on a backend fix).
if [ "${COMPSS_SKIP_SUITE:-0}" = "1" ]; then
    echo "==> [3/3] Full unittest suite SKIPPED (COMPSS_SKIP_SUITE=1)"
else
    echo "==> [3/3] Full unittest suite under runcompss"
    runcompss \
        --pythonpath="${REPO_ROOT}" \
        --python_interpreter=python3 \
        "${REPO_ROOT}/tests/main.py" 2>&1 | tee "${TEST_LOG}"
    # PyCOMPSs exits 0 even when unittests fail (the COMPSs run itself succeeded),
    # so inspect the captured unittest summary line instead of the exit code.
    grep -qE "^OK" "${TEST_LOG}" || {
        echo "Unittest suite reported failures under PyCOMPSs"
        status=1
    }
fi

if [ "${status}" -eq 0 ]; then
    echo "COMPSs runs OK"
else
    echo "COMPSs runs reported failures"
fi
exit "${status}"
