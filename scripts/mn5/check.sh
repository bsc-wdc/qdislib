#!/usr/bin/env bash
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License").
#
# Before a release: Qdislib under real COMPSs on MareNostrum 5, on what CI cannot run.
#
#   scripts/mn5/check.sh engine <account> <qos>    # on a GPP login: engine options
#   scripts/mn5/check.sh gpu <account> <qos>       # on an ACC login: GPU tasks
#
# Copy the repository to MN5 first and run this from its root; it writes to
# $QDISLIB_CHECK_DIR (default: ./mn5-check) and prints the job. Each log ends with
# "RESULT: PASS" or "RESULT: FAIL ...": grep it in compss-<job>.out.
#
# COMPSS_MODULE chooses the COMPSs module. Look before trusting a version: an install
# can lack enqueue_compss (Runtime/scripts/user) or its configuration
# (Runtime/configuration/xml/resources: "Cannot parse resources.xml").

set -u
KIND="${1:?engine or gpu}"; ACCOUNT="${2:?account}"; QOS="${3:?qos}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${QDISLIB_CHECK_DIR:-$PWD/mn5-check}"
mkdir -p "$OUT"
export COMPSS_PYTHON_VERSION=3.12.1

# One at a time and not through a pipe: a piped `module load` runs in a subshell.
if [ "$KIND" = gpu ]; then
    LIB=/apps/ACC/COMPSs/LIBRARIES/QDISLIB/lib
    module load intel/2023.2.0 >/dev/null 2>&1   # ACC's python needs intel and mkl
    module load mkl >/dev/null 2>&1
    module load hdf5 >/dev/null 2>&1
    module load python/3.12.1 >/dev/null 2>&1
    module load "COMPSs/${COMPSS_MODULE:-TrunkCTClean}" >/dev/null 2>&1
    APP="$REPO/examples/verify_gpu_pycompss.py"
else
    LIB=/apps/GPP/COMPSs/LIBRARIES/QDISLIB/lib
    module load hdf5 python/3.12.1 >/dev/null 2>&1
    module load "COMPSs/${COMPSS_MODULE:-3.4.post2603}" >/dev/null 2>&1
    APP="$REPO/examples/verify_engine_pycompss.py"
fi
command -v enqueue_compss >/dev/null || { echo "no enqueue_compss in this COMPSs module"; exit 2; }
PY="$REPO:$REPO/examples:$LIB"
export PYTHONPATH="$PY:${PYTHONPATH:-}"
cat > "$OUT/env.sh" <<ENV
export QDISLIB_PROBE_FILE=$OUT/probe_\$SLURM_JOB_ID.jsonl
export OMP_NUM_THREADS=1
ENV
cd "$OUT" || exit 1
enqueue_compss --num_nodes=1 --exec_time=10 --qos="$QOS" --project_name="$ACCOUNT" \
    --lang=python --worker_working_dir="$OUT" --log_dir="$OUT" \
    --env_script="$OUT/env.sh" --pythonpath="$PY" --job_name="qd-$KIND" "$APP"
