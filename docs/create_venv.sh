#!/usr/bin/env bash

# This script is aimed at creating a virtual environment for building the documentation in local machines,
# with the objective of not messing around with the system or user packages.
# For Read the docs it will use directly the requirements.txt file.

# 1st step: pick a Python interpreter. requirements.txt pins Sphinx>=8, which needs
# Python >=3.10, but the system `python3` is often older (e.g. 3.9). Select the
# newest available >=3.10 interpreter so the venv can install the docs toolchain.
PYTHON=""
for candidate in python3.13 python3.12 python3.11 python3.10 python3; do
    command -v "$candidate" >/dev/null 2>&1 || continue
    if "$candidate" -c 'import sys; sys.exit(0 if sys.version_info[:2] >= (3, 10) else 1)'; then
        PYTHON="$candidate"
        break
    fi
done
if [ -z "$PYTHON" ]; then
    echo "ERROR: no Python >=3.10 found on PATH (Sphinx>=8 in requirements.txt needs it)." >&2
    echo "       Install Python 3.10+ (e.g. via pyenv/homebrew) and re-run." >&2
    exit 1
fi
echo "Building the docs venv with $PYTHON ($("$PYTHON" --version 2>&1))."

# 2nd step: create virtual environment
"$PYTHON" -m venv documentation-builder

# 3rd step: activate the virtual environment
source documentation-builder/bin/activate

# 4th step: install the necessary requirements
./install_dependencies.sh

# 5th step: deactivate environment
deactivate
