#!/usr/bin/env bash

# 1st step: activate the virtual environment
if [ ! -f "documentation-builder/bin/activate" ]; then
    echo "Virtual environment not found. Creating..."
    ./create_venv.sh
    echo "Virtual environment created: OK"
fi
source documentation-builder/bin/activate

make clean
sphinx-apidoc -o source ../Qdislib/ -f
make html
# make latexpdf  # Fails with latest versions, waiting to be patched.

# End step: deactivate environment
deactivate