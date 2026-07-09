#!/usr/bin/env bash

# 1st step: activate the virtual environment
if [ ! -f "documentation-builder/bin/activate" ]; then
    echo "Virtual environment not found. Creating..."
    ./create_venv.sh
    echo "Virtual environment created: OK"
fi
source documentation-builder/bin/activate

make clean
# No -f: apidoc only creates pages for NEW modules and never overwrites the
# curated source/*.rst. Module prose lives in the module docstrings (rendered by
# `.. automodule::`), and the pages carry deliberate curation (e.g. the backend
# page lists only the public backends, and Qdislib.rst drops `:members:` to avoid
# double-documenting the re-exported public API). `-f` would clobber all of that.
sphinx-apidoc -o source ../Qdislib/
make html
make latexpdf  # Fails with latest versions, waiting to be patched.

# End step: deactivate environment
deactivate
