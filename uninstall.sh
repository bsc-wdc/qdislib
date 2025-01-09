#!/usr/bin/env bash

# Required to clean in order to uninstall the package.
./clean.sh

if [ "$(uname)" == "Darwin" ]; then
  pip3 uninstall Qdislib
else
  python3 -m pip uninstall -y Qdislib
fi
