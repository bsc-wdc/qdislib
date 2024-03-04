#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
  pip3 uninstall Qdislib
else
  python3 -m pip uninstall Qdislib
fi
