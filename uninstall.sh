#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
  pip3 uninstall qdislib
else
  python3 -m pip uninstall qdislib
fi
