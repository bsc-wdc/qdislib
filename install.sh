#!/usr/bin/env bash

if [ "$(uname)" == "Darwin" ]; then
  pip3 install .
else
  python3 -m pip install .
fi
