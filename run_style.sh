#!/bin/bash -e

# Runs black code style checks on the qdislib.
python3 -m black --line-length 79 ./src/
ev=$?
if [ "$ev" -ne 0 ]; then
  echo "[ERROR] black check failed with exit value: $ev"
  echo ""
  echo "Please, run:"
  echo "    black --line-length 79 $(pwd)/src"
  echo "Then, review changes and push them again."
  echo ""
  exit $ev
fi
