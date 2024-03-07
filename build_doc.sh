#!/usr/bin/env bash

sphinx-apidoc -o docs/source src/Qdislib/ -f
cd docs
make html
cd ..
