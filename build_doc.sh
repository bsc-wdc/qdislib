#!/usr/bin/env bash

sphinx-apidoc -o docs/source Qdislib/ -f
cd docs
make html
cd ..
