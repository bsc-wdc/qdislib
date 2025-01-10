#!/usr/bin/env bash

# This script is aimed at creating a virtual environment for building the documentation in local machines,
# with the objective of not messing around with the system or user packages.
# For Read the docs it will use directly the requirements.txt file.

# 1st step: create virtual environment
python3 -m venv documentation-builder

# 2nd step: activate the virtual environment
source documentation-builder/bin/activate

# 3rd step: install the necessary requirements
./install_dependencies.sh

# 4th step: deactivate environment
deactivate
