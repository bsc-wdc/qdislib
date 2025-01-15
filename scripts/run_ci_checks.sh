#!/bin/bash -e

# Test, coverage, and style should be run from the repo's root
export CURRENT_PATH="$(dirname "$(readlink -f "$0")")"
cd ${CURRENT_PATH}

# Add dislib to the python path
export PYTHONPATH=$PYTHONPATH:${CURRENT_PATH}/..

echo "Running Black style check"
./run_style.sh

echo "Running tests"
# Run the tests in ./tests with PyCOMPSs
./run_tests.sh

echo "Running code coverage"
./run_coverage.sh
