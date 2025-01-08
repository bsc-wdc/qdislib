#!/usr/bin/env bash

cd ..
# Run the coverage of the Qdislib using the tests in ../tests (sequential)
coverage3 run --source Qdislib tests
# Create the report
coverage3 report
# Report coverage results to the CLI.
coverage3 report -m
# Upload coverage report to codecov.io
##bash <(curl -s https://codecov.io/bash) -t <token_id>
cd -
