#!/bin/bash -e

# Runs bandit on the qdislib.
bandit -c ../pyproject.toml -r ..
