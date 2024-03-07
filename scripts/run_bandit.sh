#!/bin/bash -e

# Runs bandit on the Qdislib.
bandit -c ../pyproject.toml -r ..
