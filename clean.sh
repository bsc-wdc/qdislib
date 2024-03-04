#!/usr/bin/env bash

#---------------------------------------------------
# SCRIPT GLOBAL CONSTANTS
#---------------------------------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

rm -rf "${SCRIPT_DIR}/build"
rm -rf "${SCRIPT_DIR}/src/qdislib.egg-info"
rm -rf "${SCRIPT_DIR}/dist"

find "${SCRIPT_DIR}" | grep __pycache__ | xargs rm -rf
