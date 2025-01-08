#!/usr/bin/env bash

#---------------------------------------------------
# SCRIPT GLOBAL CONSTANTS
#---------------------------------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#---------------------------------------------------
# REMOVE PACKAGE BUILDING FILES
#---------------------------------------------------

rm -rf "${SCRIPT_DIR}/build"
rm -rf "${SCRIPT_DIR}/src/Qdislib.egg-info"
rm -rf "${SCRIPT_DIR}/dist"

find "${SCRIPT_DIR}" | grep __pycache__ | xargs rm -rf

#---------------------------------------------------
# REMOVE DOCUMENTATION BUILDING FILES
#---------------------------------------------------

cd docs
make clean
cd ..

#---------------------------------------------------
# REMOVE OTHER FILES
#---------------------------------------------------

rm -f .coverage
rm -f scripts/output.log
rm -rf Qdislib.egg-info
