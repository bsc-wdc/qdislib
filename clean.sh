#!/usr/bin/env bash

#---------------------------------------------------
# SCRIPT GLOBAL CONSTANTS
#---------------------------------------------------

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#---------------------------------------------------
# REMOVE PACKAGE BUILDING FILES
#---------------------------------------------------

rm -rf "${SCRIPT_DIR}/build"
rm -rf "${SCRIPT_DIR}/Qdislib.egg-info"
rm -rf "${SCRIPT_DIR}/dist"

find "${SCRIPT_DIR}" | grep __pycache__ | xargs rm -rf

#---------------------------------------------------
# REMOVE DOCUMENTATION BUILDING FILES
#---------------------------------------------------

cd docs
./clean.sh
cd ..

#---------------------------------------------------
# REMOVE OTHER FILES
#---------------------------------------------------

rm -f "${SCRIPT_DIR}/output.log"
rm -f "${SCRIPT_DIR}/.coverage"
rm -f "${SCRIPT_DIR}/scripts/output.log"
