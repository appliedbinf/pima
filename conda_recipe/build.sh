#!/usr/bin/env bash

$PYTHON -m pip install -vv --no-deps --ignore-installed .

#ln -s ${PREFIX}/bin/pima ${PREFIX}/bin/pima.py

# Copy the [de]activate scripts to $PREFIX/etc/conda/[de]activate.d.
# This will allow them to be run on environment activation.
#CHANGE in "activate" "deactivate"; do
#mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
#cp "${RECIPE_DIR}/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
#
