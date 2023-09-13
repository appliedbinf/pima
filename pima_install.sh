#!/usr/bin/env bash

python -m pip install --no-deps dna_features_viewer si-prefix 

# upgrade python_circos from the conda version that controls the dependencies, but lacks features
python -m pip install --force-reinstall --no-deps git+https://github.com/ponnhide/pyCircos.git