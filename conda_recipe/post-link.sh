#!/usr/bin/env bash

pima_install.sh >$CONDA_PREFIX/pima_install.log 2>&1

cat >>"$PREFIX/.messages.txt" <<EOF
  ----------------------------------------------------------------------
  
  ${PKG_NAME} version ${PKG_VERSION}-${PKG_BUILDNUM} has been 
  
  successfully installed!
  
  ---------------------------------------------------------------------- 
EOF

