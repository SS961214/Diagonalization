#!/bin/sh

if [ -n "$MAGMA_PATH" ]; then exit 0;

if [ $(ldconfig -p | grep libmagma | wc -l) -gt 0 ]; then
  cp basic.mk Makefile
else
  if [ $(locate libmagma | wc -l) -gt 0 ]; then
    MAGMA=$(locate libmagma | head -n 1)
    MAGMA=${MAGMA%%/libmagma*}
    echo "MAGMA_PATH=-L${MAGMA}" > Makefile
    echo basic.mk >> Makefile
  else
    COLOR_f=$$(printf '\033[36m');
    COLOR_b=$$(printf '\033[m');
    echo "$${COLOR_f}Can't find MAGMA library. Please specify its location by setting MAGMA_PATH.$${COLOR_b}"
  fi
fi
