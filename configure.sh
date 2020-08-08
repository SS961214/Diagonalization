#!/bin/sh

echo "Searching for MAGMA library..."
if [ -n "$MAGMA_PATH" ]; then exit 0;fi

echo 'Executing "ldconfig -p | grep libmagma"...'
if [ $(ldconfig -p | grep libmagma | wc -l) -gt 0 ]; then
    cp basic.mk Makefile
    echo "Makefile was successfully created."
    exit 0
fi

echo 'Executing "locate libmagma"...'
if [ $(locate libmagma | wc -l) -gt 0 ]; then
    MAGMA=$(locate libmagma | head -n 1)
    MAGMA=${MAGMA%%/libmagma*}
    echo "MAGMA_PATH=-L${MAGMA}" > Makefile
    cat basic.mk >> Makefile
    echo "Makefile was successfully created."
else
    COLOR_f=$(printf '\033[31m')
    COLOR_b=$(printf '\033[m')
    echo "$${COLOR_f}Can't find MAGMA library. Please specify its location by setting MAGMA_PATH.$${COLOR_b}"
fi
