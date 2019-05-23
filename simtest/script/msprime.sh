#!/bin/bash

[ -z "$SCRAPP_SIM_CURDIR" ] && echo "SCRAPP_SIM_CURDIR empty! Aborting" && exit

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

OUT=${SCRAPP_SIM_CURDIR}/tree

mkdir -p ${OUT}
rm ${OUT}/* 2> /dev/null

${BASE}/script/trees.py --out-dir ${OUT} "$@"
