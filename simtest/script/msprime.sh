#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

OUT=${BASE}/tree

mkdir -p ${OUT}
rm ${OUT}/*

${BASE}/script/trees.py --species 30 --out-dir ${OUT} "$@"
