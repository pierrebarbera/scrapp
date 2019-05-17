#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

OUT=${BASE}/tree

mkdir -p ${OUT}
rm ${OUT}/*

${BASE}/script/trees.py 30 ${OUT} "$@"
