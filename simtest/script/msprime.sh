#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
TEMPL=${BASE}/script/template/control.txt

TREE=${BASE}/tree/

OUT=${BASE}/tree

mkdir -p ${OUT}
rm ${OUT}/*

${BASE}/script/trees.py 10 ${OUT}