#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
TEMPL=${BASE}/script/template/control.txt

TREE=${BASE}/tree/true_tree.newick

OUT=${BASE}/msa

mkdir -p ${OUT}
rm ${OUT}/*

cd ${OUT}

seq-gen -q -l 1000 -m GTR -g 4 -a 1.0 -or "$@" < ${TREE} > full_TRUE.phy

${BASE}/script/split_msa.py
