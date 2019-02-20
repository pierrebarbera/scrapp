#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
TEMPL=${BASE}/script/template/control.txt

TREE=${BASE}/tree/tree.newick

OUT=${BASE}/msa

mkdir -p ${OUT}
rm ${OUT}/*

cd ${OUT}

awk 'NR==FNR { a[n++]=$0; next } /%%TREE%%/ {gsub("%%TREE%%",a[0])}1' ${TREE} ${TEMPL} > control.txt

indelible
