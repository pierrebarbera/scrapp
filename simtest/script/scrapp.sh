#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=${BASE}/tree/eval.raxml.bestModel
JPLACE=${BASE}/placed/epa_result.jplace

OUT=${BASE}/delimit

mkdir -p ${OUT}
rm -r ${OUT}/*

${BASE}/../scrapp.py  --jplace ${JPLACE} --alignment ${QRY} --work-dir ${OUT} --parallel threads --min-weight 0.5 --num-threads 4 "$@"
