#!/bin/bash
[ -z "$SCRAPP_SIM_CURDIR" ] && echo "SCRAPP_SIM_CURDIR empty! Aborting" && exit

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

REF=${SCRAPP_SIM_CURDIR}/msa/reference.fasta
QRY=${SCRAPP_SIM_CURDIR}/msa/query.fasta
TREE=${SCRAPP_SIM_CURDIR}/tree/reference.newick
MODEL=${SCRAPP_SIM_CURDIR}/tree/eval.raxml.bestModel
JPLACE=${SCRAPP_SIM_CURDIR}/placed/epa_result.jplace

OUT=${SCRAPP_SIM_CURDIR}/delimit

mkdir -p ${OUT}
rm -r ${OUT}/* 2> /dev/null

${BASE}/../scrapp.py  --jplace ${JPLACE} --alignment ${QRY} --work-dir ${OUT} --parallel threads --min-weight 0.5 "$@"
