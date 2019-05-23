#!/bin/bash
[ -z "$SCRAPP_SIM_CURDIR" ] && echo "SCRAPP_SIM_CURDIR empty! Aborting" && exit

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)

REF=${SCRAPP_SIM_CURDIR}/msa/reference.fasta
QRY=${SCRAPP_SIM_CURDIR}/msa/query.fasta
TREE=${SCRAPP_SIM_CURDIR}/tree/reference.newick
# TREE=${BASE}/tree/eval.raxml.bestTree
MODEL=${SCRAPP_SIM_CURDIR}/tree/eval.raxml.bestModel
# MODEL="GTR+G"

OUT=${SCRAPP_SIM_CURDIR}/placed

mkdir -p ${OUT}
rm ${OUT}/* 2> /dev/null

cd $OUT
epa-ng --tree ${TREE} --msa ${REF} --query ${QRY} --model ${MODEL} --filter-acc-lwr 0.99 --filter-max 10 "$@"
cd -
