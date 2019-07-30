#!/bin/bash
BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
[ -z "$SCRAPP_SIM_CURDIR" ] && echo "SCRAPP_SIM_CURDIR empty! Aborting" && exit

TREE=${SCRAPP_SIM_CURDIR}/tree/true_tree.newick

OUT=${SCRAPP_SIM_CURDIR}/msa

mkdir -p ${OUT}
rm ${OUT}/* 2> /dev/null

cd ${OUT}

seq-gen -q -m GTR -g 4 -a 1.0 -or "$@" < ${TREE} > full_TRUE.phy

${BASE}/script/split_msa.py ${SCRAPP_SIM_CURDIR}
