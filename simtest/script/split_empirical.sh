#!/bin/bash

[ -z "$SCRAPP_SIM_CURDIR" ] && echo "SCRAPP_SIM_CURDIR empty! Aborting" && exit

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
TREE=${BASE}/real_data/tree.newick
MSA=${BASE}/real_data/tax_cons_border.fasta

OUT=${SCRAPP_SIM_CURDIR}/tree

mkdir -p ${OUT}
rm ${OUT}/* 2> /dev/null

${BASE}/script/real_tree_prune.py --out-dir ${OUT} --msa ${MSA} --tree ${TREE} "$@"

# move result to appropriate folder
mkdir -p ${SCRAPP_SIM_CURDIR}/msa
REF=${SCRAPP_SIM_CURDIR}/msa/reference.fasta
QRY=${SCRAPP_SIM_CURDIR}/msa/query.fasta
mv ${OUT}/reference.fasta ${REF}
mv ${OUT}/query.fasta ${QRY}
