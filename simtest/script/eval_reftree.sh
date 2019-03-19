#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick

OUT=${BASE}/tree

echo "start at `date`"

cd $OUT
raxml-ng --eval --tree ${TREE} --msa ${REF} --model GTR+G+FO --prefix ${OUT}/eval --threads 10 "$@"
cd -

echo "end at `date`"
