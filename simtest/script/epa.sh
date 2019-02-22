#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=GTR+G+FO

OUT=${BASE}/placed

mkdir -p $OUT
rm $OUT/*

echo "start at `date`"

cd $OUT
epa-ng --tree ${TREE} --msa ${REF} --query ${QRY} --model ${MODEL} --filter-acc-lwr 0.99 --filter-max 10 "$@"
cd -

echo "end at `date`"
