#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
REF=${BASE}/msa/reference.fasta
QRY=${BASE}/msa/query.fasta
TREE=${BASE}/tree/reference.newick
MODEL=${BASE}/tree/eval.raxml.bestModel
JPLACE=${BASE}/placed/epa_result.jplace
SC=${BASE}/script

NUM_THREADS=4

set -e

echo "start at `date`"

cd ${SC}

echo "generate the tree..."
./msprime.sh
echo "tree done!"

# generate the sequences and split into query and ref set
echo "generate the sequences..."
./seqgen.sh
echo "sequences done!"

# infer model params
echo "infer model params..."
./eval_reftree.sh --threads ${NUM_THREADS}
echo "model params done!"

# run placement
echo "place..."
./epa.sh --threads ${NUM_THREADS}
echo "placement done!"

# run scrapp
echo "running scrapp..."
./scrapp.sh --num-threads ${NUM_THREADS}
# ./scrapp.sh --ref-align-outgrouping ${REF}
echo "scrapp done!"

# print statistic
echo "printing statistic..."
./stat.py
./compare_species_counts ../delimit/summary.newick ../tree/annot_reference.newick

echo "end at `date`"
