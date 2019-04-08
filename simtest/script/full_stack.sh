#!/bin/bash
DATE=`date '+%Y-%m-%d-%H:%M'`

BASE=$(cd `dirname "${BASH_SOURCE[0]}"`/.. && pwd)
SC=${BASE}/script

set -e

echo "start at `date`"

cd $SC

echo "generate the tree..."
./msprime.sh
echo "tree done!"

# generate the sequences and split into query and ref set
echo "generate the sequences..."
./seqgen.sh
echo "sequences done!"

# infer model params
echo "infer model params..."
./eval_reftree.sh
echo "model params done!"

# run placement
echo "place..."
./epa.sh
echo "placement done!"

# run scrapp
echo "running scrapp..."
./scrapp.sh
echo "scrapp done!"

# print statistic
echo "printing statistic..."
./stat.py
./compare_species_counts ../delimit/summary.newick ../tree/annot_reference.newick

echo "end at `date`"
